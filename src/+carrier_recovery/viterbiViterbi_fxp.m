function [v, ThetaPU] = viterbiViterbi_fxp(x, NPol, NTaps, VVFilter, ...
                                                Pilots, BlockLen, StepSize, ...
                                                PilotThreshold, CordicIts, T) %#codegen
%vitERBIVITERBI_FXP  Fixed-point Viterbi-Viterbi carrier phase recovery
%                        with step-based phase update and optional pilot-aided
%                        cycle-slip correction.
%
%   [v, ThetaPU] = viterbiViterbi_fxp(x, NPol, NTaps, VVFilter, ...
%                      Pilots, BlockLen, StepSize, PilotThreshold, CordicIts, T)
%
%   Inputs
%     x         - input signal [N x NPol] (fi or double, complex)
%     NPol      - number of polarisations (double scalar)
%     NTaps     - one-sided VV filter half-length; window = 2*NTaps+1 (double)
%     VVFilter  - VV filter coefficients [(2*NTaps+1) x 1] (fi or double, real)
%     Pilots    - pilot symbols, one per block [NBlocks x NPol] (complex fi or double)
%     BlockLen  - block length in symbols (double scalar)
%     StepSize  - phase update interval in symbols (double scalar, 1..BlockLen)
%                 The estimator, unwrapper and pilot correction fire once every
%                 StepSize symbols, aligned to the start of each block.
%                 The phase is held constant between updates.
%                   StepSize = 1        -> symbol-by-symbol (full bandwidth)
%                   StepSize = BlockLen -> one update per block (minimum bandwidth)
%                 Pilot symbols are treated as regular data by the VV estimator.
%     PilotThreshold - threshold for pilot-based cycle-slip correction in radians (double scalar)
%     CordicIts - number of iterations for CORDIC operations (double scalar)
%                 Defaults to 'fixed16'.
%
%   Outputs
%     v       - phase-corrected signal [N x NPol], type T.x
%     ThetaPU - phase estimate [N x NPol], type T.theta
%               At step positions: newly estimated, unwrapped, pilot-corrected.
%               Between steps: held from the most recent step.
%
%   Fixed-point types table T must supply:
%     T.x     - input / output signal type
%     T.w     - filter coefficient type
%     T.theta - phase / angle type  (must accommodate ±pi)
%     T.acc   - accumulator type    (pilot correlation sums)
%
%   Codegen notes
%     - No convmtx: tap-delay indexing throughout.
%     - ThetaML is not pre-allocated; estimation and unwrapping are merged
%       into a single loop and skipped entirely at non-step positions,
%       reducing computation by a factor of StepSize.
%     - cordicangle and cordicrotate outputs are explicitly cast to the
%       intended fi type immediately after each call (CORDIC ignores fimath).
%     - ThetaPrev is updated only at step positions; the unwrapper anchor
%       therefore always reflects the last computed (not held) phase.
%     - UsePilots is a runtime branch; codegen compiles both paths.

    %% ----------------------------------------------------------------
    %  Default types table
    %% ----------------------------------------------------------------
    if nargin < 10 || isempty(T)
        T = carrier_recovery.fxp_types('fixed16');
    end

    %% ----------------------------------------------------------------
    %  Fixed-point constants
    %% ----------------------------------------------------------------
    PI_OVER2 = cast(pi/2, 'like', T.theta);
    PI_VAL   = cast(pi,   'like', T.theta);
    PI_OVER4 = cast(pi/4, 'like', T.theta);
    ZERO_TH  = cast(0, 'like', T.theta);
    QUARTER  = cast(0.25, 'like', T.theta);
    ZERO_ACC = cast(0,    'like', T.acc);
    CORDIC_ITS = coder.const(CordicIts);

    %% ----------------------------------------------------------------
    %  Dimensions
    %% ----------------------------------------------------------------
    N       = size(x, 1);
    L_filt  = 2 * NTaps + 1;
    halfL   = floor(L_filt / 2);
    NBlocks = ceil(N / BlockLen);

    %% ----------------------------------------------------------------
    %  Cast inputs to fixed-point
    %% ----------------------------------------------------------------
    x_fi = cast(x,        'like', T.x);
    w_fi = cast(VVFilter, 'like', T.w);

    %% ----------------------------------------------------------------
    %  Pre-allocate outputs
    %% ----------------------------------------------------------------
    ThetaPU = zeros(N, NPol, 'like', T.theta);
    v       = complex(zeros(N, NPol, 'like', T.x));

    %% ================================================================
    %  Pilot correlation  -->  PhiRef [NBlocks x NPol]
    %
    %  All pilot references are computed upfront before the main loop
    %  so that any step position can look up its block's reference
    %  without ordering constraints.
    %% ================================================================
    PhiRef    = zeros(NBlocks, NPol, 'like', T.theta);
    Pilots_fi = cast(Pilots, 'like', T.x);

    for pol = 1:NPol
        for blk = 1:min(NBlocks, size(Pilots, 1))
            blockStart = (blk - 1) * BlockLen + 1;
            if blockStart <= N
                rx = x_fi(blockStart, pol);

                pilot_re =  cast(real(Pilots_fi(blk, pol)), 'like', T.acc);
                pilot_im = -cast(imag(Pilots_fi(blk, pol)), 'like', T.acc);
                rx_re    =  cast(real(rx), 'like', T.acc);
                rx_im    =  cast(imag(rx), 'like', T.acc);

                corr_re = pilot_re * rx_re - pilot_im * rx_im;
                corr_im = pilot_re * rx_im + pilot_im * rx_re;

                % cordicangle ignores fimath and returns FL = (input FL - 2).
                % Cast immediately to T.theta to restore the correct
                % numerictype and SpecifyPrecision fimath.
                corr_fi          = complex(cast(corr_re, 'like', T.theta), ...
                                           cast(corr_im, 'like', T.theta));
                PhiRef(blk, pol) = cast(cordicangle(corr_fi, CORDIC_ITS), 'like', T.theta);
            end
        end
    end

    %% ================================================================
    %  Per-polarisation processing
    %
    %  Estimation, unwrapping, pilot correction, and phase hold are
    %  merged into a single pass.  At step positions the full pipeline
    %  executes; at all other positions the previous phase is replicated.
    %
    %  Step positions within each block are determined by:
    %    posInBlock = mod(i-1, BlockLen)   (0-indexed, resets each block)
    %    isStep     = mod(posInBlock, StepSize) == 0
    %  The first symbol of every block (posInBlock = 0) is always a step.
    %% ================================================================
    for pol = 1:NPol

        ThetaPrev = ZERO_TH;    % unwrapper anchor; updated only at steps

        for i = 1:N

            posInBlock = mod(i - 1, BlockLen);
            isStep     = (mod(posInBlock, StepSize) == 0);

            if isStep

                %% ------------------------------------------------
                %  VV phase estimation
                %  theta_ml = angle(sum_k w(k)*x(i+k-halfL)^4)/4 - pi/4
                %% ------------------------------------------------
                sum4_re = ZERO_TH;
                sum4_im = ZERO_TH;

                for k = 1:L_filt
                    idx = i - halfL - 1 + k;
                    if idx >= 1 && idx <= N
                        s = x_fi(idx, pol);
                    else
                        s = complex(cast(0, 'like', T.x), cast(0, 'like', T.x));
                    end

                    % s^4 via two complex squarings (avoids fi .^4)
                    s_re = cast(real(s), 'like', T.theta);
                    s_im = cast(imag(s), 'like', T.theta);

                    s2_re = s_re * s_re - s_im * s_im;
                    s2_im = cast(2, 'like', T.theta) * s_re * s_im;

                    s4_re = s2_re * s2_re - s2_im * s2_im;
                    s4_im = cast(2, 'like', T.theta) * s2_re * s2_im;

                    w_k     = cast(w_fi(k), 'like', T.theta);
                    sum4_re = sum4_re + w_k * s4_re;
                    sum4_im = sum4_im + w_k * s4_im;
                end

                sum4_fi  = complex(sum4_re, sum4_im);
                % cordicangle output FL = (input FL - 2); cast immediately.
                theta_ml = cast(cordicangle(sum4_fi, CORDIC_ITS), 'like', T.theta) ...
                           * QUARTER - PI_OVER4;

                %% ------------------------------------------------
                %  Phase unwrapping
                %  n = floor(0.5 + (ThetaPrev - theta_ml) / (pi/2))
                %% ------------------------------------------------
                diff_val = ThetaPrev - theta_ml;
                n_val    = floor(double(diff_val) / double(PI_OVER2) + 0.5);
                n_fi     = cast(n_val, 'like', T.theta);
                theta_uw = theta_ml + n_fi * PI_OVER2;

                % theta_ml in (-pi/4, pi/4], n_fi in {-1,0,1}:
                % theta_uw bounded to (-3pi/4, 3pi/4] — no fi overflow risk.

                %% ------------------------------------------------
                %  Pilot-aided cycle-slip correction
                %% ------------------------------------------------
                BlockIdx = ceil(i / BlockLen);

                % Wrap difference to (-pi, pi] for shortest-path error,
                % then threshold at ±pi/2 to identify a one-quadrant slip.
                PhaseDiff_d = mod( ...
                    double(theta_uw - PhiRef(BlockIdx, pol)) + pi, ...
                    2*pi) - pi;
                n_slip = round(PhaseDiff_d / double(PilotThreshold));
                n_slip_fi = cast(n_slip, 'like', T.theta);
                theta_uw = theta_uw - n_slip_fi * PI_OVER2;

                ThetaPU(i, pol) = theta_uw;
                ThetaPrev       = theta_uw;   % advance anchor to this step

            else

                %% ------------------------------------------------
                %  Hold: replicate the last computed phase
                %  ThetaPrev is not updated so the next step's
                %  unwrapper still anchors to the last real estimate.
                %% ------------------------------------------------
                ThetaPU(i, pol) = ThetaPrev;

            end

        end  % for i

        %% ------------------------------------------------------------
        %  Phase correction: v(i) = x(i) * exp(-j * ThetaPU(i))
        %  cordicrotate output cast to T.x to enforce SpecifyPrecision.
        %% ------------------------------------------------------------
        for i = 1:N
            % CORDIC rotation is most reliable in the principal range.
            % Reduce angle to [-pi, pi], then map to [-pi/2, pi/2]
            % using a sign flip of the input symbol for quadrant handling.
            theta_d = mod(double(-ThetaPU(i, pol)) + pi, 2*pi) - pi;
            s_in    = x_fi(i, pol);

            if theta_d > pi/2
                theta_d = theta_d - pi;
                s_in    = -s_in;
            elseif theta_d < -pi/2
                theta_d = theta_d + pi;
                s_in    = -s_in;
            end

            theta_safe = cast(theta_d, 'like', T.theta);
            if theta_safe > PI_OVER2
                theta_safe = theta_safe - PI_VAL;
                s_in       = -s_in;
            elseif theta_safe < -PI_OVER2
                theta_safe = theta_safe + PI_VAL;
                s_in       = -s_in;
            end

            v(i, pol) = cast(cordicrotate(theta_safe, s_in, CORDIC_ITS), 'like', T.x);
        end

    end  % for pol
end