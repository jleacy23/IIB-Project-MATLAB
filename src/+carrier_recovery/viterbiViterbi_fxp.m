function [v, ThetaPU] = viterbiViterbi_fxp(x, NPol, NTaps, VVFilter, ...
                                                Pilots, BlockLen, StepSize, ...
                                                PilotThreshold, CordicIts, T) %#codegen
%vitERBIVITERBI_FXP  Fixed-point Viterbi-Viterbi carrier phase recovery
%                        with block-based phase update and optional pilot-aided
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
%     StepSize  - retained for interface compatibility (double scalar, 1..BlockLen)
%                 In this implementation, one phase estimate is produced per block
%                 and applied to all symbols in that block (equivalent to
%                 StepSize = BlockLen behavior).
%                 Pilot symbols are treated as regular data by the VV estimator.
%     PilotThreshold - threshold for pilot-based cycle-slip correction in radians (double scalar)
%     CordicIts - number of iterations for CORDIC operations (double scalar)
%                 Defaults to 'fixed16'.
%
%   Outputs
%     v       - phase-corrected signal [N x NPol], type T.x
%     ThetaPU - phase estimate [N x NPol], type T.theta
%               One estimated, unwrapped, pilot-corrected phase per block,
%               replicated across all symbols in that block.
%
%   Fixed-point types table T must supply:
%     T.x     - input / output signal type
%     T.w     - filter coefficient type
%     T.theta - phase / angle type  (must accommodate ±pi)
%     T.acc   - accumulator type    (pilot correlation sums)
%
%   Codegen notes
%     - No convmtx: tap-delay indexing throughout.
%     - One VV ML phase estimate is computed per block from all symbols in
%       that block, then unwrapped and optionally pilot-corrected.
%     - cordicangle and cordicrotate outputs are explicitly cast to the
%       intended fi type immediately after each call (CORDIC ignores fimath).
%     - ThetaPrev is updated once per block; the unwrapper anchor therefore
%       always reflects the last computed block phase.
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

    % Retained for API compatibility; block mode always performs one
    % estimator/unwrap/correction update per block.
    if StepSize < 1 || StepSize > BlockLen
        error('carrier_recovery:viterbiViterbi_fxp:InvalidStepSize', ...
              'StepSize must be in [1, BlockLen].');
    end

    %% ================================================================
    %  Per-polarisation processing
    %
    %  For each block, accumulate VV ML statistics across all symbols
    %  in the block, perform one unwrapping/correction update, then
    %  replicate that phase over the entire block.
    %% ================================================================
    for pol = 1:NPol

        ThetaPrev = ZERO_TH;    % unwrapper anchor; updated once per block

        for blk = 1:NBlocks
            blockStart = (blk - 1) * BlockLen + 1;
            blockEnd   = min(blk * BlockLen, N);

            %% ------------------------------------------------
            %  Block VV phase estimation
            %  theta_ml = angle(sum_i sum_k w(k)*x(i+k-halfL)^4)/4 - pi/4
            %% ------------------------------------------------
            sum4_re_blk = ZERO_TH;
            sum4_im_blk = ZERO_TH;

            for i = blockStart:blockEnd
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

                    w_k = cast(w_fi(k), 'like', T.theta);
                    sum4_re_blk = sum4_re_blk + w_k * s4_re;
                    sum4_im_blk = sum4_im_blk + w_k * s4_im;
                end
            end

            sum4_fi  = complex(sum4_re_blk, sum4_im_blk);
            % cordicangle output FL = (input FL - 2); cast immediately.
            theta_ml = cast(cordicangle(sum4_fi, CORDIC_ITS), 'like', T.theta) ...
                       * QUARTER - PI_OVER4;

            %% ------------------------------------------------
            %  Phase unwrapping (same anchor/update rule as one step per block)
            %  n = floor(0.5 + (ThetaPrev - theta_ml) / (pi/2))
            %% ------------------------------------------------
            diff_val = ThetaPrev - theta_ml;
            n_val    = floor(double(diff_val) / double(PI_OVER2) + 0.5);
            n_fi     = cast(n_val, 'like', T.theta);
            theta_uw = theta_ml + n_fi * PI_OVER2;

            %% ------------------------------------------------
            %  Pilot-aided cycle-slip correction (per block, unchanged)
            %% ------------------------------------------------
            PhaseDiff_d = mod( ...
                double(theta_uw - PhiRef(blk, pol)) + pi, ...
                2*pi) - pi;
            n_slip = round(PhaseDiff_d / double(PilotThreshold));
            n_slip_fi = cast(n_slip, 'like', T.theta);
            theta_uw = theta_uw - n_slip_fi * PI_OVER2;

            ThetaPU(blockStart:blockEnd, pol) = theta_uw;
            ThetaPrev = theta_uw;
        end  % for blk

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