function [v, ThetaPU] = bps_fxp(z, N, NPol, M, B, BlockLen, StepSize, ...
                                     Pilots, UsePilots, PilotThreshold, CordicIts, T) %#codegen
%bps_FXP  Fixed-point Blind Phase Search (BPS) carrier phase recovery
%             with step-based phase update and optional pilot-aided
%             cycle-slip correction.
%
%   [v, ThetaPU] = bps_fxp(z, N, NPol, M, B, BlockLen, StepSize,
%                               Pilots, UsePilots, PilotThreshold, CordicIts, T)
%
%   Inputs
%     z         - input signal [Nsym x NPol] (fi or double, complex)
%     N         - one-sided BPS filter half-length; window = 2*N+1 (double)
%     NPol      - number of polarisations (double)
%     M         - QAM order, e.g. 4, 16, 64 (double)
%     B         - number of blind test phases (double, must be even)
%     BlockLen  - block length in symbols (double scalar)
%                 Pilots are taken from the first P symbols of each block.
%     StepSize  - phase update interval in symbols (double scalar, 1..BlockLen)
%                 The BPS estimator, unwrapper and pilot correction fire once
%                 every StepSize symbols, aligned to the start of each block.
%                 The phase is held constant between updates.
%                   StepSize = 1        -> symbol-by-symbol (full bandwidth)
%                   StepSize = BlockLen -> one update per block (minimum bandwidth)
%                 Pilot symbols are treated as regular data by the BPS estimator.
%     Pilots    - pilot symbols at block start [P x 1] (complex fi or double)
%                 Ignored when UsePilots = false
%     UsePilots - logical: enable pilot-aided cycle-slip correction
%     PilotThreshold - threshold for pilot-based cycle-slip correction in radians (double scalar)
%     CordicIts - number of iterations for CORDIC operations (double scalar)
%                 Defaults to 'fixed16'.
%
%   Outputs
%     v       - phase-corrected signal [Nsym x NPol], type T.x
%     ThetaPU - phase estimate [Nsym x NPol], type T.theta
%               At step positions: newly estimated, unwrapped, pilot-corrected.
%               Between steps: held from the most recent step.
%
%   Fixed-point implementation notes
%     - convmtx is replaced by explicit tap-delay window indexing.
%     - Test-phase rotations use cordicrotate, matching hardware CORDIC datapaths.
%     - The BPS metric accumulator (m_buf) is a plain double vector: it sums
%       squared distances which grow with L and SNR, making fi overflow likely.
%       Decisions (modem.slicer) are inherently double; keeping the metric in
%       double is consistent and avoids an unnecessarily wide fi accumulator.
%     - Thetas is not pre-allocated; estimation and unwrapping are merged into
%       a single loop, and the BPS metric computation is skipped entirely at
%       non-step positions, reducing computation by a factor of StepSize.
%     - ThetaPrev is updated only at step positions; the unwrapper anchor
%       therefore always reflects the last computed (not held) phase.
%
%   CORDIC type behaviour -- explicit casts are mandatory
%     cordicangle and cordicrotate IGNORE the fimath of their fi inputs.
%     Their outputs carry MATLAB's default FullPrecision fimath, not the
%     SpecifyPrecision fimath in T.  Additionally, cordicangle returns a
%     fraction length of (input FL - 2), which does not match T.theta.
%     Every CORDIC output is explicitly cast back to the intended fi type
%     immediately after the call, before any further arithmetic.

    %% ----------------------------------------------------------------
    %  Default types table
    %% ----------------------------------------------------------------
    if nargin < 10 || isempty(T)
        T = bps_fxp_types('fixed16');
    end

    %% ----------------------------------------------------------------
    %  Fixed-point constants
    %% ----------------------------------------------------------------
    PI_OVER2 = cast(pi/2, 'like', T.theta);
    ZERO_ACC = cast(0,    'like', T.acc);
    CORDIC_ITS = coder.const(CordicIts);

    %% ----------------------------------------------------------------
    %  Dimensions
    %% ----------------------------------------------------------------
    Nsym    = size(z, 1);
    L       = 2 * N + 1;       % BPS averaging window length
    halfL   = N;
    P       = length(Pilots);
    NBlocks = ceil(Nsym / BlockLen);

    %% ----------------------------------------------------------------
    %  Cast inputs to fixed-point
    %% ----------------------------------------------------------------
    z_fi      = cast(z,      'like', T.x);
    Pilots_fi = cast(Pilots, 'like', T.x);

    %% ----------------------------------------------------------------
    %  Precompute test phase angles as fi
    %
    %  Loop index b = 1..B maps to b_val = b-1-B/2 = -B/2..B/2-1,
    %  giving ThetaTest = (pi/2)*b_val/B in (-pi/4, pi/4].
    %% ----------------------------------------------------------------
    ThetaTest_fi = zeros(1, B, 'like', T.theta);

    for b = 1:B
        b_val   = (b - 1) - B/2;
        theta_d = (pi/2) * b_val / B;
        ThetaTest_fi(b) = cast(theta_d, 'like', T.theta);
    end

    %% ================================================================
    %  Pilot correlation  -->  PhiRef [NBlocks x NPol]
    %
    %  All pilot references are computed upfront before the main loop
    %  so that any step position can look up its block's reference.
    %  Correlation is accumulated across ALL polarisations (coherent
    %  combining) before taking angle(), matching the floating-point
    %  reference implementation.
    %% ================================================================
    PhiRef = zeros(NBlocks, NPol, 'like', T.theta);

    if UsePilots
        for blk = 1:NBlocks
            blockStart = (blk - 1) * BlockLen + 1;

            corr_re = ZERO_ACC;
            corr_im = ZERO_ACC;

            for pol = 1:NPol
                for p = 1:P
                    idx = blockStart + p - 1;
                    if idx >= 1 && idx <= Nsym
                        rx = z_fi(idx, pol);

                        pilot_re =  cast(real(Pilots_fi(p)), 'like', T.acc);
                        pilot_im = -cast(imag(Pilots_fi(p)), 'like', T.acc);
                        rx_re    =  cast(real(rx),           'like', T.acc);
                        rx_im    =  cast(imag(rx),           'like', T.acc);

                        corr_re = corr_re + pilot_re * rx_re - pilot_im * rx_im;
                        corr_im = corr_im + pilot_re * rx_im + pilot_im * rx_re;
                    end
                end
            end

            % cordicangle ignores fimath and returns FL = (input FL - 2).
            % Cast immediately to T.theta to restore the correct
            % numerictype and SpecifyPrecision fimath.
            corr_fi = complex(cast(corr_re, 'like', T.theta), ...
                              cast(corr_im, 'like', T.theta));
            phi     = cast(cordicangle(corr_fi, CORDIC_ITS), 'like', T.theta);

            % Same reference broadcast to all pols (coherent combining)
            for pol = 1:NPol
                PhiRef(blk, pol) = phi;
            end
        end
    end

    %% ================================================================
    %  Pre-allocate outputs
    %% ================================================================
    ThetaPU   = zeros(Nsym, NPol, 'like', T.theta);
    v         = complex(zeros(Nsym, NPol, 'like', T.x));
    m_buf     = zeros(1, B);   % double BPS metric buffer, re-used each step

    %% ================================================================
    %  Per-polarisation processing
    %
    %  Estimation, unwrapping, pilot correction, and phase hold are merged
    %  into a single pass.  At step positions the full BPS pipeline executes;
    %  at all other positions the previous phase is replicated.
    %
    %  Step positions within each block are determined by:
    %    posInBlock = mod(i-1, BlockLen)   (0-indexed, resets each block)
    %    isStep     = mod(posInBlock, StepSize) == 0
    %  The first symbol of every block is always a step (posInBlock = 0).
    %% ================================================================
    ThetaPrev = zeros(1, NPol, 'like', T.theta);   % per-pol unwrapper anchors

    for i = 1:Nsym

        posInBlock = mod(i - 1, BlockLen);
        isStep     = (mod(posInBlock, StepSize) == 0);

        for pol = 1:NPol

            if isStep

                %% ------------------------------------------------
                %  BPS phase estimation at this step position
                %% ------------------------------------------------
                for b = 1:B
                    m_buf(b) = 0.0;
                end

                for k = 1:L
                    idx = i - halfL - 1 + k;
                    if idx >= 1 && idx <= Nsym
                        s = z_fi(idx, pol);
                    else
                        s = complex(cast(0, 'like', T.x), cast(0, 'like', T.x));
                    end

                    for b = 1:B
                        % cordicrotate ignores fimath; cast output to T.x.
                        % Angle negated and cast to T.theta before CORDIC call.
                        neg_theta = cast(-ThetaTest_fi(b), 'like', T.theta);
                        s_rot_fi  = cast(cordicrotate(neg_theta, s, CORDIC_ITS), 'like', T.x);

                        s_rot_d  = complex(double(real(s_rot_fi)), ...
                                           double(imag(s_rot_fi)));
                        s_dec    = modem.slicer(s_rot_d, M);

                        err_re   = real(s_rot_d) - real(s_dec);
                        err_im   = imag(s_rot_d) - imag(s_dec);
                        m_buf(b) = m_buf(b) + err_re*err_re + err_im*err_im;
                    end
                end

                [~, best_b]  = min(m_buf);
                theta_est    = ThetaTest_fi(best_b);

                %% ------------------------------------------------
                %  Phase unwrapping
                %  n = floor(0.5 - (theta_est - ThetaPrev) / (pi/2))
                %% ------------------------------------------------
                diff_val = theta_est - ThetaPrev(pol);
                n_val    = floor(0.5 - double(diff_val) / double(PI_OVER2));
                n_fi     = cast(n_val, 'like', T.theta);
                theta_uw = theta_est + n_fi * PI_OVER2;

                % theta_est in (-pi/4, pi/4], n_fi in {-1,0,1}:
                % theta_uw bounded to (-3pi/4, 3pi/4] — no fi overflow risk.

                %% ------------------------------------------------
                %  Pilot-aided cycle-slip correction
                %% ------------------------------------------------
                if UsePilots
                    BlockIdx = ceil(i / BlockLen);

                    % Wrap difference to (-pi, pi] for shortest-path error,
                    % then threshold at ±pi/2 to identify a one-quadrant slip.
                    PhaseDiff_d = mod( ...
                        double(theta_uw - PhiRef(BlockIdx, pol)) + pi, ...
                        2*pi) - pi;
                    n_slip = round(PhaseDiff_d / PilotThreshold);
                    n_slip_fi = cast(n_slip, 'like', T.theta);
                    theta_uw = theta_uw - n_slip_fi * PI_OVER2;
                end

                ThetaPU(i, pol)   = theta_uw;
                ThetaPrev(pol)    = theta_uw;  % advance anchor to this step

            else

                %% ------------------------------------------------
                %  Hold: replicate the last computed phase.
                %  ThetaPrev(pol) not updated so the next step's
                %  unwrapper anchors to the last real estimate.
                %% ------------------------------------------------
                ThetaPU(i, pol) = ThetaPrev(pol);

            end

        end  % for pol
    end  % for i

    %% ================================================================
    %  Phase correction: v(i,pol) = z(i,pol) * exp(-j * ThetaPU(i,pol))
    %  cordicrotate output cast to T.x to enforce SpecifyPrecision fimath.
    %% ================================================================
    for i = 1:Nsym
        for pol = 1:NPol
            neg_theta = cast(-ThetaPU(i, pol), 'like', T.theta);
            v(i, pol) = cast(cordicrotate(neg_theta, z_fi(i, pol), CORDIC_ITS), 'like', T.x);
        end
    end
end