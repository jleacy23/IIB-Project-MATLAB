function [v, ThetaPU] = bps_fxp(z, N, NPol, M, B, BlockLen, StepSize, ...
                                     Pilots, PilotThreshold, CordicIts, T) %#codegen
%bps_FXP  Fixed-point Blind Phase Search (BPS) carrier phase recovery
%             with step-based phase update and optional pilot-aided
%             cycle-slip correction.
%
%   [v, ThetaPU] = bps_fxp(z, N, NPol, M, B, BlockLen, StepSize,
%                               Pilots, PilotThreshold, CordicIts, T)
%
%   Inputs
%     z         - input signal [Nsym x NPol] (fi or double, complex)
%     N         - one-sided BPS filter half-length; window = 2*N+1 (double)
%     NPol      - number of polarisations (double)
%     M         - QAM order, e.g. 4, 16, 64 (double)
%     B         - number of blind test phases (double, must be even)
%     BlockLen  - block length in symbols (double scalar)
%     StepSize  - phase update interval in symbols (double scalar, 1..BlockLen)
%                 The BPS estimator, unwrapper and pilot correction fire once
%                 every StepSize symbols, aligned to the start of each block.
%                 The phase is held constant between updates.
%                   StepSize = 1        -> symbol-by-symbol (full bandwidth)
%                   StepSize = BlockLen -> one update per block (minimum bandwidth)
%                 Pilot symbols are treated as regular data by the BPS estimator.
%     Pilots    - pilot symbols, one per block [NBlocks x NPol] (complex fi or double)
%     PilotThreshold - threshold for pilot-based cycle-slip correction in radians (double scalar)
%     CordicIts - number of CORDIC iterations for the pilot-reference angle,
%                 every blind test-phase rotation, and the final per-symbol
%                 de-rotation (double scalar).  Angular resolution is
%                 ~atan(2^-CordicIts); equals the swept precision in the sweep.
%
%   Outputs
%     v       - phase-corrected signal [Nsym x NPol], type T.x
%     ThetaPU - phase estimate [Nsym x NPol], type T.theta
%               At step positions: newly estimated, unwrapped, pilot-corrected.
%               Between steps: held from the most recent step.
%
%   Fixed-point implementation notes
%     - convmtx is replaced by explicit tap-delay window indexing.
%     - Test-phase rotations use CORDIC rotation (cordic.rotate).
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
%   Fixed-point cast notes
%     All angles and rotations go through CORDIC (cordic.vectoring /
%     cordic.rotate): the pilot-reference angle, each blind test-phase
%     rotation in the BPS metric, and the final per-symbol de-rotation.
%     The BPS metric itself (squared slicer distance) stays in double, as it
%     accumulates across taps and would otherwise need a very wide fi type;
%     each candidate's CORDIC-rotated sample is cast to double only to feed
%     modem.slicer.

    %% ----------------------------------------------------------------
    %  Default types table
    %% ----------------------------------------------------------------
    if nargin < 11 || isempty(T)
        T = carrier_recovery.fxp_types('fixed16');
    end

    %% ----------------------------------------------------------------
    %  Fixed-point constants
    %% ----------------------------------------------------------------
    PI_OVER2 = cast(pi/2, 'like', T.theta);
    ZERO_ACC = cast(0,    'like', T.acc);

    %% ----------------------------------------------------------------
    %  Dimensions
    %% ----------------------------------------------------------------
    Nsym    = size(z, 1);
    L       = 2 * N + 1;       % BPS averaging window length
    halfL   = N;
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

    for blk = 1:min(NBlocks, size(Pilots, 1))
        blockStart = (blk - 1) * BlockLen + 1;
        if blockStart <= Nsym
            for pol = 1:NPol
                rx = z_fi(blockStart, pol);

                pilot_re =  cast(real(Pilots_fi(blk, pol)), 'like', T.acc);
                pilot_im = -cast(imag(Pilots_fi(blk, pol)), 'like', T.acc);
                rx_re    =  cast(real(rx), 'like', T.acc);
                rx_im    =  cast(imag(rx), 'like', T.acc);

                corr_re = pilot_re * rx_re - pilot_im * rx_im;
                corr_im = pilot_re * rx_im + pilot_im * rx_re;

                PhiRef(blk, pol) = cordic.vectoring(corr_re, corr_im, CordicIts, T);
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
                        [sr, si] = cordic.rotate(real(s), imag(s), ...
                                                 -ThetaTest_fi(b), CordicIts, T);
                        s_rot_d  = complex(double(sr), double(si));
                        s_dec    = modem.slicer(s_rot_d);

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
                BlockIdx = ceil(i / BlockLen);

                % Wrap difference to (-pi, pi] for shortest-path error,
                % then threshold at ±pi/2 to identify a one-quadrant slip.
                PhaseDiff_d = mod( ...
                    double(theta_uw - PhiRef(BlockIdx, pol)) + pi, ...
                    2*pi) - pi;
                n_slip = round(PhaseDiff_d / PilotThreshold);
                n_slip_fi = cast(n_slip, 'like', T.theta);
                theta_uw = theta_uw - n_slip_fi * PI_OVER2;

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
    %  Phase correction: v(i,pol) = z(i,pol) * exp(-j * ThetaPU(i,pol)) via CORDIC
    %% ================================================================
    for i = 1:Nsym
        for pol = 1:NPol
            [vr, vi] = cordic.rotate(real(z_fi(i, pol)), imag(z_fi(i, pol)), ...
                                     -ThetaPU(i, pol), CordicIts, T);
            v(i, pol) = complex(vr, vi);
        end
    end
end