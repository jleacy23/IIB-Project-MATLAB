function [v, ThetaPU] = cr_bps_fxp(z, N, NPol, M, B, BlockLen, ...
                                     Pilots, UsePilots, BlockBased, T) %#codegen
%CR_BPS_FXP  Fixed-point Blind Phase Search (BPS) carrier phase recovery.
%
%   [v, ThetaPU] = cr_bps_fxp(z, N, NPol, M, B, BlockLen,
%                               Pilots, UsePilots, BlockBased, T)
%
%   Inputs
%     z          - input signal [Nsym x NPol] (fi or double, complex)
%     N          - one-sided filter half-length; window = 2*N+1 (double)
%     NPol       - number of polarisations (double)
%     M          - QAM order, e.g. 4, 16, 64 (double)
%     B          - number of blind test phases (double, must be even)
%     BlockLen   - block length in symbols (double)
%     Pilots     - pilot symbols at block start [P x 1] (complex fi or double)
%     UsePilots  - logical: enable pilot-aided cycle-slip correction
%     BlockBased - logical: hold phase constant over each BlockLen-symbol block
%     T          - (optional) fixed-point types table from cr_bps_fxp_types.
%                  Defaults to 'fixed16'.
%
%   Outputs
%     v          - phase-corrected signal [Nsym x NPol], type T.x
%     ThetaPU    - unwrapped phase estimate [Nsym x NPol], type T.theta
%
%   Fixed-point implementation notes
%     - convmtx is replaced by explicit tap-delay window indexing.
%     - Test-phase rotations use cordicrotate(-ThetaTest(b), s), matching
%       the CORDIC rotation datapath used in hardware implementations.
%       Test phase angles are precomputed once as fi scalars; cordicrotate
%       is called once per (sample, test-phase) pair inside the inner loop.
%     - Symbol decisions are made in double precision:
%         (1) rotate sample via cordicrotate  (2) cast to complex double
%         (3) qam_slicer(s_rot, M)           (4) compute squared error in double
%       qam_slicer implements nearest-neighbour QAM decisions using pure
%       arithmetic (scale / round-to-odd / clamp / rescale), with no lookup
%       tables and no dependency on Communications Toolbox.
%
%   CORDIC type behaviour -- explicit casts are mandatory
%     cordicangle and cordicrotate IGNORE the fimath of their fi inputs.
%     Their outputs carry MATLAB's default FullPrecision fimath, not the
%     SpecifyPrecision fimath in T.  Additionally, cordicangle returns a
%     fraction length of (input FL - 2), which does not match T.theta.
%     Every CORDIC output is therefore explicitly cast back to the intended
%     fi type immediately after the call, before any further arithmetic.
%     Without this, subsequent fi operations would silently revert to
%     FullPrecision (allowing unconstrained bit growth) and the type mismatch
%     from cordicangle would propagate undetected.

    %% ----------------------------------------------------------------
    %  Default types table
    %% ----------------------------------------------------------------
    if nargin < 10 || isempty(T)
        T = cr_bps_fxp_types('fixed16');
    end

    %% ----------------------------------------------------------------
    %  Fixed-point constants
    %% ----------------------------------------------------------------
    PI_OVER2 = cast(pi/2, 'like', T.theta);
    ZERO_TH  = cast(0,    'like', T.theta);
    ONE_TH   = cast(1,    'like', T.theta);
    ZERO_ACC = cast(0,    'like', T.acc);

    %% ----------------------------------------------------------------
    %  Dimensions
    %% ----------------------------------------------------------------
    Nsym    = size(z, 1);
    L       = 2 * N + 1;
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
    %  Loop index b = 1..B  maps to  b_val = b - 1 - B/2  = -B/2..B/2-1
    %  giving ThetaTest = (pi/2) * b_val / B  in (-pi/4, pi/4].
    %% ----------------------------------------------------------------
    ThetaTest_fi = zeros(1, B, 'like', T.theta);

    for b = 1:B
        b_val   = (b - 1) - B/2;
        theta_d = (pi/2) * b_val / B;
        ThetaTest_fi(b) = cast(theta_d, 'like', T.theta);
    end

    %% ----------------------------------------------------------------
    %  Pilot correlation  -->  PhiRef [NBlocks x NPol]
    %
    %  Correlation is accumulated across ALL polarisations before taking
    %  angle() (coherent combining), matching the original floating-point
    %  behaviour.
    %% ----------------------------------------------------------------
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

            % cordicangle ignores fimath and returns FL = (input FL - 2),
            % which differs from T.theta's FL.  Cast immediately to bring
            % the result into T.theta's numerictype and SpecifyPrecision
            % fimath before it is stored or used in any further fi arithmetic.
            corr_fi = complex(cast(corr_re, 'like', T.theta), ...
                              cast(corr_im, 'like', T.theta));
            phi = cast(cordicangle(corr_fi), 'like', T.theta);

            for pol = 1:NPol
                PhiRef(blk, pol) = phi;
            end
        end
    end

    %% ================================================================
    %  Step 1: BPS phase estimation
    %% ================================================================
    Thetas = zeros(Nsym, NPol, 'like', T.theta);
    m_buf  = zeros(1, B);          % double metric buffer, re-used each symbol

    for i = 1:Nsym
        for pol = 1:NPol

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
                    % cordicrotate ignores fimath; its output carries default
                    % FullPrecision fimath with T.x's numerictype.  Cast
                    % immediately to T.x to restore SpecifyPrecision fimath
                    % before the double conversion.  The angle is negated and
                    % explicitly cast to T.theta before the CORDIC call so
                    % its type is unambiguous regardless of fi negation rules.
                    neg_theta = cast(-ThetaTest_fi(b), 'like', T.theta);
                    s_rot_fi  = cast(cordicrotate(neg_theta, s), 'like', T.x);

                    s_rot_d  = complex(double(real(s_rot_fi)), double(imag(s_rot_fi)));
                    s_dec    = qam_slicer(s_rot_d, M);

                    err_re   = real(s_rot_d) - real(s_dec);
                    err_im   = imag(s_rot_d) - imag(s_dec);
                    m_buf(b) = m_buf(b) + err_re * err_re + err_im * err_im;
                end
            end

            [~, best_b]    = min(m_buf);
            Thetas(i, pol) = ThetaTest_fi(best_b);
        end
    end

    %% ================================================================
    %  Step 2: Phase unwrapping + pilot cycle-slip correction
    %           + optional block-constant phase
    %
    %  Unwrap formula:
    %    n        = floor(0.5 - (Thetas(i) - ThetaPrev) / (pi/2))
    %    theta_uw = Thetas(i) + n * (pi/2)
    %% ================================================================
    ThetaPU   = zeros(Nsym, NPol, 'like', T.theta);
    ThetaPrev = zeros(1, NPol,    'like', T.theta);

    for i = 1:Nsym
        for pol = 1:NPol

            diff_val = Thetas(i, pol) - ThetaPrev(pol);
            n_val    = floor(0.5 - double(diff_val) / double(PI_OVER2));
            n_fi     = cast(n_val, 'like', T.theta);
            theta_uw = Thetas(i, pol) + n_fi * PI_OVER2;

            % theta_uw stays within the fi range (±8 rad for fixed16 FL=12)
            % because Thetas(i) ∈ (-pi/4, pi/4] and n_val is 0 or ±1,
            % so theta_uw never exceeds ~±3pi/4 < ±8.  No wrap needed.

            if UsePilots
                BlockIdx = ceil(i / BlockLen);

                % theta_uw is unwrapped; PhiRef is wrapped to (-pi, pi]
                % by cordicangle.  Wrap the difference to (-pi, pi] to
                % obtain the shortest-path phase error, then threshold at
                % ±pi/2 to detect a cycle slip of one quadrant.
                PhaseDiff_d = mod(double(theta_uw - PhiRef(BlockIdx, pol)) ...
                                  + pi, 2*pi) - pi;

                n_slip = ZERO_TH;
                if PhaseDiff_d > pi/2
                    n_slip =  ONE_TH;
                elseif PhaseDiff_d < -pi/2
                    n_slip = -ONE_TH;
                end
                theta_uw = theta_uw - n_slip * PI_OVER2;
            end

            if BlockBased && mod(i, BlockLen) ~= 1 && i > 1
                theta_uw = ThetaPU(i - 1, pol);
            end

            ThetaPU(i, pol) = theta_uw;
            ThetaPrev(pol)  = theta_uw;
        end
    end

    %% ================================================================
    %  Step 3: Phase correction
    %    v(i, pol) = z(i, pol) * exp(-j * ThetaPU(i, pol))
    %
    %  cordicrotate output is cast explicitly to T.x to enforce
    %  SpecifyPrecision fimath rather than relying on the implicit
    %  recast performed by assignment to v.
    %% ================================================================
    v = complex(zeros(Nsym, NPol, 'like', T.x));
    for i = 1:Nsym
        for pol = 1:NPol
            neg_theta = cast(-ThetaPU(i, pol), 'like', T.theta);
            v(i, pol) = cast(cordicrotate(neg_theta, z_fi(i, pol)), 'like', T.x);
        end
    end
end