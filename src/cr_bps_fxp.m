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
%     M          - QAM order, e.g. 4, 16, 64 (double, compile-time constant)
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
%       tables and no dependency on Communications Toolbox.  M may be a
%       runtime variable; no coder.const is required.

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
    PI_VAL   = cast(pi,   'like', T.theta);
    TWO_PI   = cast(2*pi, 'like', T.theta);
    ZERO_TH  = cast(0,    'like', T.theta);
    ONE_TH   = cast(1,    'like', T.theta);
    ZERO_ACC = cast(0,    'like', T.acc);

    %% ----------------------------------------------------------------
    %  Dimensions
    %% ----------------------------------------------------------------
    Nsym    = size(z, 1);
    L       = 2 * N + 1;            % BPS averaging window length
    halfL   = N;                     % = floor(L/2)
    P       = length(Pilots);
    NBlocks = ceil(Nsym / BlockLen);

    %% ----------------------------------------------------------------
    %  Cast inputs to fixed-point
    %% ----------------------------------------------------------------
    z_fi      = cast(z,       'like', T.x);
    Pilots_fi = cast(Pilots,  'like', T.x);

    %% ----------------------------------------------------------------
    %  Precompute test phase angles as fi
    %
    %  Original:  b_vec = -B/2 : B/2-1
    %             ThetaTest = (pi/2) * b_vec / B
    %
    %  Loop index b = 1..B  maps to  b_val = b - 1 - B/2  = -B/2..B/2-1
    %
    %  Rotation is applied per sample via cordicrotate(-ThetaTest_fi(b), s)
    %  rather than a complex multiply, matching hardware CORDIC datapaths.
    %% ----------------------------------------------------------------
    ThetaTest_fi = zeros(1, B, 'like', T.theta);

    for b = 1:B
        b_val   = (b - 1) - B/2;            % -B/2 .. B/2-1
        theta_d = (pi/2) * b_val / B;
        ThetaTest_fi(b) = cast(theta_d, 'like', T.theta);
    end

    %% ----------------------------------------------------------------
    %  Pilot correlation  -->  PhiRef [NBlocks x NPol]
    %
    %  Replicates original behaviour: correlation is accumulated across
    %  ALL polarisations before taking angle(), yielding one shared
    %  reference phase per block (coherent combining across pols).
    %% ----------------------------------------------------------------
    PhiRef = zeros(NBlocks, NPol, 'like', T.theta);

    if UsePilots
        for blk = 1:NBlocks
            blockStart = (blk - 1) * BlockLen + 1;

            % Accumulate conj(Pilots) .* z  across pols in T.acc precision
            corr_re = ZERO_ACC;
            corr_im = ZERO_ACC;

            for pol = 1:NPol
                for p = 1:P
                    idx = blockStart + p - 1;
                    if idx >= 1 && idx <= Nsym
                        rx = z_fi(idx, pol);

                        % conj(pilot) * rx in acc precision
                        pilot_re =  cast(real(Pilots_fi(p)), 'like', T.acc);
                        pilot_im = -cast(imag(Pilots_fi(p)), 'like', T.acc);  % conjugate
                        rx_re    =  cast(real(rx),           'like', T.acc);
                        rx_im    =  cast(imag(rx),           'like', T.acc);

                        corr_re = corr_re + pilot_re * rx_re - pilot_im * rx_im;
                        corr_im = corr_im + pilot_re * rx_im + pilot_im * rx_re;
                    end
                end
            end

            % Extract angle via CORDIC (replaces floating-point angle())
            corr_fi = complex(cast(corr_re, 'like', T.theta), ...
                              cast(corr_im, 'like', T.theta));
            phi     = cordicangle(corr_fi);

            % Broadcast same reference to all pols (matches original)
            for pol = 1:NPol
                PhiRef(blk, pol) = phi;
            end
        end
    end

    %% ================================================================
    %  Step 1: BPS phase estimation
    %
    %  For each output symbol i and polarisation:
    %    (a) Form an L-sample window centred on i (zero-padded at edges).
    %    (b) For each of B test phases, rotate the L samples (fi multiply),
    %        cast to double, make QAM decisions, accumulate squared error.
    %    (c) Choose the test phase with minimum metric -> Thetas(i, pol).
    %
    %  The metric accumulator m_buf is a plain double [1 x B] vector.
    %  Keeping this in double:
    %    - prevents fi overflow from summing L squared-magnitude differences
    %    - is consistent with the decisions (which are inherently double)
    %    - avoids the need for a very wide T.acc just for this sum
    %% ================================================================
    Thetas  = zeros(Nsym, NPol, 'like', T.theta);
    m_buf   = zeros(1, B);          % double metric buffer, re-used each symbol

    for i = 1:Nsym
        for pol = 1:NPol

            % Reset metric buffer for this symbol/pol
            for b = 1:B
                m_buf(b) = 0.0;
            end

            % Accumulate squared-error metric over the L-sample window
            for k = 1:L
                idx = i - halfL - 1 + k;   % tap index into z (1-based)
                if idx >= 1 && idx <= Nsym
                    s = z_fi(idx, pol);
                else
                    s = complex(cast(0, 'like', T.x), cast(0, 'like', T.x));
                end

                for b = 1:B
                    % Rotate sample by test phase via CORDIC.
                    % cordicrotate computes exp(-j*theta)*s using the
                    % CORDIC algorithm, matching hardware rotation datapaths.
                    % The result type follows the input fi type T.x.
                    s_rot_fi = cordicrotate(-ThetaTest_fi(b), s);

                    % Cast rotated fi sample to double for the decision step.
                    % qam_slicer implements nearest-neighbour QAM decisions
                    % using pure arithmetic — no lookup tables, no coder.const
                    % requirement on M.
                    s_rot_d = complex(double(real(s_rot_fi)), double(imag(s_rot_fi)));
                    s_dec   = qam_slicer(s_rot_d, M);

                    % Squared Euclidean distance, accumulated in double
                    err_re     = real(s_rot_d) - real(s_dec);
                    err_im     = imag(s_rot_d) - imag(s_dec);
                    m_buf(b)   = m_buf(b) + err_re * err_re + err_im * err_im;
                end
            end

            % Best test phase index (argmin over B)
            [~, best_b]    = min(m_buf);
            Thetas(i, pol) = ThetaTest_fi(best_b);
        end
    end

    %% ================================================================
    %  Step 2: Phase unwrapping + pilot cycle-slip correction
    %           + optional block-constant phase
    %
    %  Unwrap formula (matches original exactly):
    %    n        = floor(0.5 - (Thetas(i) - ThetaPrev) / (pi/2))
    %    theta_uw = Thetas(i) + n * (pi/2)
    %% ================================================================
    ThetaPU   = zeros(Nsym, NPol, 'like', T.theta);
    ThetaPrev = zeros(1, NPol,    'like', T.theta);   % running phase per pol

    for i = 1:Nsym
        for pol = 1:NPol

            % Standard pi/2-ambiguity unwrap
            diff_val = Thetas(i, pol) - ThetaPrev(pol);
            n_val    = floor(0.5 - double(diff_val) / double(PI_OVER2));
            n_fi     = cast(n_val, 'like', T.theta);
            theta_uw = Thetas(i, pol) + n_fi * PI_OVER2;

            % Wrap result to [-pi, pi] to prevent fi overflow
            wrap_n   = cast(floor(double(theta_uw + PI_VAL) / double(TWO_PI)), ...
                            'like', T.theta);
            theta_uw = theta_uw - TWO_PI * wrap_n;

            % Pilot-aided cycle-slip correction (per polarisation)
            if UsePilots
                BlockIdx  = ceil(i / BlockLen);
                PhaseDiff = theta_uw - PhiRef(BlockIdx, pol);

                n_slip = ZERO_TH;
                if PhaseDiff > PI_OVER2
                    n_slip =  ONE_TH;
                elseif PhaseDiff < -PI_OVER2
                    n_slip = -ONE_TH;
                end
                theta_uw = theta_uw - n_slip * PI_OVER2;
            end

            % Block-constant phase: repeat first-sample estimate across block
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
    %  via CORDIC rotation (replaces floating-point exp(-j*.) multiply).
    %% ================================================================
    v = complex(zeros(Nsym, NPol, 'like', T.x));
    for i = 1:Nsym
        for pol = 1:NPol
            v(i, pol) = cordicrotate(-ThetaPU(i, pol), z_fi(i, pol));
        end
    end
end