function [v, ThetaPU] = cr_viterbiViterbi_fxp(x, NPol, NTaps, VVFilter, ...
                                                Pilots, L, UsePilots, BlockBased, T) %#codegen
%CR_VITERBIVITERBI_FXP  Fixed-point Viterbi-Viterbi carrier phase recovery
%                        with optional pilot-aided cycle-slip correction and
%                        block-constant phase output.
%
%   [v, ThetaPU] = cr_viterbiViterbi_fxp(x, NPol, NTaps, VVFilter, ...
%                      Pilots, L, UsePilots, BlockBased, T)
%
%   Inputs
%     x          - input signal [N x NPol] (fi or double, complex)
%     NPol       - number of polarisations (double scalar)
%     NTaps      - one-sided filter half-length (double scalar)
%                  Filter window = 2*NTaps+1 samples
%     VVFilter   - VV filter coefficients [(2*NTaps+1) x 1] (fi or double, real)
%     Pilots     - pilot symbols used at block start [P x 1] (complex)
%                  Ignored when UsePilots = false (pass [] or zeros)
%     L          - block length in symbols (double scalar)
%     UsePilots  - logical: enable pilot-aided cycle-slip correction
%     BlockBased - logical: hold phase constant over each L-symbol block
%     T          - (optional) fixed-point types table from
%                  cr_viterbiViterbi_fxp_types.  Defaults to 'fixed16'.
%
%   Outputs
%     v          - phase-corrected signal [N x NPol], same type as T.x
%     ThetaPU    - unwrapped (and optionally pilot-corrected) phase
%                  estimate [N x NPol], same type as T.theta
%
%   Fixed-point types table T must supply:
%     T.x      - input / output signal type
%     T.w      - filter coefficient type
%     T.theta  - phase / angle type  (must accommodate ±pi)
%     T.acc    - accumulator type    (used for pilot correlation sums)
%
%   Codegen notes
%     - No convmtx: phase estimation uses explicit tap-delay indexing.
%     - cordicangle replaces angle() for phase extraction.
%     - cordicrotate replaces exp(-j*theta)*x for phase correction.
%     - All if-branches on UsePilots / BlockBased are runtime-legal so
%       codegen compiles both paths without coder.const wrapping.
%     - PhiRef and ThetaPU are pre-allocated before all loops.

    %% ----------------------------------------------------------------
    %  Default types table
    %% ----------------------------------------------------------------
    if nargin < 9 || isempty(T)
        T = cr_viterbiViterbi_fxp_types('fixed16');
    end

    %% ----------------------------------------------------------------
    %  Fixed-point constants
    %% ----------------------------------------------------------------
    PI_VAL   = cast(pi,    'like', T.theta);
    PI_OVER2 = cast(pi/2,  'like', T.theta);
    PI_OVER4 = cast(pi/4,  'like', T.theta);
    TWO_PI   = cast(2*pi,  'like', T.theta);
    QUARTER  = cast(0.25,  'like', T.theta);
    ZERO_TH  = cast(0,     'like', T.theta);
    ONE_TH   = cast(1,     'like', T.theta);
    ZERO_ACC = cast(0,     'like', T.acc);

    %% ----------------------------------------------------------------
    %  Dimensions
    %% ----------------------------------------------------------------
    N      = size(x, 1);
    L_filt = 2 * NTaps + 1;
    halfL  = floor(L_filt / 2);
    P      = length(Pilots);        % number of pilot symbols per block
    NBlocks = ceil(N / L);

    %% ----------------------------------------------------------------
    %  Cast inputs to fixed-point
    %% ----------------------------------------------------------------
    x_fi = cast(x,        'like', T.x);
    w_fi = cast(VVFilter, 'like', T.w);

    %% ----------------------------------------------------------------
    %  Pre-allocate outputs and working arrays
    %% ----------------------------------------------------------------
    ThetaML = zeros(N, NPol, 'like', T.theta);  % raw ML phase estimate
    ThetaPU = zeros(N, NPol, 'like', T.theta);  % unwrapped / corrected phase
    v       = complex(zeros(N, NPol, 'like', T.x));

    % Pilot-derived block phase references (one scalar per block per pol)
    PhiRef  = zeros(NBlocks, NPol, 'like', T.theta);

    %% ================================================================
    %  Pilot correlation  (runs only when UsePilots == true)
    %  PhiRef(b,pol) = angle( sum_p  conj(Pilots(p)) * x(blockStart+p-1, pol) )
    %% ================================================================
    if UsePilots
        Pilots_fi = cast(Pilots, 'like', T.x);

        for pol = 1:NPol
            for b = 1:NBlocks
                blockStart = (b - 1) * L + 1;

                % Accumulate real and imaginary parts separately so we
                % stay within T.acc precision throughout the sum.
                corr_re = ZERO_ACC;
                corr_im = ZERO_ACC;

                for p = 1:P
                    idx = blockStart + p - 1;
                    if idx >= 1 && idx <= N
                        rx = x_fi(idx, pol);

                        % conj(pilot) * rx  in acc precision
                        pilot_re =  real(Pilots_fi(p));
                        pilot_im = -imag(Pilots_fi(p));   % conjugate

                        rx_re = cast(real(rx), 'like', T.acc);
                        rx_im = cast(imag(rx), 'like', T.acc);
                        pr    = cast(pilot_re, 'like', T.acc);
                        pi_c  = cast(pilot_im, 'like', T.acc);

                        % (pilot_re + j*pilot_im_conj)(rx_re + j*rx_im)
                        corr_re = corr_re + pr * rx_re - pi_c * rx_im;
                        corr_im = corr_im + pr * rx_im + pi_c * rx_re;
                    end
                end

                % Convert accumulator to theta type for cordicangle
                corr_fi = complex( cast(corr_re, 'like', T.theta), ...
                                   cast(corr_im, 'like', T.theta) );
                PhiRef(b, pol) = cordicangle(corr_fi);
            end
        end
    end

    %% ================================================================
    %  Per-polarisation VV phase estimation  +  unwrap / correction
    %% ================================================================
    for pol = 1:NPol

        %% ------------------------------------------------------------
        %  Step 1 – 4th-power FIR phase estimation
        %  ThetaML(i) = angle( sum_k  w(k) * x(i+k-halfL)^4 ) / 4 - pi/4
        %% ------------------------------------------------------------
        for i = 1:N
            sum4_re = ZERO_TH;
            sum4_im = ZERO_TH;

            for k = 1:L_filt
                idx = i - halfL - 1 + k;
                if idx >= 1 && idx <= N
                    s = x_fi(idx, pol);
                else
                    s = complex(cast(0, 'like', T.x), cast(0, 'like', T.x));
                end

                % s^4 via two complex squarings (avoids ^4 on fi)
                s_re = cast(real(s), 'like', T.theta);
                s_im = cast(imag(s), 'like', T.theta);

                % s^2
                s2_re = s_re * s_re - s_im * s_im;
                s2_im = cast(2, 'like', T.theta) * s_re * s_im;

                % s^4 = (s^2)^2
                s4_re = s2_re * s2_re - s2_im * s2_im;
                s4_im = cast(2, 'like', T.theta) * s2_re * s2_im;

                % Weighted accumulation  (weight in T.theta precision)
                w_k = cast(w_fi(k), 'like', T.theta);
                sum4_re = sum4_re + w_k * s4_re;
                sum4_im = sum4_im + w_k * s4_im;
            end

            sum4_fi   = complex(sum4_re, sum4_im);
            theta4    = cordicangle(sum4_fi);                 % in (-pi, pi]
            ThetaML(i, pol) = cast(theta4, 'like', T.theta) * QUARTER - PI_OVER4;
        end

        %% ------------------------------------------------------------
        %  Step 2 – Phase unwrapping + pilot cycle-slip correction
        %           + optional block-constant phase
        %% ------------------------------------------------------------
        ThetaPrev = ZERO_TH;

        for i = 1:N

            %-- Standard quadrant-ambiguity unwrap --------------------
            diff_val = ThetaPrev - ThetaML(i, pol);
            n_val    = floor(double(diff_val) / double(PI_OVER2) + 0.5);
            n_fi     = cast(n_val, 'like', T.theta);
            theta_uw = ThetaML(i, pol) + n_fi * PI_OVER2;

            % Wrap to [-pi, pi] to prevent fixed-point overflow
            wrap_n   = cast(floor(double(theta_uw + PI_VAL) / double(TWO_PI)), ...
                            'like', T.theta);
            theta_uw = theta_uw - TWO_PI * wrap_n;

            %-- Pilot-aided cycle-slip correction ---------------------
            if UsePilots
                BlockIdx  = ceil(i / L);          % 1-based block index
                PhaseDiff = theta_uw - PhiRef(BlockIdx, pol);

                % Determine integer number of pi/2 slips
                n_slip = ZERO_TH;
                if PhaseDiff > PI_OVER2
                    n_slip = ONE_TH;
                elseif PhaseDiff < -PI_OVER2
                    n_slip = -ONE_TH;
                end

                theta_uw = theta_uw - n_slip * PI_OVER2;
            end

            %-- Block-constant phase (hold first sample of block) -----
            if BlockBased && mod(i, L) ~= 1 && i > 1
                % Overwrite with the phase already stored for this block.
                % ThetaPU(i-1) was either the block's first-sample estimate
                % (if that was also block-held) or the previous symbol's
                % corrected phase — either way we replicate it.
                theta_uw = ThetaPU(i - 1, pol);
            end

            % Store and update running reference
            ThetaPU(i, pol) = theta_uw;
            ThetaPrev       = theta_uw;
        end

        %% ------------------------------------------------------------
        %  Step 3 – Phase correction via CORDIC rotation
        %           v(i) = x(i) * exp(-j * ThetaPU(i))
        %% ------------------------------------------------------------
        for i = 1:N
            v(i, pol) = cordicrotate(-ThetaPU(i, pol), x_fi(i, pol));
        end

    end  % for pol
end