function v = cr_viterbiViterbi_fxp(x, NPol, NTaps, VVFilter, T) %#codegen
%CR_VITERBIVITERBI_FXP  Fixed-point Viterbi-Viterbi carrier phase recovery.
%
%   v = cr_viterbiViterbi_fxp(x, NPol, NTaps, VVFilter, T)
%
%   Inputs
%     x            - input signal [samples x NPol] (fi or double, complex)
%     NPol         - number of polarisations (double scalar)
%     NTaps        - number of past/future symbols (double scalar)
%     VVFilter     - VV filter coefficients [L_filt x 1] (fi or double, real)
%     T            - (optional) fixed-point types table from
%                    cr_viterbiViterbi_fxp_types.  If omitted, uses 'fixed16'.
%
%   The types table T must supply prototype fi objects for:
%     T.x      - input / output signal type
%     T.w      - filter coefficient type
%     T.theta  - phase angle type (must fit \pm pi)
%     T.acc    - accumulator type
%
%   Implementation notes for codegen:
%     - convmtx is replaced with explicit tap-delay indexing.
%     - cordicangle extracts the phase of the complex 4th-power sum.
%     - cordicrotate applies exp(-j*theta) for phase correction.
%     - Phase unwrapping keeps a running phase, bounded to [-pi, pi].
%     - All arithmetic uses SpecifyPrecision fimath so every product
%       and sum is truncated to the same WL/FL — no bit growth.

    %% Default types table
    if nargin < 5 || isempty(T)
        T = cr_viterbiViterbi_fxp_types('fixed16');
    end

    %% Constants (cast to fi)
    PI_VAL    = cast(pi,     'like', T.theta);
    PI_OVER2  = cast(pi/2,  'like', T.theta);
    PI_OVER4  = cast(pi/4,  'like', T.theta);
    TWO_PI    = cast(2*pi,  'like', T.theta);
    QUARTER   = cast(0.25,  'like', T.theta);
    ZERO_TH   = cast(0,     'like', T.theta);

    %% Dimensions
    N       = size(x, 1);
    L_filt  = 2 * NTaps + 1;
    halfL   = floor(L_filt / 2);

    %% Cast inputs to fixed-point
    x_fi      = cast(x,        'like', T.x);
    w_fi      = cast(VVFilter, 'like', T.w);

    %% Allocate outputs and working arrays
    ThetaML = zeros(N, NPol, 'like', T.theta);
    v       = complex(zeros(N, NPol, 'like', T.x));

    %% ====================================================================
    %  Per-polarisation processing
    %  ====================================================================
    for pol = 1:NPol

        % --- Phase estimation via 4th-power + FIR filter ---
        %  For each output sample i, accumulate the weighted complex
        %  4th-power sum, then extract its angle via cordicangle.
        %    sum4 = sum_k  w(k) * xBlock(k)^4
        %    ThetaML4 = cordicangle(sum4)
        %    ThetaML  = ThetaML4 / 4 - pi/4

        for i = 1:N
            sum4 = complex(ZERO_TH, ZERO_TH);

            for k = 1:L_filt
                % Index into zero-padded signal
                idx = i - halfL - 1 + k;
                if idx >= 1 && idx <= N
                    s = x_fi(idx, pol);
                else
                    s = complex(cast(0, 'like', T.x));
                end

                % s^4 = (s^2)^2  (two complex multiplications)
                s2 = s * s;
                s4 = s2 * s2;

                % Weighted complex accumulation
                sum4 = sum4 + cast(w_fi(k), 'like', T.theta) * cast(s4, 'like', T.theta);
            end

            % cordicangle extracts angle from complex fi in (-pi, pi]
            theta4 = cordicangle(sum4);

            % ThetaML = theta4 / 4 - pi/4
            ThetaML(i, pol) = cast(theta4, 'like', T.theta) * QUARTER - PI_OVER4;
        end

        % --- Phase unwrapping ---
        ThetaPrev = ZERO_TH;
        for i = 1:N
            diff_val  = ThetaPrev - ThetaML(i, pol);
            n_val     = floor(double(diff_val) / double(PI_OVER2) + 0.5);
            n_fi      = cast(n_val, 'like', T.theta);
            theta_uw  = ThetaML(i, pol) + n_fi * PI_OVER2;

            % Wrap to [-pi, pi] to keep bounded in fixed-point
            wrap_n   = cast(floor(double(theta_uw + PI_VAL) / double(TWO_PI)), 'like', T.theta);
            theta_uw = theta_uw - TWO_PI * wrap_n;

            ThetaPrev = theta_uw;

            % Phase correction: v = x * exp(-j * theta) via CORDIC rotate
            v(i, pol) = cordicrotate(-theta_uw, x_fi(i, pol));
        end
    end
end
