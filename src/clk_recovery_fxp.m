function Out = clk_recovery_fxp(In, PSType, NSymb, ki, kp, T) %#codegen
%CLK_RECOVERY_FXP  Fixed-point clock recovery (DPLL: interpolator + TED + loop filter + NCO).
%
%   Out = clk_recovery_fxp(In, PSType, NSymb, ki, kp, T)
%
%   Inputs
%     In     - input signal at 2 Sa/symbol (column vector, one pol)
%     PSType - pulse-shaping type: 'NRZ' (Gardner TED) or 'Nyquist' (MGTED)
%     NSymb  - number of transmitted symbols (output limited to NSymb*2)
%     ki     - integral constant of the loop filter
%     kp     - proportional constant of the loop filter
%     T      - fixed-point types table from clk_recovery_fxp_types
%
%   Output
%     Out    - clock-recovered signal (column vector, 2 Sa/symbol)
%
%   Codegen / fixed-point notes
%     - All arithmetic uses types from T (cast … 'like' …).
%     - while loop replaced with for + early break for codegen.
%     - convmtx-free: direct sample indexing.
%     - fimath with SpecifyPrecision prevents bit growth.

    if nargin < 6 || isempty(T)
        T = clk_recovery_fxp_types('fixed32');
    end

    LIn = int32(length(In));

    % --- Cast scalar loop-filter and interpolator constants ---
    ki_fi = cast(ki, 'like', T.mu);
    kp_fi = cast(kp, 'like', T.mu);

    % Interpolator polynomial coefficients (Farrow cubic)
    c1_6  = cast(1/6,  'like', T.coeff);
    c1_2  = cast(1/2,  'like', T.coeff);
    c1_3  = cast(1/3,  'like', T.coeff);
    cOne  = cast(1,    'like', T.coeff);

    % --- State variables ---
    Etamn = cast(0.5,  'like', T.mu);
    Wk    = cast(1,    'like', T.mu);
    LF_I  = cast(1,    'like', T.acc);
    mun   = cast(0,    'like', T.mu);

    % --- Pre-allocate output (worst case: same length as input) ---
    maxOut = int32(NSymb * 2);
    Out    = complex(zeros(maxOut, 1, 'like', T.x));

    % --- Seed first three output samples ---
    nSeed = min(int32(3), min(LIn, maxOut));
    for idx = int32(1):nSeed
        Out(idx) = cast(In(idx), 'like', T.x);
    end

    n  = int32(3);     % output index (next to write)
    mn = int32(3);     % input base-point index

    % --- Pad input by one sample for the cubic look-ahead In(mn+1) ---
    InPad = complex(zeros(LIn + 1, 1, 'like', T.x));
    for idx = int32(1):LIn
        InPad(idx) = cast(In(idx), 'like', T.x);
    end
    % InPad(LIn+1) is already zero

    % --- Main DPLL loop (bounded for codegen) ---
    for iter = int32(1):int32(LIn * 2)   %#ok<ITCM>  upper bound
        if mn > LIn
            break;
        end
        if n > maxOut
            break;
        end

        % --- Cubic interpolator (Farrow architecture) ---
        % Out(n) = In(mn-2)*(-1/6*mu^3 + 1/6*mu) ...
        %        + In(mn-1)*( 1/2*mu^3 + 1/2*mu^2 - mu) ...
        %        + In(mn)  *(-1/2*mu^3 -     mu^2 + 1/2*mu + 1) ...
        %        + In(mn+1)*( 1/6*mu^3 + 1/2*mu^2 + 1/3*mu)
        mu2 = cast(mun * mun,        'like', T.acc);
        mu3 = cast(mu2 * mun,        'like', T.acc);

        a0 = cast(-c1_6 * mu3 + c1_6 * mun,                   'like', T.acc);
        a1 = cast( c1_2 * mu3 + c1_2 * mu2 - mun,             'like', T.acc);
        a2 = cast(-c1_2 * mu3 - mu2         + c1_2 * mun + cOne, 'like', T.acc);
        a3 = cast( c1_6 * mu3 + c1_2 * mu2  + c1_3 * mun,     'like', T.acc);

        s0 = cast(InPad(mn - 2), 'like', T.x);
        s1 = cast(InPad(mn - 1), 'like', T.x);
        s2 = cast(InPad(mn),     'like', T.x);
        s3 = cast(InPad(mn + 1), 'like', T.x);

        Out(n) = cast(a0 * s0 + a1 * s1 + a2 * s2 + a3 * s3, 'like', T.x);

        % --- Timing Error Detector (odd samples only) ---
        if mod(n, int32(2)) == int32(1) && n >= int32(3)
            if strcmp(PSType, 'Nyquist')
                % Mueller–Müller-type (modified Gardner TED)
                p_prev = cast(abs(Out(n-2)), 'like', T.acc);
                p_mid  = cast(abs(Out(n-1)), 'like', T.acc);
                p_curr = cast(abs(Out(n)),   'like', T.acc);
                ek = cast(p_mid * p_mid * (p_prev * p_prev - p_curr * p_curr), ...
                     'like', T.acc);
            else
                % Gardner TED (NRZ)
                ek = cast(real(conj(Out(n-1)) * (Out(n) - Out(n-2))), ...
                     'like', T.acc);
            end

            % Loop filter (PI)
            LF_I(:) = ki_fi * ek + LF_I;
            LF_P    = cast(kp_fi * ek, 'like', T.acc);
            Wk(:)   = LF_P + LF_I;
        end

        % --- NCO: advance base-point and fractional interval ---
        diff = cast(Etamn - Wk, 'like', T.mu);
        if diff > cast(-1, 'like', T.mu) && diff < cast(0, 'like', T.mu)
            mn = mn + int32(1);
        elseif diff >= cast(0, 'like', T.mu)
            mn = mn + int32(2);
        end

        Etamn(:) = cast(Etamn - Wk, 'like', T.mu);
        % mod(Etamn, 1): keep fractional part in [0, 1)
        Etamn(:) = Etamn - cast(floor(double(Etamn)), 'like', T.mu);

        if Wk ~= cast(0, 'like', T.mu)
            mun(:) = cast(Etamn / Wk, 'like', T.mu);
        else
            mun(:) = cast(0, 'like', T.mu);
        end

        n = n + int32(1);
    end

    % --- Trim output ---
    nOut = min(n - int32(1), maxOut);
    Out  = Out(1:nOut);
end
