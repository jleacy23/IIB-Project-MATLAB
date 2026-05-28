function Out = recovery_fxp(In, NSymb, ki, kp, NLanes, T) %#codegen
%RECOVERY_FXP  Fixed-point parallelised Gardner DPLL clock recovery.
%
%   Out = recovery_fxp(In, NSymb, ki, kp, NLanes, T)
%
%   Fixed-point counterpart of clk_recovery.recovery for the NRZ Gardner
%   TED only.  The 'Nyquist' (MGTED) branch is not implemented — that path
%   would need a different accumulator profile (signal magnitude is raised
%   to the fourth power) and is left for a future variant.
%
%   Inputs
%     In     - input signal at 2 Sa/symbol [Nsamp x 1] (fi or castable)
%     NSymb  - number of transmitted symbols (output limited to NSymb*2)
%     ki, kp - integral / proportional loop-filter constants (double)
%     NLanes - parallel lanes per block
%     T      - (optional) fixed-point types table from
%              clk_recovery.recovery_fxp_types.  If omitted, uses 'fixed32'.
%
%   Output
%     Out    - clock-recovered signal [Nout x 1], T.x precision
%
%   Implementation notes for codegen:
%     - The data path (signal samples and the interpolator MAC) uses fi
%       at T.x / T.acc precision.
%     - The NCO state (Etamn, mun, Wk, LF_I) is stored at fi precision
%       (T.nco / T.lf) but the in-loop recurrence steps (modulo-1, divide
%       by Wk) are computed in double — fi cannot cleanly express the
%       division by a near-unity value without large word lengths, and
%       restoring the result to T.nco bounds it back into the design range.
%     - Integer-valued bookkeeping (mn, n, l) is kept as double for
%       portable indexing semantics.

    if nargin < 5 || isempty(NLanes)
        NLanes = 2;
    end
    if nargin < 6 || isempty(T)
        T = clk_recovery.recovery_fxp_types('fixed32');
    end

    NIn = size(In, 1);

    %% Cast input
    InX = cast(In, 'like', T.x);

    %% Initial state (matches float recovery.m)
    Etamn = cast(0.5, 'like', T.nco);
    mun   = cast(0,   'like', T.nco);
    Wk    = cast(1,   'like', T.lf);
    LF_I  = cast(1,   'like', T.lf);

    ki_fi = cast(ki, 'like', T.lf);
    kp_fi = cast(kp, 'like', T.lf);

    %% Output buffer (one polarisation, column vector) — T.x precision.
    %  The float reference silently grows Out via out-of-bounds assignment
    %  when the NCO advances n past length(In) inside a block.  Codegen
    %  doesn't grow pre-allocated arrays, so allocate with a one-block
    %  margin (NLanes + small constant) and truncate to NSymb*2 at the end.
    NOutMax = NIn + NLanes + 4;
    Out = complex(zeros(NOutMax, 1, 'like', T.x));
    if NIn >= 3
        Out(1) = InX(1);
        Out(2) = InX(2);
        Out(3) = InX(3);
    end

    %% Pre-cast Farrow polynomial constants (real, in T.coef precision)
    cA = cast(-1/6, 'like', T.coef);
    cB = cast( 1/6, 'like', T.coef);
    cC = cast( 1/2, 'like', T.coef);
    cD = cast(-1/2, 'like', T.coef);
    cE = cast( 1/3, 'like', T.coef);
    cF = cast( 1,   'like', T.coef);

    %% Temporaries
    acc = complex(cast(0, 'like', T.acc));

    mn = 3;
    n  = 3;

    while mn <= NIn
        %% Per-lane state propagated with the current block-constant Wk
        EtamnL = Etamn;
        munL   = mun;
        mnL    = mn;

        %% -------- Per-lane interpolation -------------------------------
        for l = 0:NLanes-1
            if mnL >= NIn
                In_p1 = complex(cast(0, 'like', T.x));
            else
                In_p1 = InX(mnL + 1);
            end

            outIdx = n + l;
            if outIdx > NOutMax
                % NCO has overrun the allocated output buffer: skip the
                % write (the loop will exit shortly when mn > NIn).  The
                % per-lane NCO recurrence still advances below.
            elseif mnL - 2 < 1 || mnL > NIn
                % Out of input range: hold the previous output sample
                idxHold = outIdx - 1;
                if idxHold < 1
                    idxHold = 1;
                end
                Out(outIdx) = Out(idxHold);
            else
                % Cubic Farrow coefficients (real, T.coef precision)
                mu1 = munL;
                mu2 = cast(mu1 * mu1, 'like', T.coef);
                mu3 = cast(mu2 * mu1, 'like', T.coef);

                c0 = cast(cA*mu3 + cB*mu1,            'like', T.coef);
                c1 = cast(cC*mu3 + cC*mu2 - mu1,      'like', T.coef);
                c2 = cast(cD*mu3 - mu2 + cC*mu1 + cF, 'like', T.coef);
                c3 = cast(cB*mu3 + cC*mu2 + cE*mu1,   'like', T.coef);

                acc(:) = complex(cast(0, 'like', T.acc));
                acc(:) = acc + cast(c0, 'like', T.acc) * InX(mnL - 2);
                acc(:) = acc + cast(c1, 'like', T.acc) * InX(mnL - 1);
                acc(:) = acc + cast(c2, 'like', T.acc) * InX(mnL);
                acc(:) = acc + cast(c3, 'like', T.acc) * In_p1;
                Out(outIdx) = cast(acc, 'like', T.x);
            end

            %% NCO recurrence (decision + mod-1 in double, restore to fi)
            eta_d = double(EtamnL);
            wk_d  = double(Wk);
            d     = eta_d - wk_d;

            if d > -1 && d < 0
                mnL = mnL + 1;
            elseif d >= 0
                mnL = mnL + 2;
            end

            % mod(eta - Wk, 1)  in double (positive result for any d)
            eta_next_d = d - floor(d);

            % mun_next = Etamn_next / Wk
            if wk_d ~= 0
                mun_next_d = eta_next_d / wk_d;
            else
                mun_next_d = 0;
            end

            EtamnL = cast(eta_next_d, 'like', T.nco);
            munL   = cast(mun_next_d, 'like', T.nco);
        end

        %% -------- Combined block TED (NRZ Gardner) ---------------------
        %  Mirrors the float reference's `nL <= length(Out)` guard, using
        %  NOutMax as the (codegen-known) upper bound on the Out buffer.
        ek_sum = cast(0, 'like', T.ek);
        for l = 0:NLanes-1
            nL = n + l;
            if mod(nL, 2) == 1 && nL >= 3 && nL <= NOutMax
                y_now = Out(nL);
                y_pm1 = Out(nL - 1);
                y_pm2 = Out(nL - 2);
                term  = conj(y_pm1) * (y_now - y_pm2);
                ek_sum(:) = ek_sum + cast(real(term), 'like', T.ek);
            end
        end

        %% -------- Loop filter update (once per block) -----------------
        ek_lf = cast(ek_sum, 'like', T.lf);
        LF_I(:) = LF_I + ki_fi * ek_lf;
        LF_P    = kp_fi * ek_lf;
        Wk(:)   = LF_P + LF_I;

        %% Commit block-final NCO state ---------------------------------
        mn    = mnL;
        mun   = munL;
        Etamn = EtamnL;
        n     = n + NLanes;
    end

    %% Limit output length to NSymb*2 (mirrors the float reference).
    %  When NSymb*2 <= NOutMax we truncate; otherwise we keep the full
    %  allocated buffer (the float silently returns the grown length).
    if NSymb * 2 < NOutMax
        Out = Out(1 : NSymb * 2);
    end
end
