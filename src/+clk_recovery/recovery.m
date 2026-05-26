function Out = recovery(In, PSType, NSymb, ki, kp, NLanes)
%recovery  Parallelised clock recovery using a DPLL (per-lane interpolator
%   + block-combined TED + loop filter + NCO).
%
%   Out = recovery(In, PSType, NSymb, ki, kp, NLanes)
%
%   NLanes samples are produced per block. Each lane runs its own cubic
%   Farrow interpolator with a per-lane fractional interval mun_l, but
%   all lanes within a block share the same loop-filter output Wk. The
%   per-lane (mn_l, mun_l, Etamn_l) are obtained by propagating the NCO
%   recurrence with the *current* Wk across the lanes of the block.
%
%   The timing error driving the loop filter is the sum of the TED
%   contributions of every lane in the block that lands on an odd output
%   index. The first lane's TED uses Out(n-1) and Out(n-2) from the
%   previous block (held in the output buffer).
%
%   Inputs
%     In     - input signal at 2 Sa/symbol (column vector, one polarisation)
%     PSType - pulse-shaping type: 'NRZ' (Gardner TED) or 'Nyquist' (MGTED)
%     NSymb  - number of transmitted symbols (output limited to NSymb*2)
%     ki     - integral constant of the loop filter
%     kp     - proportional constant of the loop filter
%     NLanes - number of parallel lanes per block (default 2)
%
%   Output
%     Out    - clock-recovered signal (column vector, 2 Sa/symbol)

    if nargin < 6 || isempty(NLanes)
        NLanes = 2;
    end

    Etamn = 0.5;
    Wk = 1; LF_I = Wk;
    mun = 0; n = 3;
    mn = n;
    Out = zeros(1, length(In));
    Out(1:3) = In(1:3);
    LIn = length(In);

    while mn <= LIn
        % Per-lane NCO state propagated with the current (block-constant) Wk
        EtamnL = Etamn;
        munL   = mun;
        mnL    = mn;

        % --- Per-lane interpolation -----------------------------------
        for l = 0:NLanes-1
            if mnL >= LIn
                In_p1 = 0;
            else
                In_p1 = In(mnL+1);
            end
            if mnL-2 < 1 || mnL > LIn
                % Out of input range: hold previous output sample
                Out(n+l) = Out(max(1, n+l-1));
            else
                Out(n+l) = In(mnL-2)*(-1/6*munL^3 + 1/6*munL) + ...
                           In(mnL-1)*( 1/2*munL^3 + 1/2*munL^2 - munL) + ...
                           In(mnL)  *(-1/2*munL^3 - munL^2 + 1/2*munL + 1) + ...
                           In_p1    *( 1/6*munL^3 + 1/2*munL^2 + 1/3*munL);
            end

            % Advance NCO recurrence to the next lane using current Wk
            if -1 < (EtamnL - Wk) && (EtamnL - Wk) < 0
                mnL = mnL + 1;
            elseif (EtamnL - Wk) >= 0
                mnL = mnL + 2;
            end
            EtamnL = mod(EtamnL - Wk, 1);
            munL   = EtamnL / Wk;
        end

        % --- Combined block error (TED fires on odd output indices) ---
        ek_sum = 0;
        for l = 0:NLanes-1
            nL = n + l;
            if mod(nL, 2) == 1 && nL >= 3 && nL <= length(Out)
                switch PSType
                    case 'Nyquist'
                        ek_sum = ek_sum + abs(Out(nL-1))^2 * ...
                                          (abs(Out(nL-2))^2 - abs(Out(nL))^2);
                    case 'NRZ'
                        ek_sum = ek_sum + real(conj(Out(nL-1)) * ...
                                               (Out(nL) - Out(nL-2)));
                end
            end
        end

        % --- Loop filter update (once per block) ----------------------
        LF_I = ki*ek_sum + LF_I;
        LF_P = kp*ek_sum;
        Wk   = LF_P + LF_I;

        % --- Commit block-final NCO state -----------------------------
        mn    = mnL;
        mun   = munL;
        Etamn = EtamnL;
        n     = n + NLanes;
    end

    % Limit output length to NSymb*2
    if NSymb*2 < length(Out)
        Out = Out(1:NSymb*2).';
    else
        Out = Out.';
    end
end
