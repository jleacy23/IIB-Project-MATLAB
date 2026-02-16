function Out = clk_recovery(In, PSType, NSymb, ki, kp)
%CLK_RECOVERY  Clock recovery using a DPLL (interpolator + TED + loop filter + NCO).
%
%   Out = clk_recovery(In, PSType, NSymb, ki, kp)
%
%   Inputs
%     In     - input signal at 2 Sa/symbol (column vector, one polarisation)
%     PSType - pulse-shaping type: 'NRZ' (Gardner TED) or 'Nyquist' (MGTED)
%     NSymb  - number of transmitted symbols (output limited to NSymb*2)
%     ki     - integral constant of the loop filter
%     kp     - proportional constant of the loop filter
%
%   Output
%     Out    - clock-recovered signal (column vector, 2 Sa/symbol)

    Etamn = 0.5;
    Wk = 1; LF_I = Wk;
    mun = 0; n = 3;
    mn = n;
    Out = zeros(1, length(In));
    Out(1:3) = In(1:3);
    LIn = length(In);

    while mn <= LIn
        if mn == LIn
            In(mn+1) = 0;
        end

        % Cubic interpolator with Farrow architecture
        Out(n) = In(mn-2)*(-1/6*mun^3 + 1/6*mun) + ...
                 In(mn-1)*( 1/2*mun^3 + 1/2*mun^2 - mun) + ...
                 In(mn)  *(-1/2*mun^3 - mun^2 + 1/2*mun + 1) + ...
                 In(mn+1)*( 1/6*mun^3 + 1/2*mun^2 + 1/3*mun);

        % Ts-spaced timing error (odd samples only)
        if mod(n, 2) == 1
            switch PSType
                case 'Nyquist'
                    ek = abs(Out(n-1)).^2 .* ...
                         (abs(Out(n-2)).^2 - abs(Out(n)).^2);
                case 'NRZ'
                    ek = real(conj(Out(n-1)) .* (Out(n) - Out(n-2)));
            end

            % Loop filter
            LF_I = ki*ek + LF_I;
            LF_P = kp*ek;
            Wk   = LF_P + LF_I;
        end

        % NCO — base point and fractional interval
        if -1 < (Etamn - Wk) && (Etamn - Wk) < 0
            mn = mn + 1;
        elseif (Etamn - Wk) >= 0
            mn = mn + 2;
        end
        Etamn = mod(Etamn - Wk, 1);
        mun   = Etamn / Wk;

        n = n + 1;
    end

    % Limit output length to NSymb*2
    if NSymb*2 < length(Out)
        Out = Out(1:NSymb*2).';
    else
        Out = Out.';
    end
end