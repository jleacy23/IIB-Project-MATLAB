function y = combined_cd_td_gardner_adaptive(In, SpS, NTapCD, ...
        D, L, CLambda, Rs, Rolloff, Span, ki, kp, NSymb, NLanes, AdaptOpts)
%COMBINED_CD_TD_GARDNER_ADAPTIVE  Time-domain CD + RRC matched filter as a
%   single per-polarisation FIR, then Gardner DPLL timing recovery, then
%   butterfly CMA equaliser.
%
%   y = combined_cd_td_gardner_adaptive(In, SpS, NTapCD, D, L, CLambda, ...
%           Rs, Rolloff, Span, ki, kp, NSymb, NLanes, AdaptOpts)
%
%   Inputs
%     In         - input signal [samples x 2]
%     SpS        - samples per symbol (must be 2 for Gardner)
%     NTapCD     - chromatic-dispersion FIR length
%     D, L, CLambda, Rs - dispersion / signal parameters
%     Rolloff, Span     - RRC matched filter parameters
%     ki, kp     - Gardner DPLL loop-filter gains
%     NSymb      - number of transmitted symbols
%     NLanes     - clk_recovery parallel lanes per block
%     AdaptOpts  - adaptive equaliser settings struct (NTaps, Mu, ...) -
%                  see eq_clk.apply_adaptive_eq.
%
%   The CD chirp taps and the RRC matched filter taps are convolved into a
%   single FIR per polarisation, applied via 'same' linear convolution.
%   The combined FIR is trivially parallelisable (no inter-tap dependency).

    if SpS ~= 2
        error('combined_cd_td_gardner_adaptive:SpS', ...
              'Gardner DPLL requires SpS = 2.');
    end

    NPol = size(In, 2);

    %% Build the combined CD + matched-filter FIR --------------------
    gCD = eq_clk.cd_fir_taps(D, L, CLambda, Rs, SpS, NTapCD);
    hMF = rcosdesign(Rolloff, Span, SpS, 'sqrt').';
    hMF = hMF / sum(abs(hMF).^2);

    gComb = conv(gCD, hMF);   % length(gCD) + length(hMF) - 1

    %% Apply per polarisation ---------------------------------------
    z = zeros(size(In));
    for p = 1:NPol
        z(:, p) = conv(In(:, p), gComb, 'same');
    end

    %% Gardner timing recovery (per polarisation) -------------------
    zClk = run_gardner_per_pol(z, NSymb, ki, kp, NLanes);

    %% Adaptive butterfly CMA --------------------------------------
    if NPol == 1
        zClk = [zClk, zClk];
        singlePol = true;
    else
        singlePol = false;
    end
    y = eq_clk.apply_adaptive_eq(zClk, SpS, AdaptOpts);
    if singlePol
        y = y(:, 1);
    end
end


% =====================================================================
function z = run_gardner_per_pol(x, NSymb, ki, kp, NLanes)
    NPol = size(x, 2);
    outCols = cell(1, NPol);
    minLen = inf;
    for p = 1:NPol
        out = clk_recovery.recovery(x(:, p), 'NRZ', NSymb, ki, kp, NLanes);
        outCols{p} = out;
        minLen = min(minLen, length(out));
    end
    z = zeros(minLen, NPol);
    for p = 1:NPol
        z(:, p) = outCols{p}(1:minLen);
    end
end
