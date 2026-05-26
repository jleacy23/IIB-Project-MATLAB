function y = combined_cd_fd_gardner_adaptive(In, SpS, NFFT, NOverlap, ...
        D, L, CLambda, Rs, Rolloff, ki, kp, NSymb, NLanes, AdaptOpts)
%COMBINED_CD_FD_GARDNER_ADAPTIVE  Frequency-domain CD + RRC matched filter
%   (overlap-save), then time-domain Gardner DPLL timing recovery, then
%   butterfly CMA equaliser.
%
%   y = combined_cd_fd_gardner_adaptive(In, SpS, NFFT, NOverlap, D, L, ...
%           CLambda, Rs, Rolloff, ki, kp, NSymb, NLanes, AdaptOpts)
%
%   Inputs
%     In          - input signal [samples x 2]
%     SpS         - samples per symbol (must be 2 for Gardner)
%     NFFT        - overlap-save FFT block size
%     NOverlap    - overlap-save overlap length (even, < NFFT)
%     D, L, CLambda, Rs - dispersion / signal parameters
%     Rolloff     - RRC matched-filter roll-off
%     ki, kp      - Gardner DPLL loop-filter gains
%     NSymb       - number of transmitted symbols
%     NLanes      - clk_recovery parallel lanes per block
%     AdaptOpts   - adaptive equaliser settings struct (NTaps, Mu, ...) -
%                   see eq_clk.apply_adaptive_eq.

    if SpS ~= 2
        error('combined_cd_fd_gardner_adaptive:SpS', ...
              'Gardner DPLL requires SpS = 2.');
    end

    NPol = size(In, 2);

    %% Combined frequency-domain CD + matched-filter mask -----------
    HCD = eq_clk.cd_fd_response(D, L, CLambda, Rs, SpS, NFFT);
    HMF = eq_clk.rrc_fd_response(Rolloff, NFFT, SpS);
    H   = HCD .* HMF;

    z = eq_clk.overlap_save_apply(In, H, NFFT, NOverlap);

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
