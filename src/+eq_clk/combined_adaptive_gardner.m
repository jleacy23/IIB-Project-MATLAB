function y = combined_adaptive_gardner(In, SpS, ki, kp, NSymb, NLanes, ...
        ClkFirst, AdaptOpts)
%COMBINED_ADAPTIVE_GARDNER  Adaptive butterfly equaliser + Gardner DPLL
%   timing recovery, with selectable ordering.
%
%   y = combined_adaptive_gardner(In, SpS, ki, kp, NSymb, NLanes, ...
%                                  ClkFirst, AdaptOpts)
%
%   Inputs
%     In        - input signal [samples x 2]
%     SpS       - samples per symbol (must be 2 for Gardner)
%     ki, kp    - Gardner DPLL loop-filter gains
%     NSymb     - number of transmitted symbols
%     NLanes    - clk_recovery parallel lanes per block
%     ClkFirst  - true (default): Gardner -> CMA (decimates to 1 Sa/sym)
%                 false:          CMA (T/2-spaced, no decimation) ->
%                                 Gardner (2 Sa/sym) -> decimate
%     AdaptOpts - struct of adaptive equaliser settings (see
%                 eq_clk.apply_adaptive_eq).  Must contain at least NTaps
%                 and Mu; all other adaptive_eq.equalize parameters are
%                 also routed through this struct.
%
%   When ClkFirst is true the stock adaptive_eq.equalize is used.  In the
%   swap case a local fractionally-spaced butterfly CMA helper
%   (fs_cma_local) preserves the 2 Sa/symbol rate; it consumes only the
%   NTaps and Mu fields of AdaptOpts.

    if nargin < 7 || isempty(ClkFirst)
        ClkFirst = true;
    end

    if ClkFirst
        % ---- Gardner first --------------------------------------------------
        z = run_gardner(In, SpS, NSymb, ki, kp, NLanes);
        y = eq_clk.apply_adaptive_eq(z, SpS, AdaptOpts);
    else
        % ---- Adaptive first (T/2-spaced, no decimation) --------------------
        z = fs_cma_local(In, AdaptOpts.NTaps, AdaptOpts.Mu);
        z = run_gardner(z, SpS, NSymb, ki, kp, NLanes);
        y = z(1:SpS:end, :);
    end
end


% =====================================================================
function z = run_gardner(x, SpS, NSymb, ki, kp, NLanes)
%RUN_GARDNER  Per-polarisation Gardner DPLL timing recovery.
    if SpS ~= 2
        error('combined_adaptive_gardner:SpS', ...
              'Gardner DPLL requires SpS = 2.');
    end

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


% =====================================================================
function y = fs_cma_local(x, NTaps, Mu)
%FS_CMA_LOCAL  T/2-spaced butterfly CMA, non-decimating.
%   Minimal fractionally spaced CMA used only when the adaptive stage
%   must precede the timing recovery.  Output rate equals input rate.

    NSamp = size(x, 1);
    halfN = floor(NTaps/2);
    xp = [x(end-halfN+1:end,:); x; x(1:NTaps-1-halfN,:)];

    w1V = zeros(NTaps, 1); w1V(halfN+1) = 1;
    w1H = zeros(NTaps, 1);
    w2V = zeros(NTaps, 1);
    w2H = zeros(NTaps, 1); w2H(halfN+1) = 1;

    R_CMA = 2;

    y = zeros(NSamp, 2);
    for n = 1:NSamp
        xV = xp(n:n+NTaps-1, 1);
        xH = xp(n:n+NTaps-1, 2);

        y1 = w1V' * xV + w1H' * xH;
        y2 = w2V' * xV + w2H' * xH;
        y(n, :) = [y1, y2];

        e1 = R_CMA - abs(y1)^2;
        e2 = R_CMA - abs(y2)^2;
        w1V = w1V + Mu * e1 * conj(y1) * xV;
        w1H = w1H + Mu * e1 * conj(y1) * xH;
        w2V = w2V + Mu * e2 * conj(y2) * xV;
        w2H = w2H + Mu * e2 * conj(y2) * xH;
    end
end
