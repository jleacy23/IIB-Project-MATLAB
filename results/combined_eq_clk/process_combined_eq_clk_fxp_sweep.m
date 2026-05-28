function process_combined_eq_clk_fxp_sweep(varargin)
%PROCESS_COMBINED_EQ_CLK_FXP_SWEEP  FEC SNR + energy/bit vs fxp precision.
%
%   process_combined_eq_clk_fxp_sweep()
%   process_combined_eq_clk_fxp_sweep('MatFile', path, 'FECBER', 2e-2, ...
%                                     'Node', '45nm', 'M', 4)
%
%   Loads combined_eq_clk_fxp_sweep.mat and prints:
%     1. A long-form table — one row per (block, po2, EqFL, ClkFL) with
%        FEC SNR (mean ± std across trials) and energy per bit (per
%        stage and total).
%     2. Pivot grids per (block, po2): EqFL rows x ClkFL columns,
%        first for FEC SNR (dB), then for total energy per bit (pJ).
%
%   Energy model
%     Per-symbol real-multiplication and real-addition counts come from
%     report/full/full.tex Tables tab:cd_cost (static CD + MF
%     overlap-save), tab:aeq_cost (sign-sign butterfly CMA), and
%     tab:clk_cost (Gardner / Modified-Godard).  Each stage's energy is
%     computed via energy.receiver(NAdd, NMult, E_A, E_M, M, 1, n) with
%     n = eq_wl for the static and adaptive stages, n = clk_wl for the
%     clock-recovery stage.  Per-symbol formulas already include the
%     oversampling factor eta, so Oversampling = 1 is passed in.
%
%   Name/Value options:
%     'MatFile'   - path to the .mat (default: alongside this script)
%     'FECBER'    - FEC threshold used to score designs (default 2e-2)
%     'Node'      - '14nm' (scaled estimate, default) or '45nm' (calibrated)
%     'M'         - modulation order (default 4 for QPSK)
%     'SavePlots'   - true to write the diagonal-precision plot PNG
%                     alongside the .mat (default true)
%     'GardnerFLs'  - vector of FL values to include in the Gardner
%                     energy-breakdown table at the end (default [4, 6]).
%                     EqFL == ClkFL == FL is enforced.

    here = fileparts(mfilename('fullpath'));
    p = inputParser;
    p.addParameter('MatFile', ...
        fullfile(here, 'combined_eq_clk_fxp_sweep.mat'));
    p.addParameter('FECBER', 2e-2);
    p.addParameter('Node', '14nm');
    p.addParameter('M', 4);
    p.addParameter('SavePlots',  true);
    p.addParameter('GardnerFLs', [4, 6]);
    p.parse(varargin{:});
    matFile     = p.Results.MatFile;
    fecBer      = p.Results.FECBER;
    node        = p.Results.Node;
    M           = p.Results.M;
    savePlots   = p.Results.SavePlots;
    gardnerFLs  = p.Results.GardnerFLs;

    if ~isfile(matFile)
        error('process_combined_eq_clk_fxp_sweep:missingMat', ...
              ['Could not find %s.\n', ...
               'Run runtests(''combined_eq_clk_fxp_sweep'') first.'], ...
              matFile);
    end
    S      = load(matFile);
    tbl    = S.tbl;
    P      = S.params;
    SNR_dB = S.SNR_dB_vec(:).';

    [EAdd, EMult] = energyCoefs(node);
    fprintf('\nEnergy coefficients (%s): E_A = %.2f*n fJ, E_M = %.2f*n^2 fJ\n', ...
        node, EAdd, EMult);

    tbl = sortrows(tbl, {'block_name', 'po2', 'eq_fl', 'clk_fl'});

    %% --- Per-row FEC SNR (mean/std across trials) ----------------
    nRow = height(tbl);
    fecMean = nan(nRow, 1);
    fecStd  = nan(nRow, 1);
    for k = 1:nRow
        berMat = tbl.ber{k};                  % [NTrials x NSNR]
        nTrials = size(berMat, 1);
        perTrialSnr = nan(nTrials, 1);
        for tr = 1:nTrials
            perTrialSnr(tr) = fecSnrFromBer( ...
                SNR_dB, berMat(tr, :), fecBer);
        end
        valid = isfinite(perTrialSnr);
        if any(valid)
            fecMean(k) = mean(perTrialSnr(valid));
            fecStd(k)  = std(perTrialSnr(valid));
        end
    end

    %% --- Per-row energy per bit ----------------------------------
    %  energy.receiver returns energy per bit in the same units as
    %  (EAdd, EMult) — fJ here.
    E_static = nan(nRow, 1);
    E_adapt  = nan(nRow, 1);
    E_clk    = nan(nRow, 1);
    for k = 1:nRow
        blk    = char(tbl.block_name(k));
        po2    = tbl.po2(k);
        eqWL   = tbl.eq_wl(k);
        clkWL  = tbl.clk_wl(k);
        nCD    = tbl.n_cd(k);
        nAEQ   = tbl.n_aeq(k);
        N      = P.NFFT;
        eta    = P.SpS;
        beta   = P.Rolloff;

        [NM_s, NA_s] = staticOps(N, nCD, eta, po2);
        [NM_a, NA_a] = adaptOpsSignSign(nAEQ);
        switch blk
            case 'cd_gardner_cma'
                [NM_c, NA_c] = gardnerOps();
            case 'cd_godard_cma'
                [NM_c, NA_c] = godardOps(N, nCD, eta, beta);
            otherwise
                NM_c = NaN; NA_c = NaN;
        end

        E_static(k) = energy.receiver(NA_s, NM_s, EAdd, EMult, M, 1, eqWL);
        E_adapt(k)  = energy.receiver(NA_a, NM_a, EAdd, EMult, M, 1, eqWL);
        E_clk(k)    = energy.receiver(NA_c, NM_c, EAdd, EMult, M, 1, clkWL);
    end
    E_total = E_static + E_adapt + E_clk;

    %% --- Long-form table -----------------------------------------
    fprintf(['\n--- FEC SNR & energy/bit per configuration ', ...
             '(FEC BER = %.0e, CFO = %.2f GHz, %s, M = %d) ---\n'], ...
        fecBer, S.CFO_GHz, node, M);
    fprintf('%-18s %-4s %-6s %-6s %-7s %-7s %-8s %-8s %-8s %-8s %-8s\n', ...
        'block', 'po2', 'EqWL', 'ClkWL', 'FEC',  'std', ...
        'E_stat',  'E_aeq', 'E_clk', 'E_tot', 'units');
    for k = 1:nRow
        meanStr = '  ---';
        stdStr  = '  ---';
        if isfinite(fecMean(k))
            meanStr = sprintf('%6.2f', fecMean(k));
        end
        if isfinite(fecStd(k))
            stdStr = sprintf('%6.2f', fecStd(k));
        end
        % Energies are in fJ/bit; convert to pJ/bit for display.
        fprintf('%-18s %-4d %-6d %-6d %-7s %-7s %-8.3f %-8.3f %-8.3f %-8.3f %-8s\n', ...
            char(tbl.block_name(k)), tbl.po2(k), ...
            tbl.eq_wl(k), tbl.clk_wl(k), meanStr, stdStr, ...
            E_static(k) * 1e-3, E_adapt(k) * 1e-3, ...
            E_clk(k)    * 1e-3, E_total(k) * 1e-3, 'pJ/bit');
    end

    %% --- Pivot grids per (block, po2) -----------------------------
    eqFL_vec  = unique(tbl.eq_fl);
    clkFL_vec = unique(tbl.clk_fl);
    pairs     = unique(tbl(:, {'block_name', 'po2'}), 'rows', 'stable');

    for pi = 1:height(pairs)
        blk = pairs.block_name(pi);
        p2  = pairs.po2(pi);

        % FEC SNR grid
        fprintf('\n--- %s   po2 = %d   FEC SNR (dB) ---\n', char(blk), p2);
        printPivot(tbl, fecMean, blk, p2, eqFL_vec, clkFL_vec, '%8.2f');

        % Energy per bit grid (pJ)
        fprintf('\n--- %s   po2 = %d   Total energy (pJ/bit) ---\n', ...
            char(blk), p2);
        printPivot(tbl, E_total * 1e-3, blk, p2, eqFL_vec, clkFL_vec, '%8.3f');
    end

    %% --- Diagonal-precision plot: FEC SNR vs FL (= EqFL = ClkFL) ---
    % One line per (block, po2) implementation.  Restricts the table
    % to rows where EqFL == ClkFL so the single x-axis is unambiguous.
    diagMask = tbl.eq_fl == tbl.clk_fl;
    tbDiag   = tbl(diagMask, :);
    fecDiag  = fecMean(diagMask);
    stdDiag  = fecStd(diagMask);
    eDiag    = E_total(diagMask);

    if isempty(tbDiag)
        warning('process_combined_eq_clk_fxp_sweep:noDiag', ...
            'No rows with EqFL == ClkFL — skipping diagonal plot.');
        return;
    end

    figure('Name', 'FXP sweep: FEC SNR vs precision (EqFL = ClkFL)', ...
        'Position', [80 80 720 480]);
    ax = axes; hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

    impl   = unique(tbDiag(:, {'block_name', 'po2'}), 'rows', 'stable');
    nImpl  = height(impl);
    colors = lines(nImpl);
    for ii = 1:nImpl
        blk = impl.block_name(ii);
        p2  = impl.po2(ii);
        sel = (tbDiag.block_name == blk) & (tbDiag.po2 == p2);
        sub = tbDiag(sel, :);
        % Sort by FL so the line is monotonic in precision
        [~, ord] = sort(sub.eq_fl);
        sub      = sub(ord, :);
        fSel     = fecDiag(sel); fSel = fSel(ord);
        sSel     = stdDiag(sel); sSel = sSel(ord);

        label = sprintf('%s, po2=%d', char(blk), p2);
        errorbar(ax, sub.eq_fl, fSel, sSel, ...
            '-o', 'Color', colors(ii, :), 'LineWidth', 1.4, ...
            'MarkerSize', 5, 'CapSize', 8, ...
            'DisplayName', label);
    end
    xlabel(ax, sprintf('FL (EqFL = ClkFL),  WL = %d + FL', P.NIntBits));
    ylabel(ax, 'FEC SNR (dB)');
    title(ax, sprintf(['FEC SNR vs fxp precision (CFO = %.2f GHz, ', ...
        'BER = %.0e)'], S.CFO_GHz, fecBer));
    legend(ax, 'Location', 'best', 'Interpreter', 'none');

    if savePlots
        outFile = fullfile(here, 'combined_eq_clk_fxp_fec_snr.png');
        exportgraphics(gcf, outFile, 'Resolution', 200);
        fprintf('Saved diagonal-precision FEC SNR plot to %s\n', outFile);
    end

    %% --- Complementary plot: total energy/bit vs FL ----------------
    %  Total energy is the sum of all three stages from energy.receiver
    %  (in fJ/bit); displayed here in pJ/bit.  Word length follows the
    %  same WL = NIntBits + FL mapping as the FEC SNR plot — eq_wl is
    %  used for the static and adaptive stages and clk_wl for the clock
    %  recovery stage, but on the diagonal they coincide.
    figure('Name', 'FXP sweep: total energy/bit vs precision (EqFL = ClkFL)', ...
        'Position', [80 80 720 480]);
    axE = axes; hold(axE, 'on'); grid(axE, 'on'); box(axE, 'on');
    for ii = 1:nImpl
        blk = impl.block_name(ii);
        p2  = impl.po2(ii);
        sel = (tbDiag.block_name == blk) & (tbDiag.po2 == p2);
        sub = tbDiag(sel, :);
        [~, ord] = sort(sub.eq_fl);
        sub      = sub(ord, :);
        eSel     = eDiag(sel) * 1e-3;     % fJ -> pJ
        eSel     = eSel(ord);

        label = sprintf('%s, po2=%d', char(blk), p2);
        plot(axE, sub.eq_fl, eSel, '-o', ...
            'Color', colors(ii, :), 'LineWidth', 1.4, ...
            'MarkerSize', 5, 'DisplayName', label);
    end
    xlabel(axE, sprintf('FL (EqFL = ClkFL),  WL = %d + FL', P.NIntBits));
    ylabel(axE, 'Total energy (pJ / bit)');
    title(axE, sprintf(['Energy per bit vs fxp precision (CFO = %.2f GHz, ', ...
        '%s, M = %d)'], S.CFO_GHz, node, M));
    legend(axE, 'Location', 'best', 'Interpreter', 'none');

    if savePlots
        outFile = fullfile(here, 'combined_eq_clk_fxp_energy.png');
        exportgraphics(gcf, outFile, 'Resolution', 200);
        fprintf('Saved diagonal-precision energy plot to %s\n', outFile);
    end

    %% --- Gardner energy breakdown at requested FL values ----------
    %  Same precision everywhere (EqFL == ClkFL == FL).  One row per
    %  (po2, FL) combination, with per-stage and total energy/bit plus
    %  FEC SNR mean / std.
    fprintf(['\n--- Gardner energy breakdown (EqFL = ClkFL, %s, ', ...
             'M = %d, FEC BER = %.0e) ---\n'], node, M, fecBer);
    fprintf('%-4s %-3s %-7s %-7s %-8s %-8s %-8s %-8s %-8s\n', ...
        'po2', 'FL', 'WL', 'FEC',  'std', ...
        'E_stat', 'E_aeq', 'E_clk', 'E_tot');
    isGardner = tbl.block_name == "cd_gardner_cma";
    for fl = gardnerFLs(:).'
        for p2v = [false, true]
            mask = isGardner & ...
                   (tbl.po2    == p2v) & ...
                   (tbl.eq_fl  == fl)  & ...
                   (tbl.clk_fl == fl);
            idx = find(mask, 1);
            if isempty(idx)
                fprintf('%-4d %-3d %-7s %-7s %-8s %-8s %-8s %-8s %-8s\n', ...
                    p2v, fl, '---', '---', '---', '---', '---', '---', '---');
                continue;
            end
            fStr = '  ---'; sStr = '  ---';
            if isfinite(fecMean(idx))
                fStr = sprintf('%6.2f', fecMean(idx));
            end
            if isfinite(fecStd(idx))
                sStr = sprintf('%6.2f', fecStd(idx));
            end
            fprintf(['%-4d %-3d %-7d %-7s %-7s %-8.3f %-8.3f %-8.3f ', ...
                     '%-8.3f\n'], ...
                p2v, fl, tbl.eq_wl(idx), fStr, sStr, ...
                E_static(idx) * 1e-3, E_adapt(idx) * 1e-3, ...
                E_clk(idx)    * 1e-3, E_total(idx) * 1e-3);
        end
    end
    fprintf('Energies in pJ/bit, FEC SNR in dB.\n');
end


% =====================================================================
function printPivot(tbl, vals, blk, p2, eqFL_vec, clkFL_vec, fmt)
    fprintf('%-10s', 'EqFL\ClkFL');
    for cj = 1:numel(clkFL_vec)
        fprintf(' %8d', clkFL_vec(cj));
    end
    fprintf('\n');
    for ri = 1:numel(eqFL_vec)
        fprintf('%-10d', eqFL_vec(ri));
        for cj = 1:numel(clkFL_vec)
            mask = (tbl.block_name == blk) & ...
                   (tbl.po2        == p2)  & ...
                   (tbl.eq_fl      == eqFL_vec(ri)) & ...
                   (tbl.clk_fl     == clkFL_vec(cj));
            idx = find(mask, 1);
            if isempty(idx) || ~isfinite(vals(idx))
                fprintf(' %8s', '---');
            else
                fprintf([' ' fmt], vals(idx));
            end
        end
        fprintf('\n');
    end
end


% =====================================================================
%  Energy / operation-count helpers
% =====================================================================

function [EAdd, EMult] = energyCoefs(node)
%ENERGYCOEFS  Linear (per-bit) and quadratic (per-bit^2) energy coefs,
%   in fJ, from report/full/full.tex.  '45nm' is the calibrated fit to
%   Horowitz 2014; '14nm' applies the 75% process-scaling estimate.
    switch lower(node)
        case '45nm'
            EAdd  = 3.16;
            EMult = 3.03;
        case '14nm'
            EAdd  = 0.79;
            EMult = 0.76;
        otherwise
            error('process_combined_eq_clk_fxp_sweep:badNode', ...
                'Unknown Node ''%s'' (use ''45nm'' or ''14nm'').', node);
    end
end


function [NMult, NAdd] = staticOps(N, NCD, eta, po2)
%STATICOPS  Per-symbol RM/RA for the overlap-save CD + matched filter.
%   tab:cd_cost in report/full/full.tex.  When po2 is true, every FFT
%   twiddle multiplication is replaced by bit shifts, eliminating the
%   real multiplications inside the FFT and IFFT.
    denom = N - NCD + 1;
    if po2
        NMult = eta * 4*N / denom;
    else
        NMult = eta * (4*N*log2(N) + 4*N) / denom;
    end
    NAdd = eta * (6*N*log2(N) + 2*N) / denom;
end


function [NMult, NAdd] = adaptOpsSignSign(NTaps)
%ADAPTOPSSIGNSIGN  Per-symbol RM/RA for the parallel-lane butterfly CMA
%   with the sign-sign weight update.  tab:aeq_cost in
%   report/full/full.tex (FIR outputs + error terms + sign-sign update).
    NMult = 8*NTaps + 2;             % 8N (FIR) + 2 (err)
    NAdd  = (4*NTaps + 2) + 1 + 8*NTaps;  % FIR + err + sign-sign update
end


function [NMult, NAdd] = gardnerOps()
%GARDNEROPS  Per-symbol RM/RA for Gardner TED + cubic Farrow.
    NMult = 4 + 76;   % TED + Farrow
    NAdd  = 8 + 60;
end


function [NMult, NAdd] = godardOps(N, NCD, eta, beta)
%GODARDOPS  Per-symbol RM/RA for the Modified Godard metric, CORDIC for
%   arg(S), and the FD phase ramp.  The FFT/IFFT are reused from the
%   overlap-save CD/MF stage and so are not counted here.
    denom = N - NCD + 1;
    NMult = 4*beta*N/denom + eta/denom + 8*N*eta/denom;
    NAdd  = 4*beta*N/denom + 0          + 4*N*eta/denom;
end


% =====================================================================
function snr = fecSnrFromBer(snrVec, berVec, fecBer)
%FECSNRFROMBER  Linear-on-log interpolation of the BER curve to find the
%   SNR at which BER == fecBer.  Returns NaN if the curve does not
%   bracket the FEC threshold.
    v = berVec(:).';
    v(v <= 0) = NaN;
    valid = isfinite(v);
    if nnz(valid) < 2
        snr = NaN;
        return;
    end
    x = snrVec(valid);
    y = log10(v(valid));
    target = log10(fecBer);
    if target < min(y) || target > max(y)
        snr = NaN;
        return;
    end
    for k = 1:numel(x) - 1
        if (y(k) - target) * (y(k+1) - target) <= 0
            snr = x(k) + (target - y(k)) / (y(k+1) - y(k)) * ...
                (x(k+1) - x(k));
            return;
        end
    end
    snr = NaN;
end
