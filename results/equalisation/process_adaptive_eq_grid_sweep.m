function process_adaptive_eq_grid_sweep(varargin)
%PROCESS_ADAPTIVE_EQ_GRID_SWEEP  Post-process the adaptive-EQ grid sweep.
%
%   process_adaptive_eq_grid_sweep()
%   process_adaptive_eq_grid_sweep('FL', 16, 'FECBER', 2e-2, ...
%                                  'MatFile', ...)
%
%   Loads adaptive_eq_grid_sweep.mat (next to this file by default) and
%   produces two sets of plots:
%
%   Plot 1 — BER vs SNR for the CMA equaliser at the requested gradient
%            precision (default FL = 16), one figure per network, one
%            line per NTaps (direct-update CMA only).  Each line is the
%            mean BER across trials; BER = 0 is mapped to the BER floor.
%   Plot 2 — FEC SNR vs FL at NTaps = 1, one figure per network, with one
%            line for each of the 4 combinations of {CMA, pilot} x
%            {direct, sign-sign}.  FEC SNR per trial is found by linear
%            interpolation of the BER curve at the FEC threshold; marker
%            = mean across trials, error bar = std.  Requires the sweep
%            to have logged the 'sign_only' dimension; if absent only
%            the two direct-update lines are drawn.
%   Plot 3 — Energy per bit vs FL at NTaps = 1, one figure per network,
%            same 4 lines as Plot 2.  Per-symbol RM/RA counts are taken
%            from report tab:adaptive_cost_total; energy is computed via
%            src/+energy/receiver.m with E_A = EAdd*n, E_M = EMult*n^2
%            (horowitz2014computing), matching process_cd_eq_precision_sweep.
%
%   Name/Value options
%     'MatFile'  - path to the .mat (default: adaptive_eq_grid_sweep.mat
%                  in this file's folder)
%     'FL'       - fractional bit width for the Plot 1 BER curves (default 16)
%     'FECBER'   - FEC threshold to mark (default 2e-2)
%     'EAdd'     - per-bit add-energy coefficient [fJ] (default 3.16/4)
%     'EMult'    - per-bit multiply-energy coefficient [fJ] (default 3.03/4)

    %% --- Parse args -----------------------------------------------------
    here = fileparts(mfilename('fullpath'));
    p = inputParser;
    p.addParameter('MatFile', fullfile(here, 'adaptive_eq_grid_sweep.mat'));
    p.addParameter('FL',      12);
    p.addParameter('FECBER',  2e-2);
    p.addParameter('EAdd',    3.16/4);   % fJ per add per bit (horowitz2014)
    p.addParameter('EMult',   3.03/4);   % fJ per mult per bit (horowitz2014)
    p.parse(varargin{:});
    matFile  = p.Results.MatFile;
    flTarget = p.Results.FL;
    fecBer   = p.Results.FECBER;
    EAdd     = p.Results.EAdd;
    EMult    = p.Results.EMult;

    %% --- Load ----------------------------------------------------------
    if ~isfile(matFile)
        error('process_adaptive_eq_grid_sweep:missingMat', ...
              ['Could not find %s.\n', ...
               'Run runtests(''adaptive_eq_grid_sweep'') first to produce it.'], ...
              matFile);
    end
    S      = load(matFile);
    tbl    = S.tbl;
    SNR_dB = S.SNR_dB_vec;
    params = S.params;

    %% --- BER floor: 1 error per total_bits_per_trial -------------------
    %  computeBER discards NOut transient symbols, then counts errors across
    %  both polarisations at 2 bits/symbol.  Reproduce that count from the
    %  saved run parameters so the floor is the smallest non-zero BER any
    %  single trial could possibly resolve.
    SUBFRAME_SYMS    = 3712;
    nSymPerTrial     = params.NSub * SUBFRAME_SYMS;
    mPerPol          = max(nSymPerTrial - params.NOut, 1);
    totalBitsPerTr   = mPerPol * params.N_pol * 2;
    berFloor         = 1 / totalBitsPerTr;

    %% --- Plot 1: CMA (direct), FL = flTarget, line per NTaps -----------
    sel = (tbl.fl == flTarget) & (tbl.mode == "CMA");
    if hasSignOnly(tbl)
        sel = sel & (tbl.sign_only == false);
    end
    sub = tbl(sel, :);
    if isempty(sub)
        error('process_adaptive_eq_grid_sweep:noRows', ...
              'No CMA rows at FL = %d in %s.', flTarget, matFile);
    end
    nets = unique(sub.network, 'stable');

    for ni = 1:numel(nets)
        netName = nets(ni);
        subN    = sub(sub.network == netName, :);
        L_km_   = subN.L_km(1);
        split_  = subN.splitting(1);

        figure('Name', sprintf('CMA BER vs SNR  FL=%d  Net %s', ...
                               flTarget, netName), ...
               'Position', [100 + (ni-1)*60, 100, 720, 540]);
        hold on; grid on;
        set(gca, 'YScale', 'log');

        ntapsVec = sort(unique(subN.ntaps));
        cmap     = lines(numel(ntapsVec));

        for ti = 1:numel(ntapsVec)
            nt     = ntapsVec(ti);
            rowi   = find(subN.ntaps == nt, 1);
            berMat = subN.ber{rowi};   % [NTrials x NSNR]

            % Floor every per-trial BER first so trials with zero errors
            % still contribute a finite (floor-valued) point to the mean.
            berClamped = max(berMat, berFloor);

            meanBER = mean(berClamped, 1, 'omitnan');

            plot(SNR_dB, meanBER, 'o-', ...
                 'Color', cmap(ti, :), ...
                 'MarkerSize', 5, 'LineWidth', 1.4, ...
                 'DisplayName', sprintf('N_{taps} = %d', nt));
        end

        yline(fecBer, 'r--', 'LineWidth', 1.2, ...
              'DisplayName', sprintf('FEC limit (%.0e)', fecBer));
        yline(berFloor, 'k:', 'LineWidth', 1.0, ...
              'DisplayName', sprintf('BER floor (1 / %d bits)', totalBitsPerTr));

        xlabel('SNR (dB)');
        ylabel('BER');
        title(sprintf(['Net %s   L = %d km   split %s   |   CMA   |   ', ...
                       'FL = %d (gradient)   |   mean over %d trials'], ...
                      netName, L_km_, split_, flTarget, params.NTrials));
        ylim([berFloor / 3, 1]);
        xlim([min(SNR_dB), max(SNR_dB)]);
        legend('show', 'Location', 'southwest');
    end

    %% =================================================================
    %  Plot 2: FEC SNR vs FL at NTaps = 1, one figure per network,
    %          one line per combination of {CMA, pilot} x {direct, sign-sign}
    %  =================================================================
    %  Combinations: {mode, sign_only, display name, line spec}.
    %  If the .mat lacks the sign_only column (older sweep), only the two
    %  direct-update lines are plotted.
    ntapsTarget = 1;
    combos = { ...
        'CMA',   false, 'CMA, direct',      '-o';  ...
        'CMA',   true,  'CMA, sign-sign',   '--o'; ...
        'pilot', false, 'pilot, direct',    '-s';  ...
        'pilot', true,  'pilot, sign-sign', '--s'};
    haveSign = hasSignOnly(tbl);
    if ~haveSign
        % Drop sign-sign rows since the data isn't there.
        combos = combos(~cell2mat(combos(:,2)), :);
        warning('process_adaptive_eq_grid_sweep:noSignOnly', ...
                ['adaptive_eq_grid_sweep.mat has no sign_only column; ', ...
                 'plotting direct-update lines only.  Re-run the sweep ', ...
                 'to capture the sign-sign variant.']);
    end
    nCombo = size(combos, 1);

    netsAll = unique(tbl.network, 'stable');
    for ni = 1:numel(netsAll)
        netName = netsAll(ni);
        subN = tbl((tbl.network == netName) & (tbl.ntaps == ntapsTarget), :);
        if isempty(subN), continue; end

        L_km_  = subN.L_km(1);
        split_ = subN.splitting(1);

        figure('Name', sprintf('AEQ FEC SNR vs FL  Net %s  N=%d', ...
                               netName, ntapsTarget), ...
               'Position', [220 + (ni-1)*60, 200, 760, 560]);
        hold on; grid on;

        cmap = lines(nCombo);
        for ci = 1:nCombo
            modeName = combos{ci, 1};
            signFlag = combos{ci, 2};
            labelStr = combos{ci, 3};
            lineSpec = combos{ci, 4};

            rows = subN(subN.mode == string(modeName), :);
            if haveSign
                rows = rows(rows.sign_only == signFlag, :);
            end
            if isempty(rows), continue; end
            rows = sortrows(rows, 'fl');

            flVals  = rows.fl;
            meanFEC = nan(size(flVals));
            stdFEC  = nan(size(flVals));
            for k = 1:numel(flVals)
                berMat = rows.ber{k};   % [NTrials x NSNR]
                trialSNR = fecCrossingsAll(SNR_dB, berMat, fecBer);
                meanFEC(k) = mean(trialSNR, 'omitnan');
                stdFEC(k)  = std (trialSNR, 'omitnan');
            end

            errorbar(flVals, meanFEC, stdFEC, lineSpec, ...
                     'Color', cmap(ci, :), ...
                     'MarkerSize', 6, 'LineWidth', 1.4, 'CapSize', 6, ...
                     'DisplayName', labelStr);
        end

        xlabel('Gradient FL (fractional bits)');
        ylabel(sprintf('FEC SNR @ BER = %.0e   [dB]', fecBer));
        title(sprintf(['Net %s   L = %d km   split %s   |   ', ...
                       'N_{taps} = %d   |   FEC SNR vs FL   |   ', ...
                       'mean \\pm std over %d trials'], ...
                      netName, L_km_, split_, ntapsTarget, params.NTrials));
        legend('show', 'Location', 'best');
        xticks(sort(unique(tbl.fl)));
    end

    %% =================================================================
    %  Plot 3: Energy per bit vs FL at NTaps = 1, one figure per network,
    %          same 4 combos as Plot 2.
    %  =================================================================
    %  Per-symbol RM/RA counts come from report tab:adaptive_cost_total
    %  (N = 1 tap, P = 32 lanes).  Energy is computed via
    %  energy.receiver(NA, NM, EAdd, EMult, M=4, Oversampling=1, n=FL),
    %  which already accounts for both polarisations and converts to
    %  energy per bit.
    M_qam  = 4;
    flAxis = sort(unique(tbl.fl));
    for ni = 1:numel(netsAll)
        netName = netsAll(ni);
        subN = tbl((tbl.network == netName) & (tbl.ntaps == ntapsTarget), :);
        if isempty(subN), continue; end

        L_km_  = subN.L_km(1);
        split_ = subN.splitting(1);

        figure('Name', sprintf('AEQ Energy/bit vs FL  Net %s  N=%d', ...
                               netName, ntapsTarget), ...
               'Position', [340 + (ni-1)*60, 280, 760, 560]);
        hold on; grid on;
        set(gca, 'YScale', 'log');

        cmap = lines(nCombo);
        for ci = 1:nCombo
            modeName = combos{ci, 1};
            signFlag = combos{ci, 2};
            labelStr = combos{ci, 3};
            lineSpec = combos{ci, 4};

            [NM, NA] = reportOpsPerSymbol(modeName, signFlag);
            if isnan(NM), continue; end

            E_per_bit = arrayfun(@(n) ...
                energy.receiver(NA, NM, EAdd, EMult, M_qam, 1, n), flAxis);

            plot(flAxis, E_per_bit, lineSpec, ...
                 'Color', cmap(ci, :), ...
                 'MarkerSize', 6, 'LineWidth', 1.4, ...
                 'DisplayName', sprintf('%s  (RM=%.3g, RA=%.3g)', ...
                                        labelStr, NM, NA));
        end

        xlabel('Gradient FL (fractional bits)');
        ylabel('Energy per bit  [fJ]');
        title(sprintf(['Net %s   L = %d km   split %s   |   ', ...
                       'N_{taps} = %d   |   ', ...
                       'E\\_A = %.2f\\cdot n   E\\_M = %.2f\\cdot n^2 fJ'], ...
                      netName, L_km_, split_, ntapsTarget, EAdd, EMult));
        legend('show', 'Location', 'best');
        xticks(flAxis);
    end

    %% =================================================================
    %  Summary table: pilot/direct @ FL=2  vs  CMA/sign-sign @ FL=4
    %  =================================================================
    %  One MATLAB table per network printed to the command window, with
    %  energy/bit and FEC SNR (mean +/- std across trials) for the two
    %  recommended low-precision operating points.
    summaryPicks = { ...
        'pilot', false, 2, 'pilot, direct  @ FL=2'; ...
        'CMA',   true,  4, 'CMA, sign-sign @ FL=4'};
    if ~haveSign
        summaryPicks = summaryPicks(~cell2mat(summaryPicks(:,2)), :);
    end

    for ni = 1:numel(netsAll)
        netName = netsAll(ni);
        subN = tbl((tbl.network == netName) & (tbl.ntaps == ntapsTarget), :);
        if isempty(subN), continue; end

        nPick      = size(summaryPicks, 1);
        Combo      = strings(nPick, 1);
        FL         = zeros(nPick, 1);
        RM         = zeros(nPick, 1);
        RA         = zeros(nPick, 1);
        Energy_fJ  = nan(nPick, 1);
        FEC_SNR_dB = nan(nPick, 1);
        FEC_Std_dB = nan(nPick, 1);

        for pi = 1:nPick
            modeName = summaryPicks{pi, 1};
            signFlag = summaryPicks{pi, 2};
            flPick   = summaryPicks{pi, 3};
            label    = summaryPicks{pi, 4};

            rows = subN(subN.mode == string(modeName) & subN.fl == flPick, :);
            if haveSign
                rows = rows(rows.sign_only == signFlag, :);
            end
            Combo(pi) = string(label);
            FL(pi)    = flPick;
            [nm, na]  = reportOpsPerSymbol(modeName, signFlag);
            RM(pi)    = nm;
            RA(pi)    = na;
            if isempty(rows), continue; end

            Energy_fJ(pi)  = energy.receiver(na, nm, EAdd, EMult, M_qam, 1, flPick);
            trialSNR       = fecCrossingsAll(SNR_dB, rows.ber{1}, fecBer);
            FEC_SNR_dB(pi) = mean(trialSNR, 'omitnan');
            FEC_Std_dB(pi) = std (trialSNR, 'omitnan');
        end

        T = table(Combo, FL, RM, RA, Energy_fJ, FEC_SNR_dB, FEC_Std_dB);
        fprintf('\n=== Summary  Net %s  L = %d km  split %s  (N_taps = %d) ===\n', ...
                netName, subN.L_km(1), subN.splitting(1), ntapsTarget);
        disp(T);
    end
end


%% =====================================================================
%  Local helpers
%  =====================================================================

function [NM, NA] = reportOpsPerSymbol(modeName, signOnly)
% REPORTOPSPERSYMBOL  Per-symbol RM/RA counts from report
% tab:adaptive_cost_total (N = 1 tap, P = 32 lanes).
%
%       Combination       | RM    | RA
%       CMA, direct       | 20    | 15
%       CMA, sign-sign    | 10    | 15
%       Pilot, direct     |  8.25 |  6.3125
%       Pilot, sign-sign  |  8    |  6.3125
    NM = NaN; NA = NaN;
    switch lower(string(modeName))
        case "cma"
            if signOnly, NM = 10;   NA = 15;
            else,        NM = 20;   NA = 15;
            end
        case "pilot"
            if signOnly, NM =  8;    NA = 6.3125;
            else,        NM =  8.25; NA = 6.3125;
            end
    end
end

function tf = hasSignOnly(tbl)
% HASSIGNONLY  True if the results table carries the sign_only column.
    tf = any(strcmp(tbl.Properties.VariableNames, 'sign_only'));
end

function snrs = fecCrossingsAll(snrDb, berPerTrial, fecBer)
% FECCROSSINGSALL  Per-trial FEC SNR crossings (linear interpolation).
    NTr  = size(berPerTrial, 1);
    snrs = nan(NTr, 1);
    for t = 1:NTr
        snrs(t) = fecCrossing(snrDb, berPerTrial(t, :), fecBer);
    end
end

function xCross = fecCrossing(x, y, yLimit)
% FECCROSSING  SNR at which BER curve y crosses yLimit (linear interp).
    xCross = NaN;
    n = numel(x);
    if n < 2, return; end
    exact = find(y == yLimit, 1, 'first');
    if ~isempty(exact)
        xCross = x(exact);
        return;
    end
    for i = 1:(n-1)
        if (y(i) - yLimit) * (y(i+1) - yLimit) < 0
            xCross = x(i) + (yLimit - y(i)) * (x(i+1) - x(i)) / ...
                     (y(i+1) - y(i));
            return;
        end
    end
end
