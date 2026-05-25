function process_adaptive_eq_grid_sweep(varargin)
%PROCESS_ADAPTIVE_EQ_GRID_SWEEP  Post-process the adaptive-EQ grid sweep.
%
%   process_adaptive_eq_grid_sweep()
%   process_adaptive_eq_grid_sweep('FL', 12, 'FECBER', 2e-2, 'MatFile', ...)
%
%   Loads adaptive_eq_grid_sweep.mat and prints, for each mode
%   (CMA, pilot) and each network, a table of FEC SNR (mean +/- std
%   across trials) at the requested gradient precision (default FL = 12)
%   for every combination of update variant (direct / sign-sign) and
%   filter length (NTaps).
%
%   Name/Value options
%     'MatFile'  - path to the .mat (default: adaptive_eq_grid_sweep.mat
%                  in this file's folder)
%     'FL'       - fractional bit width to report (default 12)
%     'FECBER'   - FEC threshold (default 2e-2)

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

    %% --- Variant set ---------------------------------------------------
    if hasSignOnly(tbl)
        variants = [false, true];
    else
        variants = false;
        warning('process_adaptive_eq_grid_sweep:noSignOnly', ...
                ['adaptive_eq_grid_sweep.mat has no sign_only column; ', ...
                 'reporting direct-update rows only.']);
    end

    %% --- FEC SNR table per (mode, network) @ FL = flTarget -------------
    modes = unique(tbl.mode, 'stable');
    for mi = 1:numel(modes)
        modeName = modes(mi);
        sel = (tbl.mode == modeName) & (tbl.fl == flTarget);
        sub = tbl(sel, :);
        if isempty(sub)
            warning('process_adaptive_eq_grid_sweep:noRows', ...
                    'No %s rows at FL = %d in %s.', ...
                    modeName, flTarget, matFile);
            continue;
        end

        nets = unique(sub.network, 'stable');
        for ni = 1:numel(nets)
            netName = nets(ni);
            subN    = sub(sub.network == netName, :);
            L_km_   = subN.L_km(1);
            split_  = subN.splitting(1);

            ntapsVec = sort(unique(subN.ntaps));
            nRows    = numel(ntapsVec) * numel(variants);

            Variant    = strings(nRows, 1);
            NTaps      = zeros(nRows, 1);
            FEC_SNR_dB = nan(nRows, 1);
            FEC_Std_dB = nan(nRows, 1);

            k = 0;
            for vi = 1:numel(variants)
                signFlag = variants(vi);
                vTag     = ternaryStr(signFlag, 'sign-sign', 'direct');
                for ti = 1:numel(ntapsVec)
                    k = k + 1;
                    Variant(k) = string(vTag);
                    NTaps(k)   = ntapsVec(ti);

                    rows = subN(subN.ntaps == ntapsVec(ti), :);
                    if hasSignOnly(tbl)
                        rows = rows(rows.sign_only == signFlag, :);
                    end
                    if isempty(rows), continue; end

                    berMat   = rows.ber{1};
                    trialSNR = fecCrossingsAll(SNR_dB, berMat, fecBer);
                    FEC_SNR_dB(k) = mean(trialSNR, 'omitnan');
                    FEC_Std_dB(k) = std (trialSNR, 'omitnan');
                end
            end

            T = table(Variant, NTaps, FEC_SNR_dB, FEC_Std_dB);
            fprintf(['\n=== %s FEC SNR @ FL = %d  Net %s  L = %d km  ', ...
                     'split %s  (FEC BER = %.0e, mean +/- std over %d trials) ===\n'], ...
                    modeName, flTarget, netName, L_km_, split_, ...
                    fecBer, params.NTrials);
            disp(T);
        end
    end

    %% --- Plot: CMA (sign-sign) FEC SNR vs FL at NTaps in {1, 3} --------
    %  One line per (network, NTaps): mean across trials with std error
    %  bars.  Requires the sweep to have logged the sign_only column.
    if ~hasSignOnly(tbl)
        warning('process_adaptive_eq_grid_sweep:noSignOnlyPlot', ...
                ['adaptive_eq_grid_sweep.mat has no sign_only column; ', ...
                 'skipping CMA sign-sign FEC SNR vs FL plot.']);
        return;
    end
    cmaTaps = [1, 3];
    selC = (tbl.mode == "CMA") & (tbl.sign_only == true) & ...
           ismember(tbl.ntaps, cmaTaps);
    subC = tbl(selC, :);
    if isempty(subC)
        warning('process_adaptive_eq_grid_sweep:noCMASignRows', ...
                'No CMA sign-sign rows with NTaps in [%s] in %s.', ...
                num2str(cmaTaps), matFile);
        return;
    end

    netsC   = unique(subC.network, 'stable');
    markers = {'o', 's', '^', 'd'};
    flAxis  = sort(unique(subC.fl));

    for ni = 1:numel(netsC)
        netName = netsC(ni);
        subN    = subC(subC.network == netName, :);
        if isempty(subN), continue; end
        L_km_  = subN.L_km(1);
        split_ = subN.splitting(1);

        figure('Name', sprintf('CMA (sign-sign) FEC SNR vs FL  Net %s', netName), ...
               'Position', [260 + (ni-1)*60, 220, 760, 560]);
        hold on; grid on;
        cmap = lines(numel(cmaTaps));

        for ti = 1:numel(cmaTaps)
            nt   = cmaTaps(ti);
            rows = subN(subN.ntaps == nt, :);
            if isempty(rows), continue; end
            rows = sortrows(rows, 'fl');

            flVals  = rows.fl;
            meanFEC = nan(size(flVals));
            stdFEC  = nan(size(flVals));
            for k = 1:numel(flVals)
                trialSNR   = fecCrossingsAll(SNR_dB, rows.ber{k}, fecBer);
                meanFEC(k) = mean(trialSNR, 'omitnan');
                stdFEC(k)  = std (trialSNR, 'omitnan');
            end

            errorbar(flVals, meanFEC, stdFEC, ['-' markers{mod(ti-1,4)+1}], ...
                     'Color', cmap(ti, :), ...
                     'MarkerSize', 6, 'LineWidth', 1.4, 'CapSize', 6, ...
                     'DisplayName', sprintf('N_{taps} = %d', nt));
        end

        xlabel('Gradient FL (fractional bits)');
        ylabel(sprintf('FEC SNR @ BER = %.0e   [dB]', fecBer));
        title(sprintf(['CMA (sign-sign)   Net %s   L = %d km   split %s   |   ', ...
                       'FEC SNR vs FL   |   mean \\pm std over %d trials'], ...
                      netName, L_km_, split_, params.NTrials));
        legend('show', 'Location', 'best');
        xticks(flAxis);
    end

    %% --- Energy + FEC SNR: CMA (sign-sign) vs FL ----------------------
    %  Per-symbol RM/RA for CMA sign-sign from report tab:adaptive_cost_cma
    %  divided by 2P (per the per-symbol convention in section 3.x):
    %    RM = (16N + 4)/2 = 8N + 2
    %    RA = (24N + 6)/2 = 12N + 3
    %  Energy is per bit via src/+energy/receiver.m (M = 4, Oversampling = 1,
    %  n = FL).  FEC SNR (mean +/- std) is per network from the BER cube.
    %  One table per network: rows = FL, columns = energy and FEC SNR for
    %  each NTaps value.
    M_qam   = 4;
    NTaps_e = [1, 3];
    RM_e    = 8 * NTaps_e + 2;
    RA_e    = 12 * NTaps_e + 3;

    FL         = flAxis(:);
    energyCols = nan(numel(FL), numel(NTaps_e));
    for ti = 1:numel(NTaps_e)
        energyCols(:, ti) = arrayfun(@(n) ...
            energy.receiver(RA_e(ti), RM_e(ti), EAdd, EMult, M_qam, 1, n), FL);
    end

    fprintf(['\n=== CMA (sign-sign) energy per bit and FEC SNR vs FL  ', ...
             '(E_A = %.3f fJ, E_M = %.3f fJ per bit) ===\n'], EAdd, EMult);
    fprintf('Per-symbol ops: ');
    for ti = 1:numel(NTaps_e)
        fprintf('N=%d -> RM=%d, RA=%d   ', NTaps_e(ti), RM_e(ti), RA_e(ti));
    end
    fprintf('\n');

    for ni = 1:numel(netsC)
        netName = netsC(ni);
        subN    = subC(subC.network == netName, :);
        if isempty(subN), continue; end
        L_km_  = subN.L_km(1);
        split_ = subN.splitting(1);

        fecMean = nan(numel(FL), numel(NTaps_e));
        fecStd  = nan(numel(FL), numel(NTaps_e));
        for ti = 1:numel(NTaps_e)
            rows = subN(subN.ntaps == NTaps_e(ti), :);
            for k = 1:numel(FL)
                rk = rows(rows.fl == FL(k), :);
                if isempty(rk), continue; end
                trialSNR     = fecCrossingsAll(SNR_dB, rk.ber{1}, fecBer);
                fecMean(k, ti) = mean(trialSNR, 'omitnan');
                fecStd (k, ti) = std (trialSNR, 'omitnan');
            end
        end

        % Interleave columns so each tap count's energy/FEC sits together.
        cols     = nan(numel(FL), 1 + 3 * numel(NTaps_e));
        colNames = cell(1, 1 + 3 * numel(NTaps_e));
        cols(:, 1)  = FL;
        colNames{1} = 'FL';
        for ti = 1:numel(NTaps_e)
            base = 1 + (ti-1)*3;
            cols(:, base + 1) = energyCols(:, ti);
            cols(:, base + 2) = fecMean(:, ti);
            cols(:, base + 3) = fecStd(:, ti);
            colNames{base + 1} = sprintf('Energy_N%d_fJ',     NTaps_e(ti));
            colNames{base + 2} = sprintf('FEC_SNR_N%d_dB',    NTaps_e(ti));
            colNames{base + 3} = sprintf('FEC_Std_N%d_dB',    NTaps_e(ti));
        end

        T = array2table(cols, 'VariableNames', colNames);
        fprintf(['\n--- Net %s  L = %d km  split %s  ', ...
                 '(FEC BER = %.0e, mean +/- std over %d trials) ---\n'], ...
                netName, L_km_, split_, fecBer, params.NTrials);
        disp(T);
    end
end


%% =====================================================================
%  Local helpers
%  =====================================================================

function s = ternaryStr(cond, ifTrue, ifFalse)
    if cond, s = ifTrue; else, s = ifFalse; end
end

function tf = hasSignOnly(tbl)
    tf = any(strcmp(tbl.Properties.VariableNames, 'sign_only'));
end

function snrs = fecCrossingsAll(snrDb, berPerTrial, fecBer)
    NTr  = size(berPerTrial, 1);
    snrs = nan(NTr, 1);
    for t = 1:NTr
        snrs(t) = fecCrossing(snrDb, berPerTrial(t, :), fecBer);
    end
end

function xCross = fecCrossing(x, y, yLimit)
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
