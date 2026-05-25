function process_combined_eq_sweep(varargin)
%PROCESS_COMBINED_EQ_SWEEP  Post-process the combined CD + AEQ sweep.
%
%   process_combined_eq_sweep()
%   process_combined_eq_sweep('FECBER', 2e-2, 'MatFile', ...)
%
%   Loads combined_eq_sweep.mat and prints, for each network, a table of
%   FEC SNR (mean +/- std across trials) for every swept combination of
%   (cd_config, aeq_mode, aeq_sign_only, aeq_ntaps, cd_fl, aeq_fl).
%   FEC SNR per trial is found by linear interpolation of the BER curve
%   at the FEC threshold.
%
%   Name/Value options
%     'MatFile' - path to the .mat (default: combined_eq_sweep.mat in
%                 this file's folder)
%     'FECBER'  - FEC threshold (default 2e-2)

    %% --- Parse args -----------------------------------------------------
    here = fileparts(mfilename('fullpath'));
    p = inputParser;
    p.addParameter('MatFile', fullfile(here, 'combined_eq_sweep.mat'));
    p.addParameter('FECBER',  2e-2);
    p.parse(varargin{:});
    matFile = p.Results.MatFile;
    fecBer  = p.Results.FECBER;

    %% --- Load ----------------------------------------------------------
    if ~isfile(matFile)
        error('process_combined_eq_sweep:missingMat', ...
              ['Could not find %s.\n', ...
               'Run runtests(''combined_eq_sweep'') first to produce it.'], ...
              matFile);
    end
    S      = load(matFile);
    tbl    = S.tbl;
    SNR_dB = S.SNR_dB_vec;
    params = S.params;

    %% --- FEC SNR table per network -------------------------------------
    nets = unique(tbl.network, 'stable');
    for ni = 1:numel(nets)
        netName = nets(ni);
        subN    = tbl(tbl.network == netName, :);
        if isempty(subN), continue; end

        nR = height(subN);
        CD_Cfg     = strings(nR, 1);
        AEQ_Mode   = strings(nR, 1);
        Variant    = strings(nR, 1);
        NTaps      = zeros(nR, 1);
        CD_FL      = zeros(nR, 1);
        AEQ_FL     = zeros(nR, 1);
        FEC_SNR_dB = nan(nR, 1);
        FEC_Std_dB = nan(nR, 1);

        for k = 1:nR
            CD_Cfg(k)   = subN.cd_config(k);
            AEQ_Mode(k) = subN.aeq_mode(k);
            Variant(k)  = ternaryStr(subN.aeq_sign_only(k), ...
                                     'sign-sign', 'direct');
            NTaps(k)    = subN.aeq_ntaps(k);
            CD_FL(k)    = subN.cd_fl(k);
            AEQ_FL(k)   = subN.aeq_fl(k);

            trialSNR      = fecCrossingsAll(SNR_dB, subN.ber{k}, fecBer);
            FEC_SNR_dB(k) = mean(trialSNR, 'omitnan');
            FEC_Std_dB(k) = std (trialSNR, 'omitnan');
        end

        T = table(CD_Cfg, AEQ_Mode, Variant, NTaps, CD_FL, AEQ_FL, ...
                  FEC_SNR_dB, FEC_Std_dB);
        T = sortrows(T, {'CD_Cfg', 'AEQ_Mode', 'Variant', ...
                         'NTaps', 'CD_FL', 'AEQ_FL'});

        L_km_  = subN.L_km(1);
        split_ = subN.splitting(1);
        fprintf(['\n=== Combined CD + AEQ FEC SNR  Net %s  L = %d km  ', ...
                 'split %s  (FEC BER = %.0e, mean +/- std over %d trials) ===\n'], ...
                netName, L_km_, split_, fecBer, params.NTrials);
        disp(T);
    end

    %% --- Plot: FEC SNR vs AEQ FL per (network, CD config) --------------
    %  One figure per (network, cd_config) pair, one line per AEQ NTaps.
    %  Sign-sign CMA only (the only AEQ variant the sweep retains).
    cdCfgs   = unique(tbl.cd_config, 'stable');
    ntapsVec = sort(unique(tbl.aeq_ntaps));
    flAxis   = sort(unique(tbl.aeq_fl));
    markers  = {'o', 's', '^', 'd'};
    figIdx   = 0;

    for ni = 1:numel(nets)
        netName = nets(ni);
        subN    = tbl(tbl.network == netName, :);
        if isempty(subN), continue; end
        L_km_  = subN.L_km(1);
        split_ = subN.splitting(1);

        for ci = 1:numel(cdCfgs)
            cdName = cdCfgs(ci);
            subC   = subN(subN.cd_config == cdName, :);
            if isempty(subC), continue; end

            figIdx = figIdx + 1;
            figure('Name', sprintf('Combined FEC SNR vs AEQ FL  Net %s  %s', ...
                                    netName, cdName), ...
                   'Position', [260 + (figIdx-1)*60, 220, 760, 560]);
            hold on; grid on;
            cmap = lines(numel(ntapsVec));

            for ti = 1:numel(ntapsVec)
                nt   = ntapsVec(ti);
                rows = subC(subC.aeq_ntaps == nt, :);
                if isempty(rows), continue; end
                rows = sortrows(rows, 'aeq_fl');

                flVals  = rows.aeq_fl;
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

            xlabel('AEQ gradient FL (fractional bits)');
            ylabel(sprintf('FEC SNR @ BER = %.0e   [dB]', fecBer));
            title(sprintf(['%s   Net %s   L = %d km   split %s   |   ', ...
                           'FEC SNR vs AEQ FL   |   mean \\pm std over %d trials'], ...
                          cdName, netName, L_km_, split_, params.NTrials));
            legend('show', 'Location', 'best');
            xticks(flAxis);
        end
    end
end


%% =====================================================================
%  Local helpers
%  =====================================================================

function s = ternaryStr(cond, ifTrue, ifFalse)
    if cond, s = ifTrue; else, s = ifFalse; end
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
