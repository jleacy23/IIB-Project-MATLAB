function process_cd_eq_precision_sweep(varargin)
%PROCESS_CD_EQ_PRECISION_SWEEP  Post-process the CD-equaliser precision sweep.
%
%   process_cd_eq_precision_sweep()
%   process_cd_eq_precision_sweep('FECBER', 2e-2, 'MatFile', ...)
%
%   Loads cd_eq_precision_sweep.mat (next to this file by default) and
%   produces two sets of plots:
%
%   Plot 1 — FEC SNR vs FL, one figure per network, one line per CD
%            equaliser implementation (time_domain, overlap_save,
%            overlap_save_po2).  FEC SNR is found by linear interpolation
%            of each per-trial BER curve at the FEC threshold; the marker
%            is the mean across trials, the error bar the std.
%   Plot 2 — Energy per bit vs FL, one figure per network, one line per
%            config.  Computed via src/+energy/receiver.m using the
%            per-symbol RM / RA counts from report tab:cd_cost_eval
%            (E_A = EAdd*n, E_M = EMult*n^2 per horowitz2014computing).
%
%   Name/Value options
%     'MatFile' - path to the .mat (default: cd_eq_precision_sweep.mat
%                 in this file's folder)
%     'FECBER'  - FEC threshold to mark (default 2e-2)
%     'EAdd'    - per-bit add-energy coefficient [fJ] (default 3.16)
%     'EMult'   - per-bit multiply-energy coefficient [fJ] (default 3.03)

    %% --- Parse args -----------------------------------------------------
    here = fileparts(mfilename('fullpath'));
    p = inputParser;
    p.addParameter('MatFile', fullfile(here, 'cd_eq_precision_sweep.mat'));
    p.addParameter('FECBER',  2e-2);
    p.addParameter('EAdd',    3.16/4);   % fJ per add per bit (horowitz2014)
    p.addParameter('EMult',   3.03/4);   % fJ per mult per bit (horowitz2014)
    p.parse(varargin{:});
    matFile = p.Results.MatFile;
    fecBer  = p.Results.FECBER;
    EAdd    = p.Results.EAdd;
    EMult   = p.Results.EMult;

    %% --- Load -----------------------------------------------------------
    if ~isfile(matFile)
        error('process_cd_eq_precision_sweep:missingMat', ...
              ['Could not find %s.\n', ...
               'Run runtests(''cd_eq_precision_sweep'') first to produce it.'], ...
              matFile);
    end
    S      = load(matFile);
    tbl    = S.tbl;
    SNR_dB = S.SNR_dB_vec;
    params = S.params;

    nets = unique(tbl.network, 'stable');

    %% =================================================================
    %  Plot 1: FEC SNR vs FL, one figure per network, line per config
    %  =================================================================
    %  For each per-trial BER curve, interpolate the SNR at which it
    %  crosses fecBer.  Marker = mean across trials, error bar = std.
    %  NaN trials (curve never crosses) are omitted from mean/std.
    cfgOrder = {'time_domain', 'overlap_save', 'overlap_save_po2'};
    for ni = 1:numel(nets)
        netName = nets(ni);
        subN = tbl(tbl.network == netName, :);
        if isempty(subN), continue; end

        L_km_  = subN.L_km(1);
        split_ = subN.splitting(1);

        figure('Name', sprintf('CD FEC SNR vs FL  Net %s', netName), ...
               'Position', [220, 200, 720, 540]);
        hold on; grid on;

        nCfg = numel(cfgOrder);
        cmap = lines(nCfg);
        for ci = 1:nCfg
            cfgName = cfgOrder{ci};
            rows = subN(subN.config == string(cfgName), :);
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

            labelStr = strrep(cfgName, '_', '\_');
            errorbar(flVals, meanFEC, stdFEC, ...
                     'o-', 'Color', cmap(ci, :), ...
                     'MarkerSize', 6, 'LineWidth', 1.4, 'CapSize', 6, ...
                     'DisplayName', labelStr);
        end

        xlabel('Gradient FL (fractional bits)');
        ylabel(sprintf('FEC SNR @ BER = %.0e   [dB]', fecBer));
        title(sprintf(['Net %s   L = %d km   split %s   |   ', ...
                       'FEC SNR vs FL   |   mean \\pm std over %d trials'], ...
                      netName, L_km_, split_, params.NTrials));
        legend('show', 'Location', 'best');
        xticks(sort(unique(tbl.fl)));
    end

    %% =================================================================
    %  Table: Energy/bit and FEC SNR for selected operating points
    %  =================================================================
    %  Selected operating points (minimum FL that meets FEC threshold):
    %    time_domain      @ FL = 8
    %    overlap_save     @ FL = 4
    %    overlap_save_po2 @ FL = 4
    selPoints = { ...
        'time_domain',      8; ...
        'overlap_save',     4; ...
        'overlap_save_po2', 4};

    M_qam = 4;
    colW  = 90;
    fprintf('\n%s\n', repmat('=', 1, colW));
    fprintf('%-12s  %-22s  %4s  %16s  %13s  %10s\n', ...
            'Network', 'Config', 'FL', 'Energy/bit (fJ)', 'FEC SNR (dB)', 'Std (dB)');
    fprintf('%s\n', repmat('-', 1, colW));

    for ni = 1:numel(nets)
        netName = nets(ni);
        subN    = tbl(tbl.network == netName, :);
        if isempty(subN), continue; end
        L_km_ = subN.L_km(1);

        for si = 1:size(selPoints, 1)
            cfgName = selPoints{si, 1};
            fl_sel  = selPoints{si, 2};

            [NM, NA] = reportOpsPerSymbol(cfgName, L_km_);
            if isnan(NM)
                E_bit = NaN;
            else
                E_bit = energy.receiver(NA, NM, EAdd, EMult, M_qam, 1, fl_sel);
            end

            rows = subN(subN.config == string(cfgName) & subN.fl == fl_sel, :);
            if isempty(rows)
                meanFEC = NaN;  stdFEC = NaN;
            else
                trialSNR = fecCrossingsAll(SNR_dB, rows.ber{1}, fecBer);
                meanFEC  = mean(trialSNR, 'omitnan');
                stdFEC   = std (trialSNR, 'omitnan');
            end

            fprintf('%-12s  %-22s  %4d  %16.3f  %13.3f  %10.3f\n', ...
                    netName, cfgName, fl_sel, E_bit, meanFEC, stdFEC);
        end
        fprintf('%s\n', repmat('-', 1, colW));
    end
    fprintf('%s\n', repmat('=', 1, colW));
end


%% =====================================================================
%  Local helpers
%  =====================================================================

function [NM, NA] = reportOpsPerSymbol(cfgName, L_km)
% REPORTOPSPERSYMBOL  Per-symbol RM/RA counts from report tab:cd_cost_eval.
%
%   Values from report/full/full.tex (evaluated with the N_CD of
%   tab:cd_taps at 2 samples per symbol):
%
%       L (km) | N   | TD RM | TD RA | OS RM | OS RA | OS-po2 RM | OS-po2 RA
%       20     | 32  |  48   |  44   |  56.9 |  75.9 |   9.48    |   75.9
%       80     | 128 | 168   | 164   |  75.9 | 104.3 |   9.48    |  104.3
    NM = NaN; NA = NaN;
    switch cfgName
        case 'time_domain'
            if     L_km == 20, NM =  48; NA =  44;
            elseif L_km == 80, NM = 168; NA = 164;
            end
        case 'overlap_save'
            if     L_km == 20, NM = 56.9; NA =  75.9;
            elseif L_km == 80, NM = 75.9; NA = 104.3;
            end
        case 'overlap_save_po2'
            if     L_km == 20, NM = 9.48; NA =  75.9;
            elseif L_km == 80, NM = 9.48; NA = 104.3;
            end
    end
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
