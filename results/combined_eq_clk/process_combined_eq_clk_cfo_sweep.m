function process_combined_eq_clk_cfo_sweep(varargin)
%PROCESS_COMBINED_EQ_CLK_CFO_SWEEP  Plot FEC SNR vs CFO and report 3 GHz table.
%
%   process_combined_eq_clk_cfo_sweep()
%   process_combined_eq_clk_cfo_sweep('MatFile', path, ...
%       'FECBER', 2e-2, 'CfoReport', 3, 'SavePlots', true)
%
%   Loads combined_eq_clk_cfo_sweep.mat and produces:
%     1. A single plot of FEC SNR vs CFO, one line per implementation
%        (block, NCD, NTaps, po2, ki, kp).
%     2. A console table of FEC SNR (mean +/- std over trials) at the
%        worst-case CFO (default 3 GHz).
%
%   Name/Value options:
%     'MatFile'   - path to the .mat (default: alongside this script)
%     'FECBER'    - FEC threshold (default 2e-2)
%     'CfoReport' - CFO (GHz) at which to print the summary table
%                   (default 3).  Must be present in CFO_GHz_vec.
%     'SavePlots' - true to write a PNG alongside the .mat (default true)

    here = fileparts(mfilename('fullpath'));
    p = inputParser;
    p.addParameter('MatFile', ...
        fullfile(here, 'combined_eq_clk_cfo_sweep.mat'));
    p.addParameter('FECBER',    2e-2);
    p.addParameter('CfoReport', 3);
    p.addParameter('SavePlots', true);
    p.parse(varargin{:});
    matFile   = p.Results.MatFile;
    fecBer    = p.Results.FECBER;
    cfoReport = p.Results.CfoReport;
    savePlots = p.Results.SavePlots;

    if ~isfile(matFile)
        error('process_combined_eq_clk_cfo_sweep:missingMat', ...
              ['Could not find %s.\n', ...
               'Run runtests(''combined_eq_clk_sweep/', ...
               'test_cfo_sweep'') first.'], matFile);
    end
    S       = load(matFile);
    tbl     = S.tbl;
    SNR_dB  = S.SNR_dB_vec(:).';
    CFO_vec = S.CFO_GHz_vec(:).';
    NCFO    = numel(CFO_vec);
    NCFG    = height(tbl);

    %% --- Per-trial FEC SNR (mean +/- std across trials) ------------
    fecSnrMean = nan(NCFG, NCFO);
    fecSnrStd  = nan(NCFG, NCFO);
    for ci = 1:NCFG
        berCube = tbl.ber{ci};                  % [NTrials x NSNR x NCFO]
        nTrials = size(berCube, 1);
        for ic = 1:NCFO
            berSlice = squeeze(berCube(:, :, ic));   % [NTrials x NSNR]
            perTrial = nan(nTrials, 1);
            for tr = 1:nTrials
                perTrial(tr) = fecSnrFromBer( ...
                    SNR_dB, berSlice(tr, :), fecBer);
            end
            valid = isfinite(perTrial);
            if any(valid)
                fecSnrMean(ci, ic) = mean(perTrial(valid));
                fecSnrStd(ci, ic)  = std(perTrial(valid));
            end
        end
    end

    figure('Name', 'CFO sweep: FEC SNR vs CFO', ...
        'Position', [80 80 720 480]);
    ax = axes; hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    colors = lines(NCFG);
    for ci = 1:NCFG
        errorbar(ax, CFO_vec, fecSnrMean(ci, :), fecSnrStd(ci, :), ...
            '-o', 'Color', colors(ci, :), 'LineWidth', 1.4, ...
            'MarkerSize', 5, 'CapSize', 8, ...
            'DisplayName', cfgLabel(tbl, ci));
    end
    xlabel(ax, 'CFO (GHz)');
    ylabel(ax, 'FEC SNR (dB)');
    title(ax, sprintf('FEC SNR vs CFO at BER = %.0e (L = %d km)', ...
        fecBer, tbl.l_km(1)));
    legend(ax, 'Location', 'eastoutside', 'Interpreter', 'none');

    if savePlots
        outFile = fullfile(here, 'combined_eq_clk_cfo_fec_snr.png');
        exportgraphics(gcf, outFile, 'Resolution', 200);
        fprintf('Saved FEC SNR vs CFO plot to %s\n', outFile);
    end

    %% --- Table of FEC SNR (mean +/- std) at cfoReport --------------
    ic = find(abs(CFO_vec - cfoReport) < 1e-9, 1);
    if isempty(ic)
        warning('process_combined_eq_clk_cfo_sweep:cfoMissing', ...
            ['CfoReport = %.3f GHz not in CFO_GHz_vec [', ...
             repmat('%g ', 1, NCFO), ']; skipping table.'], ...
            cfoReport, CFO_vec);
        return;
    end

    fprintf('\n--- FEC SNR at CFO = %g GHz (FEC BER = %.0e) ---\n', ...
        cfoReport, fecBer);
    fprintf('%-18s %-5s %-7s %-5s %-9s %-9s %-10s %-10s\n', ...
        'block', 'NCD', 'NTaps', 'po2', 'ki', 'kp', ...
        'FEC SNR', 'std');
    for ci = 1:NCFG
        meanStr = '   ---';
        stdStr  = '   ---';
        if isfinite(fecSnrMean(ci, ic))
            meanStr = sprintf('%7.2f', fecSnrMean(ci, ic));
        end
        if isfinite(fecSnrStd(ci, ic))
            stdStr = sprintf('%7.2f', fecSnrStd(ci, ic));
        end
        fprintf('%-18s %-5d %-7d %-5d %-9.1e %-9.1e %-10s %-10s\n', ...
            char(tbl.block_name(ci)), tbl.n_cd(ci), tbl.n_aeq(ci), ...
            tbl.po2(ci), tbl.ki(ci), tbl.kp(ci), meanStr, stdStr);
    end
end


% =====================================================================
function s = cfgLabel(tbl, ci)
    s = sprintf('%s, NCD=%d, NTaps=%d, po2=%d', ...
        char(tbl.block_name(ci)), tbl.n_cd(ci), tbl.n_aeq(ci), ...
        tbl.po2(ci));
end


% =====================================================================
function snr = fecSnrFromBer(snrVec, berVec, fecBer)
%FECSNRFROMBER  Linear-on-log interp of the BER curve for FEC SNR.
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
