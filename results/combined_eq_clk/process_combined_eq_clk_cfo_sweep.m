function process_combined_eq_clk_cfo_sweep(varargin)
%PROCESS_COMBINED_EQ_CLK_CFO_SWEEP  Plot the CFO-sweep results.
%
%   process_combined_eq_clk_cfo_sweep()
%   process_combined_eq_clk_cfo_sweep('MatFile', path, ...
%       'FECBER', 2e-2, 'SavePlots', true)
%
%   Loads combined_eq_clk_cfo_sweep.mat (one design per block, swept
%   over (SNR, CFO)) and produces:
%     1. BER vs SNR curves, one panel per block, one line per CFO.
%     2. FEC SNR vs CFO, one line per block.
%
%   Name/Value options:
%     'MatFile'   - path to the .mat (default: alongside this script)
%     'FECBER'    - FEC threshold (default 2e-2)
%     'SavePlots' - true to write PNGs alongside the .mat (default true)

    here = fileparts(mfilename('fullpath'));
    p = inputParser;
    p.addParameter('MatFile', ...
        fullfile(here, 'combined_eq_clk_cfo_sweep.mat'));
    p.addParameter('FECBER',    2e-2);
    p.addParameter('SavePlots', true);
    p.parse(varargin{:});
    matFile   = p.Results.MatFile;
    fecBer    = p.Results.FECBER;
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

    %% --- BER vs SNR (all blocks x CFOs on a single panel) ----------
    figure('Name', 'CFO sweep: BER vs SNR', ...
        'Position', [80 80 720 480]);
    ax = axes; hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    colors = lines(NCFG * NCFO);
    lineStyles = {'-', '--', ':', '-.'};
    iLine = 0;
    for ci = 1:NCFG
        berCube = tbl.ber{ci};
        ls = lineStyles{mod(ci - 1, numel(lineStyles)) + 1};
        for ic = 1:NCFO
            iLine = iLine + 1;
            berSlice = squeeze(berCube(:, :, ic));
            meanBer  = mean(berSlice, 1, 'omitnan');
            minBer   = min(berSlice, [], 1, 'omitnan');
            maxBer   = max(berSlice, [], 1, 'omitnan');

            meanBer(meanBer < 1e-6) = 1e-6;
            minBer(minBer   < 1e-6) = 1e-6;
            maxBer(maxBer   < 1e-6) = 1e-6;

            valid = isfinite(meanBer);
            col   = colors(iLine, :);
            fill(ax, [SNR_dB(valid), fliplr(SNR_dB(valid))], ...
                     [minBer(valid),  fliplr(maxBer(valid))], col, ...
                     'FaceAlpha', 0.08, 'EdgeColor', 'none', ...
                     'HandleVisibility', 'off');
            plot(ax, SNR_dB(valid), meanBer(valid), ls, ...
                 'Marker', 'o', 'Color', col, 'LineWidth', 1.3, ...
                 'MarkerSize', 4, ...
                 'DisplayName', sprintf('%s, CFO=%g GHz', ...
                    char(tbl.block_name(ci)), CFO_vec(ic)));
        end
    end
    yline(ax, fecBer, '--', sprintf('FEC %.0e', fecBer), ...
        'LabelHorizontalAlignment', 'left');
    set(ax, 'YScale', 'log');
    xlabel(ax, 'SNR (dB)');
    ylabel(ax, 'BER');
    title(ax, sprintf('CFO sweep BER vs SNR (L = %d km)', tbl.l_km(1)));
    legend(ax, 'Location', 'eastoutside', 'Interpreter', 'none');
    ylim(ax, [1e-6, 0.5]);

    if savePlots
        outFile = fullfile(here, 'combined_eq_clk_cfo_ber_vs_snr.png');
        exportgraphics(gcf, outFile, 'Resolution', 200);
        fprintf('Saved BER vs SNR plot to %s\n', outFile);
    end

    %% --- FEC SNR vs CFO, one line per block ------------------------
    figure('Name', 'CFO sweep: FEC SNR vs CFO', ...
        'Position', [80 80 600 400]);
    ax = axes; hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    blockColors = lines(NCFG);
    fecSnrMat   = nan(NCFG, NCFO);
    for ci = 1:NCFG
        berCube = tbl.ber{ci};
        for ic = 1:NCFO
            berSlice = squeeze(berCube(:, :, ic));
            meanBer  = mean(berSlice, 1, 'omitnan');
            fecSnrMat(ci, ic) = fecSnrFromBer(SNR_dB, meanBer, fecBer);
        end
        plot(ax, CFO_vec, fecSnrMat(ci, :), '-o', ...
            'Color', blockColors(ci, :), 'LineWidth', 1.4, ...
            'MarkerSize', 5, ...
            'DisplayName', char(tbl.block_name(ci)));
    end
    xlabel(ax, 'CFO (GHz)');
    ylabel(ax, 'FEC SNR (dB)');
    title(ax, sprintf('FEC SNR vs CFO at BER = %.0e', fecBer));
    legend(ax, 'Location', 'best', 'Interpreter', 'none');

    if savePlots
        outFile = fullfile(here, 'combined_eq_clk_cfo_fec_snr.png');
        exportgraphics(gcf, outFile, 'Resolution', 200);
        fprintf('Saved FEC SNR vs CFO plot to %s\n', outFile);
    end

    %% --- Console summary -------------------------------------------
    for ic = 1:NCFO
        fprintf('\n--- Mean BER vs SNR  (CFO = %g GHz) ---\n', CFO_vec(ic));
        header = sprintf('%-26s %-5s %-7s', 'block', 'NCD', 'NTaps');
        for si = 1:numel(SNR_dB)
            header = [header, sprintf(' %7.1f', SNR_dB(si))]; %#ok<AGROW>
        end
        fprintf('%s\n', header);
        for k = 1:NCFG
            row = sprintf('%-26s %-5d %-7d', ...
                char(tbl.block_name(k)), tbl.n_cd(k), tbl.n_aeq(k));
            m = mean(tbl.ber{k}(:, :, ic), 1, 'omitnan');
            for si = 1:numel(SNR_dB)
                row = [row, sprintf(' %7.1e', m(si))]; %#ok<AGROW>
            end
            fprintf('%s\n', row);
        end
    end

    fprintf('\n--- FEC SNR (dB) vs CFO ---\n');
    header = sprintf('%-26s', 'block');
    for ic = 1:NCFO
        header = [header, sprintf(' %7.2f', CFO_vec(ic))]; %#ok<AGROW>
    end
    fprintf('%s\n', header);
    for k = 1:NCFG
        row = sprintf('%-26s', char(tbl.block_name(k)));
        for ic = 1:NCFO
            v = fecSnrMat(k, ic);
            if isfinite(v)
                row = [row, sprintf(' %7.2f', v)]; %#ok<AGROW>
            else
                row = [row, sprintf(' %7s', '---')]; %#ok<AGROW>
            end
        end
        fprintf('%s\n', row);
    end
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
