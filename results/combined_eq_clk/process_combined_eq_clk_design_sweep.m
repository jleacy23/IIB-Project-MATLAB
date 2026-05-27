function process_combined_eq_clk_design_sweep(varargin)
%PROCESS_COMBINED_EQ_CLK_DESIGN_SWEEP  Print the design-sweep FEC SNR table.
%
%   process_combined_eq_clk_design_sweep()
%   process_combined_eq_clk_design_sweep('MatFile', path, 'FECBER', 2e-2)
%
%   Loads combined_eq_clk_design_sweep.mat and prints one row per
%   (block, NCD, NTaps, po2, ki, kp) configuration with the FEC SNR
%   mean and standard deviation across Monte-Carlo trials.
%
%   Name/Value options:
%     'MatFile' - path to the .mat (default: alongside this script)
%     'FECBER'  - FEC threshold used to score designs (default 2e-2)

    here = fileparts(mfilename('fullpath'));
    p = inputParser;
    p.addParameter('MatFile', ...
        fullfile(here, 'combined_eq_clk_design_sweep.mat'));
    p.addParameter('FECBER', 2e-2);
    p.parse(varargin{:});
    matFile = p.Results.MatFile;
    fecBer  = p.Results.FECBER;

    if ~isfile(matFile)
        error('process_combined_eq_clk_design_sweep:missingMat', ...
              ['Could not find %s.\n', ...
               'Run runtests(''combined_eq_clk_sweep/', ...
               'test_design_sweep'') first.'], matFile);
    end
    S      = load(matFile);
    tbl    = S.tbl;
    SNR_dB = S.SNR_dB_vec(:).';

    tbl = sortrows(tbl, {'block_name', 'n_cd', 'n_aeq', 'po2'});

    %% --- Per-trial FEC SNR statistics -----------------------------
    nRow = height(tbl);
    fecMean = nan(nRow, 1);
    fecStd  = nan(nRow, 1);
    for k = 1:nRow
        berMat = tbl.ber{k};               % [NTrials x NSNR]
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

    %% --- Print table ----------------------------------------------
    fprintf('\n--- FEC SNR per configuration (FEC BER = %.0e) ---\n', ...
        fecBer);
    fprintf('%-18s %-5s %-7s %-5s %-9s %-9s %-10s %-10s\n', ...
        'block', 'NCD', 'NTaps', 'po2', 'ki', 'kp', ...
        'FEC SNR', 'std');
    for k = 1:nRow
        meanStr = '   ---';
        stdStr  = '   ---';
        if isfinite(fecMean(k))
            meanStr = sprintf('%7.2f', fecMean(k));
        end
        if isfinite(fecStd(k))
            stdStr = sprintf('%7.2f', fecStd(k));
        end
        fprintf('%-18s %-5d %-7d %-5d %-9.1e %-9.1e %-10s %-10s\n', ...
            char(tbl.block_name(k)), tbl.n_cd(k), tbl.n_aeq(k), ...
            tbl.po2(k), tbl.ki(k), tbl.kp(k), meanStr, stdStr);
    end
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
