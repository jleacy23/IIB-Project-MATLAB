function energy_snr_results()
%ENERGY_SNR_RESULTS  Plot carrier-recovery RSNR vs fixed-point precision.
%
%   Re-plots the precision-sweep results saved by bit_width_full
%   (bit_width_full_precision_sweep.mat) without re-running the sweep.  Both
%   sweeps are drawn on a single graph as RSNR (mean +/- std) vs fractional
%   bits:
%       * Pilots-only phase precision  (differential Kay held at high precision)
%       * Differential-Kay frequency precision (pilots-only held at high precision)
%
%   The .mat is produced by:  runtests('bit_width_full')
%   It stores: FL_vec, rsnr_cr_mean, rsnr_cr_std, rsnr_fr_mean, rsnr_fr_std.
%
%   (Energy-vs-precision analysis now lives in precision_energy.m.)
%
%   Run:  energy_snr_results

    here    = fileparts(mfilename('fullpath'));
    matFile = fullfile(here, 'bit_width_full_precision_sweep.mat');
    if ~isfile(matFile)
        error('energy_snr_results:noData', ...
            ['No sweep results found at\n  %s\n' ...
             'Run  runtests(''bit_width_full'')  first.'], matFile);
    end

    d  = load(matFile);
    fl = d.FL_vec(:);

    fig = figure('Name', 'CR RSNR vs precision', 'Color', 'w', ...
        'Position', [300, 300, 820, 520]);
    hold on; grid on;

    errorbar(fl, d.rsnr_cr_mean(:), d.rsnr_cr_std(:), '-o', ...
        'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'auto', 'CapSize', 6, ...
        'DisplayName', 'Pilots-only precision (differential Kay at high precision)');

    % Drop the 2-bit point from the differential-Kay (frequency) sweep.
    keepFR = fl ~= 2;
    errorbar(fl(keepFR), d.rsnr_fr_mean(keepFR), d.rsnr_fr_std(keepFR), '-s', ...
        'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'auto', 'CapSize', 6, ...
        'DisplayName', 'Differential Kay precision (pilots-only at high precision)');

    xlabel('Fractional bits', 'FontSize', 11);
    ylabel('RSNR [dB]', 'FontSize', 11);
    title('Carrier recovery RSNR vs fixed-point precision', 'FontSize', 11);
    legend('Location', 'best', 'Interpreter', 'none', 'FontSize', 11);
    xlim([min(fl) - 1, max(fl) + 1]);

    exportgraphics(fig, fullfile(here, 'precision_vs_rsnr.png'));
    fprintf('Saved combined precision-vs-RSNR plot to %s\n', ...
        fullfile(here, 'precision_vs_rsnr.png'));
end
