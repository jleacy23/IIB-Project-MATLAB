function tbl = energy_snr_results(ignore_fr_energy)
% ENERGY_SNR_RESULTS  Receiver DSP energy and FEC SNR analysis.
%
%   For each FR+CR algorithm pair, returns a summary table reporting the
%   precision combo with the lowest mean FEC SNR (best DSP performance)
%   and its associated Rx energy.  Also produces a CR-precision
%   comparison plot for pilots-only vs Viterbi-Viterbi with FR fixed to
%   R&B blind.
%
%   ignore_fr_energy (default false) — when true, drops the frequency
%       recovery energy from the receiver total, modelling the limit
%       where FR cost is amortised over infinite symbols.

if nargin < 1 || isempty(ignore_fr_energy)
    ignore_fr_energy = false;
end

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src'));

%% Load grid-sweep results -----------------------------------------------
d = load(fullfile(fileparts(mfilename('fullpath')), ...
                  'bit_width_full_grid_sweep.mat'));
sweep = d.tbl;

if ignore_fr_energy
    E_rx_fJ = sweep.energy_cr_fJ;                      % FR cost amortised away
else
    E_rx_fJ = sweep.energy_fr_fJ + sweep.energy_cr_fJ; % total Rx DSP energy [fJ/bit]
end

%% Per-row trial statistics ---------------------------------------------
[mean_snr_row, std_snr_row] = row_snr_stats(sweep.fec_snr_trials);

%% Per-algorithm-combo summary (best-FEC-SNR precision combo) ----------
tbl = build_summary_table(sweep, E_rx_fJ, mean_snr_row, std_snr_row);

fprintf('\n=== Lowest-FEC-SNR precision combo per algorithm pair ===\n');
disp(tbl);

%% Plots -----------------------------------------------------------------
plot_cr_precision_comparison(sweep, "differential_kay", ...
    ["pilots_only", "viterbi_viterbi"], ...
    ["Pilots only (PO)", "Viterbi-Viterbi (V&V)"]);

plot_cr_energy_saving(sweep, "differential_kay", ...
    "viterbi_viterbi", "pilots_only");

end


function [m, s] = row_snr_stats(trials)
    N = numel(trials);
    m = nan(N, 1);
    s = nan(N, 1);
    for r = 1:N
        t = trials{r};
        t = t(~isnan(t));
        if ~isempty(t)
            m(r) = mean(t);
            s(r) = std(t);
        end
    end
end


function out = build_summary_table(sweep, E_rx_fJ, mean_snr_row, std_snr_row)

    combos = unique(sweep(:, {'fr_algo', 'cr_algo'}), 'rows');
    NC = height(combos);

    v_fr  = strings(NC, 1);  v_cr   = strings(NC, 1);
    v_flFR = nan(NC, 1);     v_flCR = nan(NC, 1);   v_bd  = nan(NC, 1);
    v_Efr = nan(NC, 1);      v_Ecr  = nan(NC, 1);   v_E   = nan(NC, 1);
    v_mSN = nan(NC, 1);      v_sSN  = nan(NC, 1);

    for ci = 1:NC
        fr = combos.fr_algo(ci);
        cr = combos.cr_algo(ci);
        idx = find(sweep.fr_algo == fr & sweep.cr_algo == cr);

        % Precision combo with minimum mean FEC SNR (best DSP performance)
        snr_sub = mean_snr_row(idx);
        [best, ki] = min(snr_sub);
        if ~isfinite(best)
            continue
        end
        r = idx(ki);

        v_fr(ci)   = fr;
        v_cr(ci)   = cr;
        v_flFR(ci) = sweep.fl_fr(r);
        v_flCR(ci) = sweep.fl_cr(r);
        v_bd(ci)   = sweep.blind_d(r);
        v_Efr(ci)  = sweep.energy_fr_fJ(r);
        v_Ecr(ci)  = sweep.energy_cr_fJ(r);
        v_E(ci)    = E_rx_fJ(r);
        v_mSN(ci)  = mean_snr_row(r);
        v_sSN(ci)  = std_snr_row(r);
    end

    out = table(v_fr, v_cr, v_flFR, v_flCR, v_bd, ...
                v_Efr, v_Ecr, v_E, v_mSN, v_sSN, ...
        'VariableNames', {'fr_algo', 'cr_algo', 'fl_fr', 'fl_cr', 'blind_d', ...
        'fr_energy_fJ', 'cr_energy_fJ', 'rx_energy_fJ', ...
        'mean_fec_snr_dB', 'std_fec_snr_dB'});
end


function plot_cr_precision_comparison(sweep, fr_target, cr_list, cr_labels)
% Compare FEC SNR vs CR fractional bit width across CR algorithms while
% holding the FR algorithm fixed.  FR precision is pinned to its highest
% value in the sweep so the curve isolates the CR-precision effect.

    fl_fr_max = max(sweep.fl_fr(sweep.fr_algo == fr_target));

    figure('Name', sprintf('CR precision sweep (FR = %s)', ...
        abbrevAlgo(fr_target)), 'Color', 'w', ...
        'Position', [260, 260, 900, 500]);
    ax = axes; hold(ax, 'on');

    N = numel(cr_list);
    cmap = lines(N);
    h    = gobjects(N, 1);

    for ai = 1:N
        cr = cr_list(ai);
        mask = sweep.fr_algo == fr_target & sweep.cr_algo == cr & ...
               sweep.fl_fr == fl_fr_max;
        rows = find(mask);
        if isempty(rows)
            continue
        end

        fl_cr_vals = unique(sweep.fl_cr(rows));
        m = nan(numel(fl_cr_vals), 1);
        s = nan(numel(fl_cr_vals), 1);
        for fi = 1:numel(fl_cr_vals)
            sub = rows(sweep.fl_cr(rows) == fl_cr_vals(fi));
            trials = vertcat(sweep.fec_snr_trials{sub});
            trials = trials(~isnan(trials));
            if ~isempty(trials)
                m(fi) = mean(trials);
                s(fi) = std(trials);
            end
        end

        good = ~isnan(m);
        h(ai) = errorbar(ax, fl_cr_vals(good), m(good), s(good), '-o', ...
            'Color', cmap(ai, :), 'LineWidth', 1.4, 'MarkerSize', 5, ...
            'MarkerFaceColor', cmap(ai, :), 'CapSize', 6);
    end

    grid(ax, 'on');
    xlabel(ax, 'CR fractional bit width fl_{CR}', 'FontSize', 11);
    ylabel(ax, 'FEC SNR [dB]', 'FontSize', 11);
    title(ax, sprintf('FEC SNR vs CR precision (FR = %s, fl_{FR} = %d)', ...
        abbrevAlgo(fr_target), fl_fr_max), 'FontSize', 11);
    legend(ax, h(isgraphics(h)), cr_labels(isgraphics(h)), ...
        'Location', 'best', 'FontSize', 14, 'Interpreter', 'none');
end


function plot_cr_energy_saving(sweep, fr_target, cr_baseline, cr_alt)
% Plots CR DSP energy saved by switching from CR_BASELINE to CR_ALT as a
% function of CR fractional bit width, with FR fixed to its
% highest-precision setting (same selection as the FEC SNR plot).
%   saving(fl_cr) = energy_cr(baseline) - energy_cr(alt)

    fl_fr_max = max(sweep.fl_fr(sweep.fr_algo == fr_target));

    e_base = cr_energy_curve(sweep, fr_target, cr_baseline, fl_fr_max);
    e_alt  = cr_energy_curve(sweep, fr_target, cr_alt,      fl_fr_max);

    [fl_cr_vals, saving] = align_and_diff(e_base, e_alt);

    figure('Name', sprintf('CR energy saving (FR = %s)', ...
        abbrevAlgo(fr_target)), 'Color', 'w', ...
        'Position', [280, 280, 900, 500]);
    ax = axes;

    plot(ax, fl_cr_vals, saving, '-o', ...
        'LineWidth', 1.6, 'MarkerSize', 6, 'MarkerFaceColor', 'auto');

    grid(ax, 'on');
    xlabel(ax, 'CR fractional bit width fl_{CR}', 'FontSize', 11);
    ylabel(ax, 'CR energy saved [fJ/bit]', 'FontSize', 11);
    title(ax, sprintf('CR energy saved switching %s \\rightarrow %s (FR = %s, fl_{FR} = %d)', ...
        abbrevAlgo(cr_baseline), abbrevAlgo(cr_alt), ...
        abbrevAlgo(fr_target), fl_fr_max), 'FontSize', 11);
    legend(ax, sprintf('%s \\rightarrow %s', ...
        abbrevAlgo(cr_baseline), abbrevAlgo(cr_alt)), ...
        'Location', 'best', 'FontSize', 14);
end


function curve = cr_energy_curve(sweep, fr_target, cr, fl_fr_max)
    mask = sweep.fr_algo == fr_target & sweep.cr_algo == cr & ...
           sweep.fl_fr == fl_fr_max;
    rows = find(mask);
    fl_cr_vals = unique(sweep.fl_cr(rows));
    e = nan(numel(fl_cr_vals), 1);
    for fi = 1:numel(fl_cr_vals)
        sub = rows(sweep.fl_cr(rows) == fl_cr_vals(fi));
        e(fi) = mean(sweep.energy_cr_fJ(sub));
    end
    curve = struct('fl_cr', fl_cr_vals, 'e', e);
end


function [fl_cr_vals, diff] = align_and_diff(a, b)
    fl_cr_vals = intersect(a.fl_cr, b.fl_cr);
    diff = nan(numel(fl_cr_vals), 1);
    for i = 1:numel(fl_cr_vals)
        ia = find(a.fl_cr == fl_cr_vals(i), 1);
        ib = find(b.fl_cr == fl_cr_vals(i), 1);
        diff(i) = a.e(ia) - b.e(ib);
    end
end


function s = abbrevAlgo(name)
    map = {'fft_search',       'R&B';       ...
           'fft_search_blind', 'R&B blind'; ...
           'differential_kay', 'DK';        ...
           'viterbi_viterbi',  'V&V';       ...
           'pilots_only',      'PO'};
    idx = strcmp(map(:, 1), name);
    if any(idx)
        s = map{idx, 2};
    else
        s = name;
    end
end
