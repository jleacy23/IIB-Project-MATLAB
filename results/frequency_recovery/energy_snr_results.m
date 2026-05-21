function tbl = energy_snr_results(ignore_fr_energy)
% ENERGY_SNR_RESULTS  Receiver DSP energy and FEC SNR analysis.
%
%   For each FR+CR algorithm pair, returns a summary table reporting the
%   precision combo with the lowest mean FEC SNR (best DSP performance)
%   and its associated Rx energy.  Produces:
%     - CR-precision comparison (pilots-only vs V&V, FR = DK)
%     - CR energy saving switching V&V to pilots-only
%     - FR-precision comparison (DK data-aided, R&B data-aided,
%       R&B blind D=512), CR = pilots-only
%     - FR energy saving switching to DK, amortised over 3712 symbols
%
%   ignore_fr_energy (default false) — when true, drops the frequency
%       recovery energy from the receiver total in the summary table,
%       modelling the limit where FR cost is amortised over infinite
%       symbols.

if nargin < 1 || isempty(ignore_fr_energy)
    ignore_fr_energy = false;
end

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src'));

%% Load grid-sweep results -----------------------------------------------
d = load(fullfile(fileparts(mfilename('fullpath')), ...
                  'bit_width_full_grid_sweep.mat'));
sweep = d.tbl;

%% Energy model — edit here to explore different estimates without
%% rerunning the grid sweep.  Mirrors bit_width_full constants.
EP = struct( ...
    'EAdd_fJ',      3.16/4, ...
    'EMult_fJ',     3.03/4, ...
    'M',            4,    ...
    'Oversampling', 1,    ...
    'TrainingLen',  11,   ...
    'FR_Nfft',      512,  ...
    'SubframeLen',  3712, ...
    'BlockLen',     32);

[sweep.energy_fr_fJ, sweep.energy_cr_fJ] = recompute_energies(sweep, EP);

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
% CR comparison (FR fixed to differential Kay) ---------------------------
plot_cr_precision_comparison(sweep, "differential_kay", ...
    ["pilots_only", "viterbi_viterbi"], ...
    ["Pilots only (PO)", "Viterbi-Viterbi (V&V)"]);

plot_cr_energy_saving(sweep, "differential_kay", ...
    "viterbi_viterbi", "pilots_only");

% FR comparison (CR fixed to pilots-only) --------------------------------
plot_fr_precision_comparison(sweep, "pilots_only", ...
    ["differential_kay", "fft_search", "fft_search_blind"], ...
    [NaN, NaN, 512], ...
    ["DK (data-aided)", "R&B (data-aided)", "R&B (blind, D=512)"]);

plot_fr_energy_saving(sweep, "pilots_only", "differential_kay", ...
    ["fft_search", "fft_search_blind"], ...
    [NaN, 512], ...
    ["R&B (data-aided)", "R&B (blind, D=512)"]);

% Per-combo energy breakdown + summary table at low-precision points -----
plot_and_print_combo_breakdown(sweep, EP, mean_snr_row, std_snr_row);

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

    [x_vals, saving] = align_and_diff(e_base, e_alt);

    figure('Name', sprintf('CR energy saving (FR = %s)', ...
        abbrevAlgo(fr_target)), 'Color', 'w', ...
        'Position', [280, 280, 900, 500]);
    ax = axes;

    plot(ax, x_vals, saving, '-o', ...
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
    curve = struct('x', fl_cr_vals, 'y', e);
end


function plot_fr_precision_comparison(sweep, cr_target, fr_algos, fr_blind_d, fr_labels)
% Compare FEC SNR vs FR fractional bit width across FR algorithm variants
% (each variant specified by an algorithm name and optional blind_d filter,
% NaN to disable).  CR precision pinned to its highest value with CR
% algorithm fixed to CR_TARGET so the curve isolates the FR-precision
% effect.

    fl_cr_max = max(sweep.fl_cr(sweep.cr_algo == cr_target));

    figure('Name', sprintf('FR precision sweep (CR = %s)', ...
        abbrevAlgo(cr_target)), 'Color', 'w', ...
        'Position', [300, 300, 900, 500]);
    ax = axes; hold(ax, 'on');

    N = numel(fr_algos);
    cmap = lines(N);
    h    = gobjects(N, 1);

    for ai = 1:N
        rows = fr_rows(sweep, fr_algos(ai), fr_blind_d(ai), cr_target, fl_cr_max);
        if isempty(rows)
            continue
        end

        fl_fr_vals = unique(sweep.fl_fr(rows));
        m = nan(numel(fl_fr_vals), 1);
        s = nan(numel(fl_fr_vals), 1);
        for fi = 1:numel(fl_fr_vals)
            sub = rows(sweep.fl_fr(rows) == fl_fr_vals(fi));
            trials = vertcat(sweep.fec_snr_trials{sub});
            trials = trials(~isnan(trials));
            if ~isempty(trials)
                m(fi) = mean(trials);
                s(fi) = std(trials);
            end
        end

        good = ~isnan(m);
        h(ai) = errorbar(ax, fl_fr_vals(good), m(good), s(good), '-o', ...
            'Color', cmap(ai, :), 'LineWidth', 1.4, 'MarkerSize', 5, ...
            'MarkerFaceColor', cmap(ai, :), 'CapSize', 6);
    end

    grid(ax, 'on');
    xlabel(ax, 'FR fractional bit width fl_{FR}', 'FontSize', 11);
    ylabel(ax, 'FEC SNR [dB]', 'FontSize', 11);
    title(ax, sprintf('FEC SNR vs FR precision (CR = %s, fl_{CR} = %d)', ...
        abbrevAlgo(cr_target), fl_cr_max), 'FontSize', 11);
    legend(ax, h(isgraphics(h)), fr_labels(isgraphics(h)), ...
        'Location', 'best', 'FontSize', 14, 'Interpreter', 'none');
end


function plot_fr_energy_saving(sweep, cr_target, fr_target, fr_alts, ...
                               fr_alt_blind_d, fr_alt_labels)
% Plots FR DSP energy saved by switching from each FR_ALTS variant to
% FR_TARGET, as a function of FR fractional bit width.  FR energy from
% the sweep is already amortised over the CPON subframe (SubframeLen in
% recompute_energies).  CR fixed to CR_TARGET at its highest precision
% (matches the FEC SNR plot's selection).
%   saving(fl_fr) = E_fr(alt) - E_fr(target)

    fl_cr_max = max(sweep.fl_cr(sweep.cr_algo == cr_target));

    e_tgt = fr_energy_curve(sweep, fr_target, NaN, cr_target, fl_cr_max);

    figure('Name', sprintf('FR energy saving (switch to %s)', ...
        abbrevAlgo(fr_target)), 'Color', 'w', ...
        'Position', [320, 320, 900, 500]);
    ax = axes; hold(ax, 'on');

    N = numel(fr_alts);
    cmap = lines(N);
    h    = gobjects(N, 1);
    leg  = strings(N, 1);

    for ai = 1:N
        e_alt = fr_energy_curve(sweep, fr_alts(ai), fr_alt_blind_d(ai), ...
                                cr_target, fl_cr_max);
        [x_vals, saving] = align_and_diff(e_alt, e_tgt);

        h(ai) = plot(ax, x_vals, saving, '-o', ...
            'Color', cmap(ai, :), 'LineWidth', 1.6, ...
            'MarkerSize', 6, 'MarkerFaceColor', cmap(ai, :));
        leg(ai) = sprintf('%s \\rightarrow %s', ...
            fr_alt_labels(ai), abbrevAlgo(fr_target));
    end

    grid(ax, 'on');
    xlabel(ax, 'FR fractional bit width fl_{FR}', 'FontSize', 11);
    ylabel(ax, 'FR energy saved [fJ/bit]', 'FontSize', 11);
    title(ax, sprintf('FR energy saved switching to %s', ...
        abbrevAlgo(fr_target)), 'FontSize', 11);
    legend(ax, h(isgraphics(h)), leg(isgraphics(h)), ...
        'Location', 'best', 'FontSize', 14);
end


function rows = fr_rows(sweep, fr_algo, blind_d, cr_target, fl_cr_max)
    mask = sweep.fr_algo == fr_algo & sweep.cr_algo == cr_target & ...
           sweep.fl_cr == fl_cr_max;
    if ~isnan(blind_d)
        mask = mask & sweep.blind_d == blind_d;
    end
    rows = find(mask);
end


function curve = fr_energy_curve(sweep, fr_algo, blind_d, cr_target, fl_cr_max)
    rows = fr_rows(sweep, fr_algo, blind_d, cr_target, fl_cr_max);
    fl_fr_vals = unique(sweep.fl_fr(rows));
    e = nan(numel(fl_fr_vals), 1);
    for fi = 1:numel(fl_fr_vals)
        sub = rows(sweep.fl_fr(rows) == fl_fr_vals(fi));
        e(fi) = mean(sweep.energy_fr_fJ(sub));
    end
    curve = struct('x', fl_fr_vals, 'y', e);
end


function [x_vals, diff] = align_and_diff(a, b)
    x_vals = intersect(a.x, b.x);
    diff = nan(numel(x_vals), 1);
    for i = 1:numel(x_vals)
        ia = find(a.x == x_vals(i), 1);
        ib = find(b.x == x_vals(i), 1);
        diff(i) = a.y(ia) - b.y(ib);
    end
end


function plot_and_print_combo_breakdown(sweep, P, mean_snr_row, std_snr_row)
% Stacked bar chart + summary table for two FR+CR combos at the (2,2) and
% (4,4) (fl_fr, fl_cr) operating points.  Stack components are:
%   FR energy, CR energy, and the single real multiplication per symbol
%   that applies the composed CFO + phase rotation.

    fr_algos = ["differential_kay",                  "fft_search_blind"];
    cr_algos = ["pilots_only",                       "pilots_only"];
    blind_ds = [NaN,                                  512];
    labels   = ["Diff. phase + Kay + Pilots only",   "R&B (blind, D=512) + Pilots only"];
    fl_pairs = [2 2; 4 4];

    NC = numel(fr_algos);
    NP = size(fl_pairs, 1);
    N  = NC * NP;

    E          = zeros(N, 3);   % cols = [FR, CR, apply]
    fl_fr_col  = nan(N, 1);
    fl_cr_col  = nan(N, 1);
    combo_col  = strings(N, 1);
    x_pos      = nan(N, 1);
    mean_snr   = nan(N, 1);
    std_snr    = nan(N, 1);

    GROUP_GAP = 1;  % extra x-axis units between different FR+CR combos

    bi = 0;
    for ci = 1:NC
        for pi = 1:NP
            bi    = bi + 1;
            fl_fr = fl_pairs(pi, 1);
            fl_cr = fl_pairs(pi, 2);
            r = find_combo_row(sweep, fr_algos(ci), cr_algos(ci), ...
                               fl_fr, fl_cr, blind_ds(ci));

            E(bi, 1) = sweep.energy_fr_fJ(r);
            E(bi, 2) = sweep.energy_cr_fJ(r);
            E(bi, 3) = apply_mult_energy(fl_cr, P);

            fl_fr_col(bi)  = fl_fr;
            fl_cr_col(bi)  = fl_cr;
            combo_col(bi)  = labels(ci);
            x_pos(bi)      = (ci - 1) * (NP + GROUP_GAP) + pi;
            mean_snr(bi)   = mean_snr_row(r);
            std_snr(bi)    = std_snr_row(r);
        end
    end

    % ---- Bar chart -------------------------------------------------------
    figure('Name', 'Rx DSP energy breakdown', 'Color', 'w', ...
        'Position', [340, 340, 900, 500]);
    ax = axes;
    bar(ax, x_pos, E, 'stacked');

    group_centers = ((1:NC) - 1) * (NP + GROUP_GAP) + (NP + 1) / 2;
    set(ax, 'XTick', group_centers, 'XTickLabel', labels);
    xlim(ax, [min(x_pos) - 1, max(x_pos) + 1]);
    ylabel(ax, 'Energy [fJ/bit]', 'FontSize', 11);
    legend(ax, {'FR', 'CR (phase)', 'CFO+phase apply (1 mult/sym)'}, ...
        'Location', 'best', 'FontSize', 12);
    title(ax, 'Receiver DSP energy breakdown', 'FontSize', 11);
    grid(ax, 'on');

    % Per-bar (fl_fr, fl_cr) annotation above each stack
    totals = sum(E, 2);
    yl     = ylim(ax);
    for bi = 1:N
        text(ax, x_pos(bi), totals(bi) + 0.015 * yl(2), ...
            sprintf('(%d,%d)', fl_fr_col(bi), fl_cr_col(bi)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 10);
    end

    % ---- Table -----------------------------------------------------------
    T = table(combo_col, fl_fr_col, fl_cr_col, ...
              E(:, 1), E(:, 2), E(:, 3), sum(E, 2), ...
              mean_snr, std_snr, ...
        'VariableNames', {'combo', 'fl_fr', 'fl_cr', ...
        'E_fr_fJ', 'E_cr_fJ', 'E_apply_fJ', 'E_total_fJ', ...
        'mean_fec_snr_dB', 'std_fec_snr_dB'});

    fprintf('\n=== Energy breakdown @ low-precision combos ===\n');
    disp(T);
end


function r = find_combo_row(sweep, fr, cr, fl_fr, fl_cr, blind_d)
    mask = sweep.fr_algo == fr & sweep.cr_algo == cr & ...
           sweep.fl_fr  == fl_fr & sweep.fl_cr == fl_cr;
    if isnan(blind_d)
        mask = mask & isnan(sweep.blind_d);
    else
        mask = mask & sweep.blind_d == blind_d;
    end
    r = find(mask, 1);
end


function E = apply_mult_energy(fl, P)
    % Energy/bit of one real multiplication per symbol — applies the
    % composed CFO + phase rotation to every output sample.
    E = energy.receiver(0, 1, P.EAdd_fJ, P.EMult_fJ, ...
                        P.M, P.Oversampling, fl);
end


function [E_fr, E_cr] = recompute_energies(sweep, P)
    N    = height(sweep);
    E_fr = nan(N, 1);
    E_cr = nan(N, 1);
    for r = 1:N
        E_fr(r) = fr_row_energy(sweep.fr_algo(r), sweep.fl_fr(r), ...
                                sweep.blind_d(r), P);
        E_cr(r) = cr_row_energy(sweep.cr_algo(r), sweep.fl_cr(r), P);
    end
end


function E = fr_row_energy(algo, fl, blind_d, P)
    switch algo
        case "fft_search"
            [NM, NA] = fft_search_counts(P.TrainingLen, P.FR_Nfft, false);
        case "differential_kay"
            [NM, NA] = differential_kay_counts(P.TrainingLen);
        case "fft_search_blind"
            [NM, NA] = fft_search_counts(blind_d, P.FR_Nfft, true);
        otherwise
            error('energy_snr_results:unknownFR', ...
                  'Unknown FR algo: %s', algo);
    end
    NM = NM / P.SubframeLen;
    NA = NA / P.SubframeLen;
    E  = energy.receiver(NA, NM, P.EAdd_fJ, P.EMult_fJ, ...
                         P.M, P.Oversampling, fl);
end


function E = cr_row_energy(algo, fl, P)
    switch algo
        case "viterbi_viterbi"
            [NM, NA] = viterbi_counts(P.BlockLen);
        case "pilots_only"
            [NM, NA] = pilots_only_counts();
        otherwise
            error('energy_snr_results:unknownCR', ...
                  'Unknown CR algo: %s', algo);
    end
    NM = NM / P.BlockLen;
    NA = NA / P.BlockLen;
    E  = energy.receiver(NA, NM, P.EAdd_fJ, P.EMult_fJ, ...
                         P.M, P.Oversampling, fl);
end


function [NM, NA] = fft_search_counts(L, Nfft, blind)
    if blind
        NM_form = 12 * L;
        NA_form = 9  * L;
    else
        NM_form = 4 * L;
        NA_form = 3 * L;
    end
    NM_fft    = 2 * L * log2(Nfft / L) + 2 * Nfft * log2(L);
    NA_fft    = 3 * L * log2(Nfft / L) + 3 * Nfft * log2(L);
    NM_search = 2 * Nfft;
    NA_search = 2 * Nfft;
    NM_interp = 5;
    NA_interp = 4;
    NM = NM_form + NM_fft + NM_search + NM_interp;
    NA = NA_form + NA_fft + NA_search + NA_interp;
end


function [NM, NA] = differential_kay_counts(L)
    NM = 7 * L + 1;
    NA = 4 * L - 2;
end


function [NM, NA] = viterbi_counts(N)
    NM = 16 * N + 1;
    NA = 14 * N - 1;
end


function [NM, NA] = pilots_only_counts()
    NM = 1;
    NA = 1;
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
