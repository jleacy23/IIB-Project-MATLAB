classdef plot_energy_snr < matlab.unittest.TestCase
%PLOT_ENERGY_SNR  Visualise bit_width_full grid-sweep results.
%
%   Loads the long-form table produced by bit_width_full
%       (bit_width_full_grid_sweep.mat)
%   and plots minimum-energy operating points against FEC SNR for each
%   (FR algo, CR algo) combination.
%
%   For every (fr_algo, cr_algo) pair the plot shows the lower envelope
%   (Pareto front) over all (fl_fr, fl_cr) bit-width combinations
%   sampled — i.e. for each FEC SNR value, the cheapest bit-width combo
%   that achieves it.  For the blind FFT search variant the
%   optimisation also runs over BlindD.
%
%   Run with:
%       runtests('plot_energy_snr')

    %% ================================================================
    %  Parameters
    %% ================================================================
    properties (Constant)
        FrAlgos = ["fft_search", "differential_kay", "fft_search_blind"]
        CrAlgos = ["viterbi_viterbi", "pilots_only"]
        FrLabel = struct( ...
            'fft_search',       'FFT Search (DA)', ...
            'differential_kay', 'Diff. Kay (DA)', ...
            'fft_search_blind', 'FFT Search (blind)')
        CrLabel = struct( ...
            'viterbi_viterbi', 'V&V', ...
            'pilots_only',     'Pilots-only')
        GridSweepFile = 'bit_width_full_grid_sweep.mat'
    end

    %% ================================================================
    %  Test class setup
    %% ================================================================
    methods (TestClassSetup)
        function setupPath(~)
            % +energy/* lives under src/, needed for the system-energy table.
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src'));
        end
    end

    %% ================================================================
    %  Tests
    %% ================================================================
    methods (Test)

        function test_plot_energy_vs_snr(testCase)  %#ok<MANU>
            dataDir = fileparts(mfilename('fullpath'));
            T = plot_energy_snr.loadCombinedTable(dataDir);

            FrAlgos_l = plot_energy_snr.FrAlgos;
            CrAlgos_l = plot_energy_snr.CrAlgos;

            nCombos = numel(FrAlgos_l) * numel(CrAlgos_l);
            colors  = lines(nCombos);
            markers = {'o', 's', '^', 'd', 'v', 'p'};

            figure('Name', 'Energy vs FEC SNR — min-energy Pareto front', ...
                'Color', 'w', 'Position', [100 100 900 600]);
            ax = axes;
            hold(ax, 'on');
            grid(ax, 'on');

            SNR_MAX = 10;   % dB — clip the x-axis to ignore high-SNR outliers
            globalMinSnr = Inf;

            k = 0;
            for fi = 1:numel(FrAlgos_l)
                for ci = 1:numel(CrAlgos_l)
                    k  = k + 1;
                    fr = FrAlgos_l(fi);
                    cr = CrAlgos_l(ci);

                    sel = T.fr_algo == fr & T.cr_algo == cr & ~isnan(T.fec_snr_db);
                    sub = T(sel, :);
                    if isempty(sub), continue; end

                    [snr_pf, e_pf] = plot_energy_snr.paretoFront( ...
                        sub.fec_snr_db, sub.energy_total_fJ);

                    inRange = snr_pf <= SNR_MAX;
                    snr_pf  = snr_pf(inRange);
                    e_pf    = e_pf(inRange);
                    if isempty(snr_pf), continue; end
                    globalMinSnr = min(globalMinSnr, min(snr_pf));

                    label = sprintf('%s + %s', ...
                        plot_energy_snr.FrLabel.(char(fr)), ...
                        plot_energy_snr.CrLabel.(char(cr)));

                    plot(ax, snr_pf, e_pf, ...
                        'LineStyle', '-', 'Marker', markers{k}, ...
                        'MarkerSize', 6, 'LineWidth', 1.8, ...
                        'Color', colors(k, :), ...
                        'DisplayName', label);
                end
            end

            if isfinite(globalMinSnr)
                xlim(ax, [globalMinSnr, SNR_MAX]);
            end
            set(ax, 'YScale', 'log', 'FontSize', 11, 'Box', 'on');
            xlabel(ax, 'FEC SNR threshold  [dB]', 'FontSize', 12);
            ylabel(ax, 'Total energy per bit  [fJ]', 'FontSize', 12);
            title(ax,  'Min-energy bit-width combination per FR \times CR pair', ...
                'FontSize', 12);
            legend(ax, 'Location', 'best', 'FontSize', 10, 'Interpreter', 'none');
        end

        function test_plot_energy_breakdown(testCase)  %#ok<MANU>
            % For each (FR, CR) pair and a few FEC SNR targets, plot the
            % min-energy point as a stacked bar split into FR vs CR energy.
            dataDir = fileparts(mfilename('fullpath'));
            T = plot_energy_snr.loadCombinedTable(dataDir);

            FrAlgos_l  = plot_energy_snr.FrAlgos;
            CrAlgos_l  = plot_energy_snr.CrAlgos;
            snrTargets = [10, 15, 20];

            nCombos = numel(FrAlgos_l) * numel(CrAlgos_l);
            labels  = strings(nCombos, 1);
            E_fr    = nan(nCombos, numel(snrTargets));
            E_cr    = nan(nCombos, numel(snrTargets));

            k = 0;
            for fi = 1:numel(FrAlgos_l)
                for ci = 1:numel(CrAlgos_l)
                    k  = k + 1;
                    fr = FrAlgos_l(fi);
                    cr = CrAlgos_l(ci);
                    sub = T(T.fr_algo == fr & T.cr_algo == cr & ~isnan(T.fec_snr_db), :);
                    labels(k) = sprintf('%s + %s', ...
                        plot_energy_snr.FrLabel.(char(fr)), ...
                        plot_energy_snr.CrLabel.(char(cr)));

                    for ti = 1:numel(snrTargets)
                        meet = sub(sub.fec_snr_db <= snrTargets(ti), :);
                        if isempty(meet), continue; end
                        [~, idx] = min(meet.energy_total_fJ);
                        E_fr(k, ti) = meet.energy_fr_fJ(idx);
                        E_cr(k, ti) = meet.energy_cr_fJ(idx);
                    end
                end
            end

            figure('Name', 'Energy breakdown (FR vs CR) at SNR targets', ...
                'Color', 'w', 'Position', [120 120 1200 500]);
            for ti = 1:numel(snrTargets)
                subplot(1, numel(snrTargets), ti);
                bar([E_fr(:, ti), E_cr(:, ti)], 'stacked');
                set(gca, 'XTickLabel', labels, 'XTickLabelRotation', 30, ...
                    'FontSize', 9, 'Box', 'on');
                ylabel('Energy per bit  [fJ]');
                title(sprintf('Min-energy config meeting FEC SNR \\leq %d dB', snrTargets(ti)));
                legend({'FR', 'CR'}, 'Location', 'northwest');
                grid on;
            end
        end

        function test_plot_blind_d_sweep(testCase)  %#ok<MANU>
            % Pareto fronts for the blind FFT search variant, one curve per
            % BlindD value, broken out by CR algorithm.
            dataDir = fileparts(mfilename('fullpath'));
            T = plot_energy_snr.loadCombinedTable(dataDir);

            CrAlgos_l = plot_energy_snr.CrAlgos;
            figure('Name', 'Blind FFT — accuracy vs energy across BlindD', ...
                'Color', 'w', 'Position', [140 140 1100 500]);

            for ci = 1:numel(CrAlgos_l)
                subplot(1, numel(CrAlgos_l), ci);
                cr  = CrAlgos_l(ci);
                sub = T(T.fr_algo == "fft_search_blind" & T.cr_algo == cr & ...
                        ~isnan(T.fec_snr_db), :);

                bds  = unique(sub.blind_d(~isnan(sub.blind_d)));
                cmap = parula(max(numel(bds), 2));
                hold on; grid on;
                for bi = 1:numel(bds)
                    sb = sub(sub.blind_d == bds(bi), :);
                    [snr_pf, e_pf] = plot_energy_snr.paretoFront( ...
                        sb.fec_snr_db, sb.energy_total_fJ);
                    plot(snr_pf, e_pf, 'o-', 'Color', cmap(bi, :), ...
                        'LineWidth', 1.6, 'MarkerSize', 5, ...
                        'DisplayName', sprintf('D = %d', bds(bi)));
                end
                set(gca, 'YScale', 'log', 'FontSize', 11, 'Box', 'on');
                xlabel('FEC SNR threshold  [dB]');
                ylabel('Total energy per bit  [fJ]');
                title(sprintf('Blind FFT + %s', ...
                    plot_energy_snr.CrLabel.(char(cr))));
                legend('Location', 'best');
            end
        end

        function test_plot_snr_heatmaps(testCase)  %#ok<MANU>
            % For each (FR, CR) pair, heatmap of FEC SNR threshold over the
            % (fl_fr, fl_cr) grid.  For the blind FFT variant we collapse
            % BlindD by taking the best (lowest) FEC SNR per cell.
            dataDir = fileparts(mfilename('fullpath'));
            T = plot_energy_snr.loadCombinedTable(dataDir);

            FrAlgos_l = plot_energy_snr.FrAlgos;
            CrAlgos_l = plot_energy_snr.CrAlgos;
            fl_vec    = unique(T.fl_fr);

            figure('Name', 'FEC SNR threshold heatmaps', 'Color', 'w', ...
                'Position', [160 100 1200 720]);

            nFR = numel(FrAlgos_l);
            nCR = numel(CrAlgos_l);
            k = 0;
            for fi = 1:nFR
                for ci = 1:nCR
                    k  = k + 1;
                    subplot(nFR, nCR, k);
                    fr  = FrAlgos_l(fi);
                    cr  = CrAlgos_l(ci);
                    sub = T(T.fr_algo == fr & T.cr_algo == cr, :);

                    G = nan(numel(fl_vec));
                    for r = 1:numel(fl_vec)
                        for c = 1:numel(fl_vec)
                            sel  = sub.fl_fr == fl_vec(r) & sub.fl_cr == fl_vec(c);
                            vals = sub.fec_snr_db(sel);
                            vals = vals(~isnan(vals));
                            if ~isempty(vals)
                                G(r, c) = min(vals);
                            end
                        end
                    end
                    imagesc(fl_vec, fl_vec, G, 'AlphaData', ~isnan(G));
                    set(gca, 'YDir', 'normal', 'Color', [0.9 0.9 0.9]);
                    colorbar;
                    xlabel('FL_{CR}'); ylabel('FL_{FR}');
                    title(sprintf('%s + %s', ...
                        plot_energy_snr.FrLabel.(char(fr)), ...
                        plot_energy_snr.CrLabel.(char(cr))), ...
                        'FontSize', 10);
                end
            end
            sgtitle('FEC SNR threshold [dB] (lower = better, NaN = no convergence)', ...
                'FontSize', 13);
        end

        function test_plot_bitwidth_along_pareto(testCase)  %#ok<MANU>
            % For each (FR, CR) pair, show which fl_fr / fl_cr values land
            % on the Pareto front as the FEC SNR target relaxes.
            dataDir = fileparts(mfilename('fullpath'));
            T = plot_energy_snr.loadCombinedTable(dataDir);

            FrAlgos_l = plot_energy_snr.FrAlgos;
            CrAlgos_l = plot_energy_snr.CrAlgos;

            figure('Name', 'Pareto-optimal bit widths vs FEC SNR target', ...
                'Color', 'w', 'Position', [180 120 1200 720]);

            nFR = numel(FrAlgos_l);
            nCR = numel(CrAlgos_l);
            k = 0;
            for fi = 1:nFR
                for ci = 1:nCR
                    k  = k + 1;
                    subplot(nFR, nCR, k);
                    fr  = FrAlgos_l(fi);
                    cr  = CrAlgos_l(ci);
                    sub = T(T.fr_algo == fr & T.cr_algo == cr & ~isnan(T.fec_snr_db), :);
                    if isempty(sub), continue; end

                    keep  = plot_energy_snr.paretoMask(sub.fec_snr_db, sub.energy_total_fJ);
                    sub_p = sub(keep, :);
                    [~, ord] = sort(sub_p.fec_snr_db);
                    sub_p = sub_p(ord, :);

                    plot(sub_p.fec_snr_db, sub_p.fl_fr, 'o-', ...
                        'LineWidth', 1.6, 'MarkerSize', 6, ...
                        'DisplayName', 'FL_{FR}');
                    hold on; grid on;
                    plot(sub_p.fec_snr_db, sub_p.fl_cr, 's--', ...
                        'LineWidth', 1.6, 'MarkerSize', 6, ...
                        'DisplayName', 'FL_{CR}');
                    ylim([0 max(plot_energy_snr.FL_vec) + 2]);
                    xlabel('FEC SNR threshold  [dB]');
                    ylabel('Fractional bits');
                    title(sprintf('%s + %s', ...
                        plot_energy_snr.FrLabel.(char(fr)), ...
                        plot_energy_snr.CrLabel.(char(cr))), ...
                        'FontSize', 10);
                    legend('Location', 'best');
                end
            end
            sgtitle('Bit-width selection along the energy Pareto front', ...
                'FontSize', 13);
        end

        function test_print_summary_tables(testCase)  %#ok<MANU>
            % Print (and save as CSV) two summary tables per algorithm pair:
            %   * configuration achieving the lowest FEC SNR (best sensitivity)
            %   * configuration achieving the lowest total energy
            dataDir = fileparts(mfilename('fullpath'));
            T = plot_energy_snr.loadCombinedTable(dataDir);

            FrAlgos_l = plot_energy_snr.FrAlgos;
            CrAlgos_l = plot_energy_snr.CrAlgos;

            rowsSnr = {};
            rowsErg = {};
            for fi = 1:numel(FrAlgos_l)
                for ci = 1:numel(CrAlgos_l)
                    fr  = FrAlgos_l(fi);
                    cr  = CrAlgos_l(ci);
                    sub = T(T.fr_algo == fr & T.cr_algo == cr & ...
                            ~isnan(T.fec_snr_db), :);
                    if isempty(sub), continue; end
                    label = sprintf('%s + %s', ...
                        plot_energy_snr.FrLabel.(char(fr)), ...
                        plot_energy_snr.CrLabel.(char(cr)));

                    [~, idx] = min(sub.fec_snr_db);
                    rowsSnr(end+1, :) = {label, ...
                        sub.fec_snr_db(idx), sub.energy_total_fJ(idx), ...
                        sub.energy_fr_fJ(idx), sub.energy_cr_fJ(idx), ...
                        sub.fl_fr(idx), sub.fl_cr(idx), sub.blind_d(idx)}; %#ok<AGROW>

                    [~, idx] = min(sub.energy_total_fJ);
                    rowsErg(end+1, :) = {label, ...
                        sub.fec_snr_db(idx), sub.energy_total_fJ(idx), ...
                        sub.energy_fr_fJ(idx), sub.energy_cr_fJ(idx), ...
                        sub.fl_fr(idx), sub.fl_cr(idx), sub.blind_d(idx)}; %#ok<AGROW>
                end
            end

            colNames = {'pair', 'fec_snr_db', 'energy_total_fJ', ...
                'energy_fr_fJ', 'energy_cr_fJ', 'fl_fr', 'fl_cr', 'blind_d'};
            T_minSnr = cell2table(rowsSnr, 'VariableNames', colNames);
            T_minErg = cell2table(rowsErg, 'VariableNames', colNames);

            fprintf('\n=== Best sensitivity per (FR, CR) pair (min FEC SNR) ===\n');
            disp(T_minSnr);
            fprintf('=== Lowest energy per (FR, CR) pair (min total energy) ===\n');
            disp(T_minErg);

            outDir = fileparts(mfilename('fullpath'));
            writetable(T_minSnr, fullfile(outDir, 'plot_energy_snr_min_fec_snr.csv'));
            writetable(T_minErg, fullfile(outDir, 'plot_energy_snr_min_energy.csv'));
        end

        function test_print_system_energy_tables(testCase)  %#ok<MANU>
            % Minimum system energy per (FR, CR) pair under typical PON
            % link parameters, for K = 1, 10, 100 ONUs and two noise
            % scenarios (shot-noise-limited, amplified ASE+NLI).  System
            % energy is E_tx_wallplug(SNR) + K * E_rx, minimised over the
            % bit-width / BlindD grid.
            dataDir = fileparts(mfilename('fullpath'));
            T = plot_energy_snr.loadCombinedTable(dataDir);

            % --- Typical PON link parameters (1550 nm DP-QPSK at 30.5 GBd)
            lambda_m   = 1550e-9;   % carrier wavelength [m]
            B          = 30.5e9;    % bandwidth = symbol rate (Nyquist) [Hz]
            M_mod      = 4;         % QPSK per polarization
            eta_lsr    = 0.10;      % laser wall-plug efficiency

            % Shot-noise-limited PON (no in-line amplification)
            alpha_dbkm = 0.2;       % SMF loss at 1550 nm [dB/km]
            L_km_pon   = 20;        % typical PON access reach [km]

            % Amplified link (3 × 80 km SMF spans, EDFAs compensating loss)
            n_spans = 3;
            NF_db   = 5;            % EDFA noise figure [dB]
            G_db    = 16;           % compensates 80 km @ 0.2 dB/km
            C_NLI   = 5e23;         % GN-model NLI coeff [W^-2 Hz^2]

            K_values = [1, 10, 100];

            FrAlgos_l = plot_energy_snr.FrAlgos;
            CrAlgos_l = plot_energy_snr.CrAlgos;

            % Pre-compute amplified E_tx per unique SNR — independent of K
            % under the end-of-trunk splitter assumption (eq:snr_nl).
            [uniqSnr, ~, ixAmp] = unique(T.fec_snr_db);
            E_tx_amp_per = nan(size(uniqSnr));
            for u = 1:numel(uniqSnr)
                if isnan(uniqSnr(u)), continue; end
                snr_lin = 10^(uniqSnr(u) / 10);
                try
                    E_tx_amp_per(u) = energy.transmitter_amplified( ...
                        snr_lin, B, lambda_m, n_spans, NF_db, G_db, C_NLI, ...
                        M_mod, eta_lsr);
                catch ME
                    if ~strcmp(ME.identifier, 'energy:transmitter_amplified:noSolution')
                        rethrow(ME);
                    end
                end
            end
            E_tx_amp_fJ = E_tx_amp_per(ixAmp) * 1e15;
            SNR_lin_all = 10.^(T.fec_snr_db / 10);

            rowsShot = {};
            rowsAmp  = {};
            for fi = 1:numel(FrAlgos_l)
                for ci = 1:numel(CrAlgos_l)
                    fr  = FrAlgos_l(fi);
                    cr  = CrAlgos_l(ci);
                    sel = T.fr_algo == fr & T.cr_algo == cr & ~isnan(T.fec_snr_db);
                    sub = T(sel, :);
                    if isempty(sub), continue; end
                    label = sprintf('%s + %s', ...
                        plot_energy_snr.FrLabel.(char(fr)), ...
                        plot_energy_snr.CrLabel.(char(cr)));

                    sub_snr_lin   = SNR_lin_all(sel);
                    sub_E_tx_amp  = E_tx_amp_fJ(sel);

                    for ki = 1:numel(K_values)
                        K = K_values(ki);

                        % Shot noise — vectorised over rows
                        E_tx_shot_fJ = energy.transmitter_shot(sub_snr_lin, ...
                            B, lambda_m, K, alpha_dbkm, L_km_pon, M_mod, eta_lsr) * 1e15;
                        E_sys_shot   = E_tx_shot_fJ + K * sub.energy_total_fJ;
                        rowsShot(end+1, :) = plot_energy_snr.minRow( ...
                            sub, label, K, E_tx_shot_fJ, E_sys_shot); %#ok<AGROW>

                        % Amplified — E_tx independent of K, but K * E_rx isn't
                        E_sys_amp = sub_E_tx_amp + K * sub.energy_total_fJ;
                        rowsAmp(end+1, :) = plot_energy_snr.minRow( ...
                            sub, label, K, sub_E_tx_amp, E_sys_amp); %#ok<AGROW>
                    end
                end
            end

            cols = {'pair', 'K', 'fec_snr_db', 'E_rx_fJ', 'E_tx_fJ', ...
                    'E_sys_fJ', 'fl_fr', 'fl_cr', 'blind_d'};
            T_shot = cell2table(rowsShot, 'VariableNames', cols);
            T_amp  = cell2table(rowsAmp,  'VariableNames', cols);

            fprintf('\n=== Min system energy — shot-noise PON ===\n');
            fprintf('lambda=%g nm, B=%g GBd, fibre %g km @ %g dB/km, eta=%g, M=%d\n', ...
                lambda_m * 1e9, B / 1e9, L_km_pon, alpha_dbkm, eta_lsr, M_mod);
            disp(T_shot);

            fprintf('=== Min system energy — amplified link (ASE + NLI) ===\n');
            fprintf(['lambda=%g nm, B=%g GBd, %d spans, NF=%g dB, G=%g dB, ' ...
                'C_NLI=%.2g W^-2 Hz^2, eta=%g, M=%d\n'], ...
                lambda_m * 1e9, B / 1e9, n_spans, NF_db, G_db, C_NLI, eta_lsr, M_mod);
            disp(T_amp);

            outDir = fileparts(mfilename('fullpath'));
            writetable(T_shot, fullfile(outDir, 'plot_energy_snr_min_sys_shot.csv'));
            writetable(T_amp,  fullfile(outDir, 'plot_energy_snr_min_sys_amp.csv'));
        end

    end

    %% ================================================================
    %  Static helpers
    %% ================================================================
    methods (Static, Access = private)

        function T = loadCombinedTable(dataDir)
            % Loads the grid-sweep table and adds a total-energy column.
            gridFile = fullfile(dataDir, plot_energy_snr.GridSweepFile);
            assert(isfile(gridFile), 'Missing %s — run test_full_grid_sweep first.', gridFile);

            S = load(gridFile, 'tbl');
            T = S.tbl;
            T.energy_total_fJ = T.energy_fr_fJ + T.energy_cr_fJ;
        end

        function row = minRow(sub, label, K, E_tx_vec, E_sys_vec)
            % Argmin over E_sys_vec; returns NaNs if all entries are NaN
            % (e.g. amplified case when SNR exceeds NLI-limited maximum
            % for every config in the pair).
            [Esys_min, idx] = min(E_sys_vec);
            if isnan(Esys_min)
                row = {label, K, NaN, NaN, NaN, NaN, NaN, NaN, NaN};
            else
                row = {label, K, sub.fec_snr_db(idx), ...
                    sub.energy_total_fJ(idx), E_tx_vec(idx), Esys_min, ...
                    sub.fl_fr(idx), sub.fl_cr(idx), sub.blind_d(idx)};
            end
        end

        function keep = paretoMask(x, y)
            % Logical mask for the lower-left Pareto front in (x, y).
            x = x(:);  y = y(:);
            n = numel(x);
            keep = true(n, 1);
            for i = 1:n
                if ~keep(i), continue; end
                dom = (x <= x(i)) & (y <= y(i)) & ((x < x(i)) | (y < y(i)));
                if any(dom)
                    keep(i) = false;
                end
            end
        end

        function [x_pf, y_pf] = paretoFront(x, y)
            % Minimum-energy lower envelope: keep points not dominated
            % on (x = FEC SNR, y = energy) — a point i is dominated if
            % some j satisfies x(j) <= x(i), y(j) <= y(i) with at least
            % one strict inequality.  Output sorted by x ascending.
            x = x(:);  y = y(:);
            n = numel(x);
            keep = true(n, 1);
            for i = 1:n
                if ~keep(i), continue; end
                dom = (x <= x(i)) & (y <= y(i)) & ((x < x(i)) | (y < y(i)));
                if any(dom)
                    keep(i) = false;
                end
            end
            x_pf = x(keep);
            y_pf = y(keep);
            [x_pf, idx] = sort(x_pf);
            y_pf = y_pf(idx);
        end

    end
end
