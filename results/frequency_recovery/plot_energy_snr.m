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

            figure('Name', 'Energy vs FEC SNR — min-energy Pareto front', ...
                'Color', 'w', 'Position', [100 100 900 600]);
            ax = axes;
            hold(ax, 'on');
            grid(ax, 'on');

            SNR_MIN = 7;    % dB
            SNR_MAX = 10;   % dB
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

                    inRange = snr_pf >= SNR_MIN & snr_pf <= SNR_MAX;
                    snr_pf  = snr_pf(inRange);
                    e_pf    = e_pf(inRange);
                    if isempty(snr_pf), continue; end
                    globalMinSnr = min(globalMinSnr, min(snr_pf));

                    label = sprintf('%s + %s', ...
                        plot_energy_snr.FrLabel.(char(fr)), ...
                        plot_energy_snr.CrLabel.(char(cr)));

                    plot(ax, snr_pf, e_pf, ...
                        'LineStyle', '-', 'Marker', 'o', ...
                        'MarkerSize', 4, 'LineWidth', 1.8, ...
                        'Color', colors(k, :), ...
                        'DisplayName', label);
                end
            end

            xlim(ax, [SNR_MIN, SNR_MAX]);
            set(ax, 'YScale', 'log', 'FontSize', 11, 'Box', 'on');
            xlabel(ax, 'FEC SNR threshold  [dB]', 'FontSize', 12);
            ylabel(ax, 'Total energy per bit  [fJ]', 'FontSize', 12);
            title(ax,  'Min-energy bit-width combination per FR \times CR pair', ...
                'FontSize', 12);
            legend(ax, 'Location', 'best', 'FontSize', 10, 'Interpreter', 'none');
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
            n_spans  = 3;
            L_span   = 80;             % span length [km]
            alpha    = 0.2 / 4.343;    % 0.2 dB/km → nepers/km
            beta_2   = -21e-24;        % SMF GVD [s^2/km] (-21 ps^2/km)
            gamma_nl = 1.3;            % SMF non-linear coeff [W^-1 km^-1]
            NF_db    = 5;              % EDFA noise figure [dB]

            K_values = [1, 10, 100];

            FrAlgos_l = plot_energy_snr.FrAlgos;
            CrAlgos_l = plot_energy_snr.CrAlgos;

            % Pre-compute amplified E_tx for every (unique SNR, K).  K
            % here doubles as the number of WDM channels: ASE and NLI in
            % the GN model are evaluated over the full WDM bandwidth
            % B_wdm = K * B_ch, so E_tx now varies with K (unlike the
            % shot-noise case where the splitter ratio cancels).
            [uniqSnr, ~, ixAmp] = unique(T.fec_snr_db);
            E_tx_amp_per = nan(numel(uniqSnr), numel(K_values));
            for ki = 1:numel(K_values)
                B_wdm = K_values(ki) * B;
                for u = 1:numel(uniqSnr)
                    if isnan(uniqSnr(u)), continue; end
                    snr_lin = 10^(uniqSnr(u) / 10);
                    try
                        E_tx_amp_per(u, ki) = energy.transmitter_amplified( ...
                            snr_lin, B_wdm, lambda_m, n_spans, L_span, ...
                            alpha, beta_2, gamma_nl, NF_db, M_mod, eta_lsr);
                    catch ME
                        if ~strcmp(ME.identifier, 'energy:transmitter_amplified:noSolution')
                            rethrow(ME);
                        end
                    end
                end
            end
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

                    sub_snr_lin = SNR_lin_all(sel);
                    sub_ixAmp   = ixAmp(sel);

                    for ki = 1:numel(K_values)
                        K = K_values(ki);

                        % Shot noise — uses per-channel bandwidth B (one
                        % 30.5 GHz Nyquist-shaped signal per ONU after
                        % the 1:K passive split).
                        E_tx_shot_fJ = energy.transmitter_shot(sub_snr_lin, ...
                            B, lambda_m, K, alpha_dbkm, L_km_pon, M_mod, eta_lsr) * 1e15;
                        E_sys_shot   = E_tx_shot_fJ + K * sub.energy_total_fJ;
                        rowsShot(end+1, :) = plot_energy_snr.minRow( ...
                            sub, label, K, E_tx_shot_fJ, E_sys_shot); %#ok<AGROW>

                        % Amplified — uses full WDM bandwidth K*B in the
                        % GN model, so E_tx now scales with K (looked up
                        % from the precomputed (uniqSnr, K) table).
                        sub_E_tx_amp = E_tx_amp_per(sub_ixAmp, ki) * 1e15;
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
            fprintf(['lambda=%g nm, B_ch=%g GBd, %d x %g km spans, ' ...
                'alpha=%g Np/km, beta_2=%.2g s^2/km, gamma=%g W^-1 km^-1, ' ...
                'NF=%g dB, eta=%g, M=%d, B_wdm = K * B_ch\n'], ...
                lambda_m * 1e9, B / 1e9, n_spans, L_span, alpha, beta_2, ...
                gamma_nl, NF_db, eta_lsr, M_mod);
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
