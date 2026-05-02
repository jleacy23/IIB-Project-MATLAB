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

            set(ax, 'YScale', 'log', 'FontSize', 11, 'Box', 'on');
            xlabel(ax, 'FEC SNR threshold  [dB]', 'FontSize', 12);
            ylabel(ax, 'Total energy per bit  [fJ]', 'FontSize', 12);
            title(ax,  'Min-energy bit-width combination per FR \times CR pair', ...
                'FontSize', 12);
            legend(ax, 'Location', 'best', 'FontSize', 10, 'Interpreter', 'none');
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
