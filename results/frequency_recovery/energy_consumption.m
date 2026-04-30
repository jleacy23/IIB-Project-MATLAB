classdef energy_consumption < matlab.unittest.TestCase
%ENERGY_CONSUMPTION  Energy per bit vs fractional bit width.
%
%   Two test methods, each producing one figure:
%
%   test_fr_energy_vs_bitwidth — energy per bit for the three FR variants
%       used in bit_width.m (FFT search data-aided, Differential Kay
%       data-aided, FFT search blind D = BlindD), plotted against
%       fractional bit width.
%
%   test_cr_energy_vs_bitwidth — energy per bit for the two CR variants
%       (Viterbi-Viterbi, pilots-only), plotted against fractional bit
%       width.
%
%   Operation counts come from tab:fft_cost, tab:diffkay_cost,
%   tab:viterbi_cost and tab:pilot_cost in
%   report/frequency_recovery/frequency_recovery.tex.  The per-bit
%   energy conversion is performed by src/+energy/receiver.m.  n is
%   taken as the fractional bit width only (integer bits ignored).
%
%   Run with:
%       runtests('energy_consumption')
%       runtests('energy_consumption', 'ProcedureName', 'test_fr_energy_vs_bitwidth')

    %% ================================================================
    %  Parameters
    %% ================================================================
    properties (Constant)

        % Algorithm parameters (match bit_width.m)
        TrainingLen = 11
        BlindD      = 512
        FR_Nfft     = 512
        BlockLen    = 32

        % Modulation
        M            = 4   % QPSK
        Oversampling = 1

        % CPON subframe length: FR estimate amortised across this many symbols
        SubframeLen = 3712

        % Bit width sweep (n = FL only; integer bits ignored)
        FL_vec = [2, 4, 6, 8, 10, 12, 14, 16]

        % Per-bit energy coefficients (fJ) — fits to horowitz2014computing
        % E_A(n) = EAdd_fJ * n,  E_M(n) = EMult_fJ * n^2
        EAdd_fJ  = 3.16
        EMult_fJ = 3.03

    end

    %% ================================================================
    %  Test class setup
    %% ================================================================
    methods (TestClassSetup)

        function setupPath(~)
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src'));
        end

    end

    %% ================================================================
    %  Tests
    %% ================================================================
    methods (Test)

        function test_fr_energy_vs_bitwidth(testCase)
            P = testCase;

            [NM_fft_DA,    NA_fft_DA]    = energy_consumption.fft_search_counts(P.TrainingLen, P.FR_Nfft, false);
            [NM_fft_blind, NA_fft_blind] = energy_consumption.fft_search_counts(P.BlindD,      P.FR_Nfft, true);
            [NM_dk,        NA_dk]        = energy_consumption.differential_kay_counts(P.TrainingLen);

            % FR runs once per subframe; amortise across SubframeLen symbols.
            NM_fft_DA    = NM_fft_DA    / P.SubframeLen;
            NA_fft_DA    = NA_fft_DA    / P.SubframeLen;
            NM_fft_blind = NM_fft_blind / P.SubframeLen;
            NA_fft_blind = NA_fft_blind / P.SubframeLen;
            NM_dk        = NM_dk        / P.SubframeLen;
            NA_dk        = NA_dk        / P.SubframeLen;

            E_fft_DA    = energy_consumption.energySweep(NA_fft_DA,    NM_fft_DA,    P);
            E_dk        = energy_consumption.energySweep(NA_dk,        NM_dk,        P);
            E_fft_blind = energy_consumption.energySweep(NA_fft_blind, NM_fft_blind, P);

            colors = lines(3);
            figure('Name', 'FR energy vs bit width', 'Color', 'w', 'Position', [100 100 750 500]);
            ax = axes; hold(ax, 'on'); grid(ax, 'on');
            plot(ax, P.FL_vec, E_fft_DA,    '-o', 'LineWidth', 1.8, 'MarkerSize', 6, ...
                'Color', colors(1,:), 'DisplayName', 'FFT Search (data-aided)');
            plot(ax, P.FL_vec, E_dk,        '-s', 'LineWidth', 1.8, 'MarkerSize', 6, ...
                'Color', colors(2,:), 'DisplayName', 'Diff. Kay (data-aided)');
            plot(ax, P.FL_vec, E_fft_blind, '-^', 'LineWidth', 1.8, 'MarkerSize', 6, ...
                'Color', colors(3,:), 'DisplayName', sprintf('FFT Search (blind D=%d)', P.BlindD));
            set(ax, 'YScale', 'log', 'FontSize', 11, 'Box', 'on', ...
                'XTick', P.FL_vec, 'XLim', [P.FL_vec(1)-1, P.FL_vec(end)+1]);
            xlabel(ax, 'Fractional bit width  n  [bits]', 'FontSize', 12);
            ylabel(ax, 'Energy per bit  [fJ]',           'FontSize', 12);
            title(ax,  'Frequency recovery: energy per bit vs bit width', 'FontSize', 11);
            legend(ax, 'Location', 'northwest', 'FontSize', 10);
        end

        function test_cr_energy_vs_bitwidth(testCase)
            P = testCase;

            [NM_vv, NA_vv] = energy_consumption.viterbi_counts(P.BlockLen);
            [NM_po, NA_po] = energy_consumption.pilots_only_counts();

            % CR runs once per BlockLen-symbol block; amortise across the block.
            NM_vv = NM_vv / P.BlockLen;
            NA_vv = NA_vv / P.BlockLen;
            NM_po = NM_po / P.BlockLen;
            NA_po = NA_po / P.BlockLen;

            E_vv = energy_consumption.energySweep(NA_vv, NM_vv, P);
            E_po = energy_consumption.energySweep(NA_po, NM_po, P);

            colors = lines(2);
            figure('Name', 'CR energy vs bit width', 'Color', 'w', 'Position', [180 100 750 500]);
            ax = axes; hold(ax, 'on'); grid(ax, 'on');
            plot(ax, P.FL_vec, E_vv, '-o', 'LineWidth', 1.8, 'MarkerSize', 6, ...
                'Color', colors(1,:), 'DisplayName', 'Viterbi-Viterbi');
            plot(ax, P.FL_vec, E_po, '-s', 'LineWidth', 1.8, 'MarkerSize', 6, ...
                'Color', colors(2,:), 'DisplayName', 'Pilots-only');
            set(ax, 'YScale', 'log', 'FontSize', 11, 'Box', 'on', ...
                'XTick', P.FL_vec, 'XLim', [P.FL_vec(1)-1, P.FL_vec(end)+1]);
            xlabel(ax, 'Fractional bit width  n  [bits]', 'FontSize', 12);
            ylabel(ax, 'Energy per bit  [fJ]',           'FontSize', 12);
            title(ax,  'Carrier recovery: energy per bit vs bit width', 'FontSize', 11);
            legend(ax, 'Location', 'northwest', 'FontSize', 10);
        end

    end

    %% ================================================================
    %  Static helpers — operation counts from frequency_recovery.tex
    %% ================================================================
    methods (Static, Access = private)

        function [NM, NA] = fft_search_counts(L, Nfft, blind)
            % tab:fft_cost — total real ops per CFO estimate (single pol).
            % Blind operation uses the (*) values for forming z[i].
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
            % tab:diffkay_cost — total real ops per CFO estimate (single pol).
            % Sum of: forming arg(z) (twice), differential phase, first-pass
            % correction, and the Kay weighted sum.
            NM = 7 * L + 1;
            NA = 4 * L - 2;
        end

        function [NM, NA] = viterbi_counts(N)
            % tab:viterbi_cost — per block of N symbols (single pol).
            NM = 16 * N + 1;            % 12N + 4N + 1
            NA = 14 * N - 1;            % 9N + 3N + 2(N-1) + 1
        end

        function [NM, NA] = pilots_only_counts()
            % tab:pilot_cost — per block (single pol).
            NM = 1;
            NA = 1;
        end

        function E = energySweep(NA, NM, P)
            E = arrayfun( ...
                @(n) energy.receiver(NA, NM, P.EAdd_fJ, P.EMult_fJ, P.M, P.Oversampling, n), ...
                P.FL_vec);
        end

    end
end
