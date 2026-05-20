classdef max_cfo < matlab.unittest.TestCase
%MAX_CFO  Normalised MSE of CFO estimate vs true CFO across full Nyquist range.
%
%   Sweeps the true carrier frequency offset from a small fraction of the
%   symbol rate up to the symbol rate itself (normalised CFO = 1) and
%   measures the MSE of the estimate, normalised to the symbol rate, for
%   each FR algorithm:
%       - data-aided FFT search       (fft_search_fxp_mex, data_aided=true)
%       - data-aided differential Kay (differential_kay_fxp_mex)
%       - blind FFT search            (D = BlindD = 512)
%
%   Bit width is held fixed at 16 integer + 16 fractional bits and
%   max_freq = 1 is passed to every estimator so the search range covers
%   the full normalised range [-1, 1] x R_s.
%
%   Produces one figure with three lines.
%
%   Run with:
%       runtests('max_cfo')

    %% ================================================================
    %  Parameters
    %% ================================================================
    properties (Constant)

        % System
        Rs          = 30.5              % symbol rate [GBd]
        N_pol       = 2
        TrainingLen = 11

        % Monte-Carlo
        NTrials     = 100                % trials per CFO point

        % Fixed-point configuration (fixed for all runs)
        FxpConfig_FR = struct('WL', 32, 'FL', 16)   % 16 int + 16 frac
        CordicIts    = 16

        % Operating point
        SNR_dB = 15                     % [dB]
        LW_Hz  = 1000e3                 % laser linewidth [Hz]

        % CFO sweep — normalised to R_s, 0 to <0.5 in steps of 0.05
        NormCFO_vec = 0 : 0.025 : 0.475

        % FFT search parameters
        FR_Nfft       = 512
        FR_Po2Twiddle = false
        MaxFreq       = 1               % full normalised search range

        % Blind FFT search observation length
        BlindD = 512

        % Build control
        Rebuild = true

    end

    %% ================================================================
    %  Test class setup
    %% ================================================================
    methods (TestClassSetup)

        function setupPath(~)
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src'));
        end

        function seedRng(~)
            rng(42);
        end

        function buildMex(testCase)
            buildDir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'build');
            addpath(buildDir);

            B.Rs            = testCase.Rs;
            B.N_pol         = testCase.N_pol;
            B.TrainingLen   = testCase.TrainingLen;
            B.FR_Nfft       = testCase.FR_Nfft;
            B.FR_Po2Twiddle = testCase.FR_Po2Twiddle;
            B.FR_BlindD     = testCase.BlindD;
            B.FxpConfig_FR  = testCase.FxpConfig_FR;
            B.CordicIts     = testCase.CordicIts;
            B.MaxFreq       = testCase.MaxFreq;

            cfg = coder.config('mex');
            cfg.GenerateReport = false;

            if testCase.Rebuild
                fprintf('Building FR MEX objects...\n');
                build_freq_recovery_fft_search_fxp_mex(B, cfg);
                build_freq_recovery_differential_kay_fxp_mex(B, cfg);
                fprintf('FR MEX built.\n');
            end
        end

    end

    %% ================================================================
    %  Tests
    %% ================================================================
    methods (Test)

        function test_nmse_vs_max_cfo(testCase)
            P  = testCase;
            NC = length(P.NormCFO_vec);

            nmse_fft_DA    = nan(NC, 1);
            nmse_dk_DA     = nan(NC, 1);
            nmse_fft_blind = nan(NC, 1);

            T_fr = freq_recovery.fxp_types(P.FxpConfig_FR);

            for ci = 1:NC
                normCFO = P.NormCFO_vec(ci);
                fprintf('[%2d/%2d] normalised CFO = %.3f\n', ci, NC, normCFO);

                nmse_fft_DA(ci)    = max_cfo.runTrials(P, 'fft_search',       0,        normCFO, T_fr);
                nmse_dk_DA(ci)     = max_cfo.runTrials(P, 'differential_kay', 0,        normCFO, T_fr);
                nmse_fft_blind(ci) = max_cfo.runTrials(P, 'fft_search_blind', P.BlindD, normCFO, T_fr);
            end

            max_cfo.plotNMSE(P, P.NormCFO_vec, nmse_fft_DA, nmse_dk_DA, nmse_fft_blind);
        end

    end

    %% ================================================================
    %  Static helpers
    %% ================================================================
    methods (Static, Access = private)

        function nmse = runTrials(P, fr_algo, blindD, normCFO, T_fr)
            Rs_Hz  = P.Rs * 1e9;
            f_true = normCFO * Rs_Hz;

            sqErr = zeros(P.NTrials, 1);
            for tr = 1:P.NTrials
                f_est     = max_cfo.estimateFreq(P, fr_algo, blindD, f_true, T_fr);
                sqErr(tr) = ((f_est - f_true) / Rs_Hz)^2;
            end
            nmse = mean(sqErr);
        end

        function f_est = estimateFreq(P, fr_algo, blindD, f_true_Hz, T_fr)
            BITS_PER_SF = 3586 * 2 * 2;

            txBits = modem.randomBits(BITS_PER_SF);
            [symbols, ~, training, ~] = modem.modulate(txBits);

            rx = channel.lo_freq_shift(symbols, f_true_Hz / 1e6, P.Rs, 1);
            rx = channel.add_awgn(rx, P.SNR_dB);
            rx = channel.add_phase_noise(rx, P.Rs, P.LW_Hz);

            rx_fi = cast(rx,       'like', T_fr.x);
            tr_fi = cast(training, 'like', T_fr.x);

            switch fr_algo
                case 'fft_search'
                    [~, f_est] = freq_recovery.fft_search_fxp_mex( ...
                        rx_fi, tr_fi, P.Rs, P.FR_Nfft, P.FR_Po2Twiddle, ...
                        P.CordicIts, P.MaxFreq, T_fr, true, 0);
                case 'fft_search_blind'
                    [~, f_est] = freq_recovery.fft_search_fxp_mex( ...
                        rx_fi, tr_fi, P.Rs, P.FR_Nfft, P.FR_Po2Twiddle, ...
                        P.CordicIts, P.MaxFreq, T_fr, false, blindD);
                case 'differential_kay'
                    [~, f_est] = freq_recovery.differential_kay_fxp_mex( ...
                        rx_fi, tr_fi, P.Rs, P.CordicIts, T_fr, true, 0, P.MaxFreq);
                otherwise
                    error('max_cfo:unknownFR', 'Unknown FR algorithm: %s', fr_algo);
            end

            f_est = double(f_est);
        end

        function plotNMSE(P, normCFO_vec, nmse_fft_DA, nmse_dk_DA, nmse_fft_blind)
            frLabels = {'FFT Search (data-aided)', 'Diff. Kay (data-aided)', ...
                        sprintf('FFT Search (blind D=%d)', P.BlindD)};
            nmseAll  = {nmse_fft_DA, nmse_dk_DA, nmse_fft_blind};
            colors   = lines(3);
            markers  = {'o', 's', '^'};

            figure('Name', 'FR NMSE vs normalised CFO', ...
                'Position', [100, 100, 750, 500], 'Color', 'w');
            ax = axes;
            hold(ax, 'on');
            grid(ax, 'on');

            for fr = 1:3
                valid = ~isnan(nmseAll{fr}) & nmseAll{fr} > 0;
                if any(valid)
                    plot(ax, normCFO_vec(valid), nmseAll{fr}(valid), ...
                        'LineStyle', '-', 'Marker', markers{fr}, ...
                        'MarkerSize', 6, 'LineWidth', 1.8, ...
                        'Color', colors(fr, :), ...
                        'DisplayName', frLabels{fr});
                end
            end

            set(ax, 'YScale', 'log', 'FontSize', 11, 'Box', 'on');
            xlim(ax, [0, 0.45]);
            xlabel(ax, 'Normalised CFO  \Deltaf / R_s', 'FontSize', 12);
            ylabel(ax, 'Normalised MSE  E[((f_{est} - \Deltaf) / R_s)^2]', 'FontSize', 12);
            title(ax, sprintf(['CFO estimator NMSE  |  WL = %d, FL = %d,  SNR = %d dB\n' ...
                'LW = %.0f kHz,  max\\_freq = %g'], ...
                P.FxpConfig_FR.WL, P.FxpConfig_FR.FL, P.SNR_dB, ...
                P.LW_Hz/1e3, P.MaxFreq), 'FontSize', 11);
            legend(ax, 'Location', 'northwest', 'FontSize', 10);
        end

    end
end
