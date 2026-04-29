classdef bit_width_mse < matlab.unittest.TestCase
%BIT_WIDTH_MSE  MSE of frequency offset estimate vs FR fractional bit width.
%
%   Sweeps the FR fractional bit width (FL) and measures the
%   mean-squared error of the estimated frequency offset (relative to the
%   known transmitted offset DeltaF_Hz) for each FR algorithm:
%
%       - data-aided FFT search       (fft_search_fxp_mex, data_aided=true)
%       - data-aided differential Kay (differential_kay_fxp_mex)
%       - blind FFT search            (D = BlindD, Nfft = FR_Nfft)
%
%   Produces one figure with three lines, one per algorithm.
%
%   Parallelism mirrors bit_width.m: each FL iteration builds a private
%   FR MEX into a per-FL temp directory that mirrors the +freq_recovery
%   package layout.  A parfor loop then runs Monte-Carlo trials
%   concurrently — each worker addpath-es its own dir so it uses the MEX
%   built for that specific FL value.
%
%   Carrier recovery is not exercised here, so no CR MEX is built.
%
%   Run with:
%       runtests('bit_width_mse')

    %% ================================================================
    %  Parameters
    %% ================================================================
    properties (Constant)

        % System
        Rs          = 30.5          % symbol rate [GBd]
        N_pol       = 2
        TrainingLen = 11

        % Monte-Carlo
        NTrials     = 50            % trials per FL point

        % Single operating SNR (high enough that AWGN does not dominate
        % over quantisation effects, but realistic for the FEC region)
        SNR_dB      = 15            % [dB]

        % Bit width sweep — integer bits fixed; WL = IntBits + FL
        IntBits     = 16
        FL_vec      = [2, 4, 6, 8, 10, 12, 14, 16]

        % Channel conditions
        DeltaF_Hz   = 2e9           % true frequency offset [Hz]  (MSE truth)
        LW_Hz       = 1000e3        % laser linewidth [Hz]

        % FFT search parameters — Nfft constant for all runs
        FR_Nfft       = 512
        FR_Po2Twiddle = false
        MaxFreq       = 0.1

        % Blind FFT search observation length
        BlindD = 512

        % CORDIC iterations
        CordicIts = 16

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

        function setupBuildPath(~)
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'build'));
        end

    end

    %% ================================================================
    %  Tests
    %% ================================================================
    methods (Test)

        function test_mse_vs_fr_bitwidth(testCase)
            P    = testCase;
            NFL  = length(P.FL_vec);
            Prms = bit_width_mse.extractParams(testCase);

            % ---- Phase 1: Serial MEX builds into per-FL temp dirs -------
            frDirs = cell(1, NFL);
            for bi = 1:NFL
                fl = P.FL_vec(bi);
                fprintf('[FR FL = %2d] building FR MEX  (%d / %d)\n', fl, bi, NFL);
                frDirs{bi} = bit_width_mse.buildFRMex(P, struct('WL', P.IntBits + fl, 'FL', fl));
            end

            % ---- Phase 2: Parallel FL sweeps ----------------------------
            mse_fft_DA    = nan(NFL, 1);
            mse_dk_DA     = nan(NFL, 1);
            mse_fft_blind = nan(NFL, 1);

            FL_vec_b   = P.FL_vec;
            IntBits_b  = P.IntBits;
            BlindD_b   = P.BlindD;

            parfor pi = 1:NFL
                addpath(frDirs{pi});

                fl     = FL_vec_b(pi);
                T_fr_w = freq_recovery.fxp_types(struct('WL', IntBits_b + fl, 'FL', fl));

                mse_fft_DA(pi)    = bit_width_mse.runMseTrials(Prms, 'fft_search',       0,        T_fr_w);
                mse_dk_DA(pi)     = bit_width_mse.runMseTrials(Prms, 'differential_kay', 0,        T_fr_w);
                mse_fft_blind(pi) = bit_width_mse.runMseTrials(Prms, 'fft_search_blind', BlindD_b, T_fr_w);
            end

            bit_width_mse.plotMSE(P, P.FL_vec, mse_fft_DA, mse_dk_DA, mse_fft_blind);

            cellfun(@(d) rmdir(d, 's'), frDirs, 'UniformOutput', false);
        end

    end

    %% ================================================================
    %  Static helpers
    %% ================================================================
    methods (Static, Access = private)

        function tempDir = buildFRMex(P, fxp_fr)
            % Clears all loaded MEX, builds FR MEX to the default src/
            % location, then copies the result into a fresh temp directory
            % that mirrors the +freq_recovery package structure.
            clear mex %#ok<CLMEX>

            B.Rs            = P.Rs;
            B.N_pol         = P.N_pol;
            B.TrainingLen   = P.TrainingLen;
            B.FR_Nfft       = P.FR_Nfft;
            B.FR_Po2Twiddle = P.FR_Po2Twiddle;
            B.FR_BlindD     = P.BlindD;
            B.FxpConfig_FR  = fxp_fr;
            B.CordicIts     = P.CordicIts;
            B.MaxFreq       = P.MaxFreq;

            cfg = coder.config('mex');
            cfg.GenerateReport = false;

            build_freq_recovery_fft_search_fxp_mex(B, cfg);
            build_freq_recovery_differential_kay_fxp_mex(B, cfg);

            srcDir  = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src');
            tempDir = tempname;
            dstPkg  = fullfile(tempDir, '+freq_recovery');
            mkdir(dstPkg);
            ext = mexext;
            for f = {'fft_search_fxp_mex', 'differential_kay_fxp_mex'}
                copyfile( ...
                    fullfile(srcDir, '+freq_recovery', [f{1} '.' ext]), ...
                    fullfile(dstPkg,                   [f{1} '.' ext]));
            end
        end

        function Params = extractParams(testCase)
            % Plain struct so it can be broadcast into parfor without
            % serialising the full TestCase object.
            Params.Rs            = testCase.Rs;
            Params.N_pol         = testCase.N_pol;
            Params.NTrials       = testCase.NTrials;
            Params.SNR_dB        = testCase.SNR_dB;
            Params.DeltaF_Hz     = testCase.DeltaF_Hz;
            Params.LW_Hz         = testCase.LW_Hz;
            Params.FR_Nfft       = testCase.FR_Nfft;
            Params.FR_Po2Twiddle = testCase.FR_Po2Twiddle;
            Params.MaxFreq       = testCase.MaxFreq;
            Params.CordicIts     = testCase.CordicIts;
        end

        function nmse = runMseTrials(Params, fr_algo, blindD, T_fr)
            % Normalised MSE: ((f_est - DeltaF) / DeltaF)^2 averaged over trials.
            sqErr = zeros(Params.NTrials, 1);
            for tr = 1:Params.NTrials
                f_est     = bit_width_mse.estimateFreq(Params, fr_algo, blindD, T_fr);
                sqErr(tr) = ((f_est - Params.DeltaF_Hz) / Params.DeltaF_Hz)^2;
            end
            nmse = mean(sqErr);
        end

        function f_est = estimateFreq(P, fr_algo, blindD, T_fr)
            BITS_PER_SF = 3586 * 2 * 2;

            txBits = modem.randomBits(BITS_PER_SF);
            [symbols, ~, training, ~] = modem.modulate(txBits);

            rx = channel.lo_freq_shift(symbols, P.DeltaF_Hz / 1e6, P.Rs, 1);
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
                    error('bit_width_mse:unknownFR', 'Unknown FR algorithm: %s', fr_algo);
            end

            f_est = double(f_est);
        end

        function plotMSE(P, fl_vec, mse_fft_DA, mse_dk_DA, mse_fft_blind)
            frLabels = {'FFT Search (data-aided)', 'Diff. Kay (data-aided)', ...
                        sprintf('FFT Search (blind D=%d)', P.BlindD)};
            mseAll   = {mse_fft_DA, mse_dk_DA, mse_fft_blind};
            colors   = lines(3);
            markers  = {'o', 's', '^'};

            figure('Name', 'FR normalised MSE vs bit width', ...
                'Position', [100, 100, 750, 500], 'Color', 'w');
            ax = axes;
            hold(ax, 'on');
            grid(ax, 'on');

            for fr = 1:3
                valid = ~isnan(mseAll{fr}) & mseAll{fr} > 0;
                if any(valid)
                    plot(ax, fl_vec(valid), mseAll{fr}(valid), ...
                        'LineStyle', '-', 'Marker', markers{fr}, ...
                        'MarkerSize', 6, 'LineWidth', 1.8, ...
                        'Color', colors(fr, :), ...
                        'DisplayName', frLabels{fr});
                end
            end

            set(ax, 'YScale', 'log', 'FontSize', 11, 'Box', 'on', ...
                'XTick', fl_vec, 'XLim', [fl_vec(1) - 1, fl_vec(end) + 1]);
            xlabel(ax, 'FR fractional bit width  [bits]', 'FontSize', 12);
            ylabel(ax, 'Normalised MSE  E[((f_{est} - \Deltaf) / \Deltaf)^2]', 'FontSize', 12);
            title(ax, sprintf(['Frequency offset NMSE  |  IntBits = %d,  SNR = %d dB\n' ...
                '\\DeltaF = %.0f MHz,  LW = %.0f kHz'], ...
                P.IntBits, P.SNR_dB, P.DeltaF_Hz/1e6, P.LW_Hz/1e3), 'FontSize', 11);
            legend(ax, 'Location', 'northeast', 'FontSize', 10);
        end

    end
end
