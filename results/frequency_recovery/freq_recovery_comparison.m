classdef freq_recovery_comparison < matlab.unittest.TestCase
%FREQ_RECOVERY_COMPARISON  Compare Tretter-Kay, FFT-search and Fitz
%   frequency estimators over a range of offsets and SNRs.
%
%   For each SNR one figure is produced showing the mean-squared
%   frequency-estimation error (MSE, in kHz^2) versus true offset for
%   all three algorithms.  Each point is averaged over NTrials
%   independent AWGN realisations.
%
%   Run with:
%       results = runtests('freq_recovery_comparison');

    % ================================================================
    %  Parameters
    % ================================================================
    properties (Constant)
        Rs        = 30.504432          % symbol rate [GBd]  (CPON spec)
        NTrials   = 30                 % independent noise trials per point
        SNR_dB_vec  = [5, 10, 15, 20] % SNR sweep [dB]
        % Frequency offset sweep [MHz] — stay within Fitz unambiguity range
        DeltaF_vec  = -1000 : 10 : 1000  % MHz

        % ---- Fixed-point configuration ----------------------------
        FxpConfig_FR  = 'fixed32'    % 'fixed16' | 'fixed32'
        CordicIts     = 16           % CORDIC iterations
        FR_Fitz_N     = 5            % Fitz autocorrelation lag (< TrainingLen=11)
        FR_Nfft       = 512          % FFT size for fft_search (power of 2, >= 11)
        FR_Po2Twiddle = false        % power-of-2 twiddle factors in fft_fxp
    end

    methods (TestClassSetup)
        function buildFxpMex(testCase)
            %BUILDFXPMEX  Compile all three frequency-recovery MEX files once.
            srcDir   = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src');
            buildDir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'build');
            addpath(srcDir);
            addpath(buildDir);

            P = struct();
            P.N_pol         = 2;
            P.TrainingLen   = 11;   % CPON: 11 training symbols per subframe
            P.Rs            = testCase.Rs;
            P.FxpConfig_FR  = testCase.FxpConfig_FR;
            P.CordicIts     = testCase.CordicIts;
            P.FR_Fitz_N     = testCase.FR_Fitz_N;
            P.FR_Nfft       = testCase.FR_Nfft;
            P.FR_Po2Twiddle = testCase.FR_Po2Twiddle;

            cfg = coder.config('mex');
            cfg.GenerateReport     = false;
            cfg.EnableMexProfiling = false;

            fprintf('\n--- Building freq_recovery MEX files (%s) ---\n', ...
                    testCase.FxpConfig_FR);
            build_freq_recovery_tretter_kay_fxp_mex(P, cfg);
            build_freq_recovery_fitz_fxp_mex(P, cfg);
            build_freq_recovery_fft_search_fxp_mex(P, cfg);
            fprintf('--- MEX build complete ---\n\n');
        end
    end

    methods (TestMethodSetup)
        function seedRng(~)
            rng(42);
        end
    end

    % ================================================================
    %  Test
    % ================================================================
    methods (Test)

        function testMSEComparison(testCase)
            %TESTMSECOMPARISON
            %   Sweeps SNR × DeltaF, runs NTrials per point, plots MSE.

            addpath(fullfile(fileparts(mfilename('fullpath')), ...
                            '..', '..', 'src'));

            Rs_      = testCase.Rs;
            DeltaF_v = testCase.DeltaF_vec;
            SNR_v    = testCase.SNR_dB_vec;
            NT       = testCase.NTrials;
            NF       = numel(DeltaF_v);
            NSNR     = numel(SNR_v);

            % Pre-generate one symbol subframe (training fixed; symbols
            % reused across trials — only the AWGN changes)
            [symbols, training] = freq_recovery_comparison.generateSubframe();

            % Fixed-point types and cast training once (same for every trial)
            T_fr        = freq_recovery.fxp_types(testCase.FxpConfig_FR);
            training_fi = cast(training, 'like', T_fr.x);

            % Results: MSE_alg(SNR_idx, DeltaF_idx)
            MSE_tk  = zeros(NSNR, NF);
            MSE_fft = zeros(NSNR, NF);
            MSE_fz  = zeros(NSNR, NF);

            for si = 1:NSNR
                SNR_dB = SNR_v(si);
                fprintf('SNR = %d dB\n', SNR_dB);

                for fi = 1:NF
                    df       = DeltaF_v(fi);
                    true_kHz = df * 1e3;   % MHz -> kHz

                    se_tk  = zeros(NT, 1);
                    se_fft = zeros(NT, 1);
                    se_fz  = zeros(NT, 1);

                    % Apply LO shift once — noise is the only thing
                    % that changes between trials
                    rx_shifted = channel.lo_freq_shift(symbols, df, Rs_, 1);

                    for tr = 1:NT
                        rx   = channel.add_awgn(rx_shifted, SNR_dB);
                        x_fi = cast(rx, 'like', T_fr.x);

                        [~, est_tk]  = freq_recovery.tretter_kay_fxp_mex( ...
                            x_fi, training_fi, Rs_, testCase.CordicIts, T_fr);
                        [~, est_fft] = freq_recovery.fft_search_fxp_mex( ...
                            x_fi, training_fi, Rs_, testCase.FR_Nfft, ...
                            testCase.FR_Po2Twiddle, testCase.CordicIts, T_fr);
                        [~, est_fz]  = freq_recovery.fitz_fxp_mex( ...
                            x_fi, training_fi, Rs_, testCase.FR_Fitz_N, ...
                            testCase.CordicIts, T_fr);

                        se_tk(tr)  = (est_tk  - true_kHz)^2;
                        se_fft(tr) = (est_fft - true_kHz)^2;
                        se_fz(tr)  = (est_fz  - true_kHz)^2;
                    end

                    % Normalise by true offset squared; NaN at zero offset
                    if true_kHz ~= 0
                        norm = true_kHz^2;
                    else
                        norm = NaN;
                    end

                    MSE_tk(si,  fi) = mean(se_tk)  / norm;
                    MSE_fft(si, fi) = mean(se_fft) / norm;
                    MSE_fz(si,  fi) = mean(se_fz)  / norm;
                end
            end

            % ============================================================
            %  Plot — one figure per SNR
            % ============================================================
            colors = lines(3);
            algNames = {'Tretter-Kay', 'FFT search', 'Fitz'};

            for si = 1:NSNR
                figure('Name', sprintf('Freq Recovery MSE | SNR = %d dB', SNR_v(si)), ...
                       'Position', [60 + (si-1)*40, 60 + (si-1)*40, 820, 520], ...
                       'Color', 'w');

                semilogy(DeltaF_v, MSE_tk(si,:),  '-', ...
                    'Color', colors(1,:), 'LineWidth', 1.8, ...
                    'DisplayName', algNames{1});
                hold on;
                semilogy(DeltaF_v, MSE_fft(si,:), '-', ...
                    'Color', colors(2,:), 'LineWidth', 1.8, ...
                    'DisplayName', algNames{2});
                semilogy(DeltaF_v, MSE_fz(si,:),  '-', ...
                    'Color', colors(3,:), 'LineWidth', 1.8, ...
                    'DisplayName', algNames{3});
                hold off;

                grid on;
                set(gca, 'FontSize', 13, 'LineWidth', 1, 'Box', 'on');
                xlabel('True frequency offset [MHz]', 'FontSize', 14);
                ylabel('Normalised MSE [-]',             'FontSize', 14);
                legend('Location', 'best', 'FontSize', 12);
                title(sprintf('Normalised Frequency Estimation MSE  |  SNR = %d dB  |  %d trials', ....
                              SNR_v(si), NT));
            end

            % Print summary table
            fprintf('\n%-14s', 'DeltaF [MHz]');
            for si = 1:NSNR
                fprintf('  SNR=%ddB TK      FFT      Fitz  ', SNR_v(si));
            end
            fprintf('\n');
            for fi = 1:NF
                fprintf('%-14.0f', DeltaF_v(fi));
                for si = 1:NSNR
                    fprintf('  %8.1f  %8.1f  %8.1f  ', ...
                        MSE_tk(si,fi), MSE_fft(si,fi), MSE_fz(si,fi));
                end
                fprintf('\n');
            end
        end

    end

    % ================================================================
    %  Private static helpers
    % ================================================================
    methods (Static, Access = private)

        function [symbols, training] = generateSubframe()
            %GENERATESUBFRAME  Produce one CPON subframe and its training sequence.
            DATA_PER_SUBFRAME = 3586;
            Nbits = DATA_PER_SUBFRAME * 2 * 2;  % 2 pol, 2 bits/sym/pol

            bits = modem.randomBits(Nbits);
            [symbols_full, ~, training, ~] = modem.modulate(bits);

            SUBFRAME_SYMS = 3712;
            symbols = symbols_full(1:SUBFRAME_SYMS, :);  % [3712 x 2]
            % training is [11 x 2], same for every subframe
        end

    end
end
