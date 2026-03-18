classdef test_FreqRecovery < matlab.unittest.TestCase
    % Tests for frequency recovery – floating-point and fixed-point MEX.
    %
    % All tests apply a known frequency offset plus AWGN, then run the
    % algorithm under test.  The constellation before and after correction
    % is plotted for each polarisation, and the post-correction BER is
    % reported and verified against BER_THRESHOLD.
    %
    % MEX compilation
    %   Both fxp MEX binaries are compiled automatically in TestClassSetup
    %   (once per run, before any test method executes).
    %
    % Prerequisites
    %   - MATLAB Coder and Fixed-Point Designer toolboxes must be licensed.
    %   - build_freq_recovery_fft_search_fxp_mex.m and
    %     build_freq_recovery_differential_kay_fxp_mex.m must be on path.

    properties (Constant)
        % ---- Signal -------------------------------------------------
        N_pol       = 2
        Rs          = 30.5          % symbol rate [GBd]
        TrainingLen = 11            % CPON training symbols per subframe

        % ---- Channel ------------------------------------------------
        SNR_dB      = 20            % [dB]  – good SNR to isolate FR errors
        DeltaF_MHz  = 3e3           % [MHz] – frequency offset to apply

        % ---- Float fft_search  --------------------------------------
        FR_FFT_K    = 8             % zero-padding factor

        % ---- Fixed-point config -------------------------------------
        FxpConfig   = 'fixed32'
        CordicIts   = 16
        FR_Nfft     = 128           % FFT size (power of 2 >= TrainingLen)
        FR_Po2Twiddle = false
        FR_BlindD   = 64
        MaxFreq     = 1

        % ---- Pass/fail ----------------------------------------------
        BER_THRESHOLD = 0.05
    end

    % =================================================================
    methods (TestClassSetup)
    % =================================================================

        function seedRng(~)
            rng(42);
        end

        function compileMex(testCase)
            fprintf('  Compiling freq-recovery MEX binaries...\n');
            cfg = coder.config('mex');
            cfg.GenerateReport   = false;
            cfg.IntegrityChecks  = false;
            cfg.ResponsivenessChecks = false;

            % Add src/ to codegen path so +freq_recovery package is found
            srcDir = fullfile(fileparts(mfilename('fullpath')), '..', 'src');
            addpath(srcDir);

            P.N_pol         = testCase.N_pol;
            P.TrainingLen   = testCase.TrainingLen;
            P.Rs            = testCase.Rs;
            P.FR_Nfft       = testCase.FR_Nfft;
            P.FR_Po2Twiddle = testCase.FR_Po2Twiddle;
            P.FR_BlindD      = testCase.FR_BlindD;
            P.FxpConfig_FR  = testCase.FxpConfig;
            P.CordicIts     = testCase.CordicIts;
            P.MaxFreq       = testCase.MaxFreq;

            build_freq_recovery_fft_search_fxp_mex(P, cfg);
            build_freq_recovery_differential_kay_fxp_mex(P, cfg);
            fprintf('  MEX compilation complete.\n');
        end

    end

    % =================================================================
    methods (Test)
    % =================================================================

        % ---- Floating-point fft_search ------------------------------
        function testFFTSearch_Float(testCase)
            [rxSym, frSym, training, txRefBits, deltaF_est] = ...
                runScenario_FFTSearch(testCase);
            BER = computeBER(testCase, frSym, txRefBits);
            fprintf('FFT-search float: delta_f_est = %.3f MHz  BER = %.2e\n', ...
                deltaF_est/1e6, BER);
            plotBeforeAfter(testCase, rxSym, frSym, 'FFT-search (float)', BER);
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD);
        end

        % ---- Floating-point differential_kay ------------------------
        function testDifferentialKay_Float(testCase)
            [rxSym, frSym, training, txRefBits, deltaF_est] = ...
                runScenario_DifferentialKay(testCase);
            BER = computeBER(testCase, frSym, txRefBits);
            fprintf('Differential-Kay float: delta_f_est = %.3f MHz  BER = %.2e\n', ...
                deltaF_est/1e6, BER);
            plotBeforeAfter(testCase, rxSym, frSym, 'Differential-Kay (float)', BER);
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD);
        end

        % ---- Fixed-point fft_search MEX -----------------------------
        function testFFTSearch_Fxp16(testCase)
            [rxSym, frSym, training, txRefBits, deltaF_est] = ...
                runScenarioFxp_FFTSearch(testCase, testCase.FxpConfig);
            BER = computeBER(testCase, frSym, txRefBits);
            fprintf('FFT-search fxp (%s): delta_f_est = %.3f MHz  BER = %.2e\n', ...
                testCase.FxpConfig, deltaF_est/1e6, BER);
            plotBeforeAfter(testCase, rxSym, frSym, ...
                sprintf('FFT-search fxp (%s)', testCase.FxpConfig), BER);
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD);
        end

        % ---- Fixed-point differential_kay MEX -----------------------
        function testDifferentialKay_Fxp16(testCase)
            [rxSym, frSym, training, txRefBits, deltaF_est] = ...
                runScenarioFxp_DifferentialKay(testCase, testCase.FxpConfig);
            BER = computeBER(testCase, frSym, txRefBits);
            fprintf('Differential-Kay fxp (%s): delta_f_est = %.3f MHz  BER = %.2e\n', ...
                testCase.FxpConfig, deltaF_est/1e6, BER);
            plotBeforeAfter(testCase, rxSym, frSym, ...
                sprintf('Differential-Kay fxp (%s)', testCase.FxpConfig), BER);
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD);
        end

    end

    % =================================================================
    methods (Access = private)
    % =================================================================

        % ---- Shared channel builder ---------------------------------
        function [symbols, training, txRefBits, rxSym] = buildChannel(testCase)
            % Modulate random bits to get CPON-framed symbols + training.
            Nbits  = 4 * 3712 * 2;   % two subframes, 2 bits/symbol (QPSK), 2 pol
            txBits = modem.randomBits(Nbits);
            [symbols, ~, training, ~] = modem.modulate(txBits);

            txRefBits = modem.symbolsToBits(symbols);

            % Apply frequency offset then AWGN (no phase noise – isolates FR)
            rxSym = channel.lo_freq_shift(symbols, testCase.DeltaF_MHz, testCase.Rs, 1);
            rxSym = channel.add_awgn(rxSym, testCase.SNR_dB);
        end

        % ---- Float fft_search runner --------------------------------
        function [rxSym, frSym, training, txRefBits, deltaF_est] = ...
                runScenario_FFTSearch(testCase)
            [~, training, txRefBits, rxSym] = buildChannel(testCase);
            [frSym, deltaF_est] = freq_recovery.fft_search( ...
                rxSym, training, testCase.Rs, testCase.FR_FFT_K);
        end

        % ---- Float differential_kay runner --------------------------
        function [rxSym, frSym, training, txRefBits, deltaF_est] = ...
                runScenario_DifferentialKay(testCase)
            [~, training, txRefBits, rxSym] = buildChannel(testCase);
            [frSym, deltaF_est] = freq_recovery.differential_kay( ...
                rxSym, training, testCase.Rs);
        end

        % ---- Fixed-point fft_search MEX runner ----------------------
        function [rxSym, frSym, training, txRefBits, deltaF_est] = ...
                runScenarioFxp_FFTSearch(testCase, config)
            T = freq_recovery.fxp_types(config);
            [~, training, txRefBits, rxSym] = buildChannel(testCase);

            rx_fi       = cast(rxSym,    'like', T.x);
            training_fi = cast(training, 'like', T.x);

            [frSym_fi, deltaF_est] = freq_recovery.fft_search_fxp_mex( ...
                rx_fi, training_fi, testCase.Rs, ...
                double(testCase.FR_Nfft), logical(testCase.FR_Po2Twiddle), ...
                double(testCase.CordicIts), double(testCase.MaxFreq), T, true, 0);

            frSym = double(frSym_fi);
        end

        % ---- Fixed-point differential_kay MEX runner ----------------
        function [rxSym, frSym, training, txRefBits, deltaF_est] = ...
                runScenarioFxp_DifferentialKay(testCase, config)
            T = freq_recovery.fxp_types(config);
            [~, training, txRefBits, rxSym] = buildChannel(testCase);

            rx_fi       = cast(rxSym,    'like', T.x);
            training_fi = cast(training, 'like', T.x);

            [frSym_fi, deltaF_est] = freq_recovery.differential_kay_fxp_mex( ...
                rx_fi, training_fi, testCase.Rs, ...
                double(testCase.CordicIts), T);

            frSym = double(frSym_fi);
        end

        % ---- BER computation ----------------------------------------
        function BER = computeBER(testCase, frSym, txRefBits)
            decidedSyms = modem.decideSymbols(frSym);
            rxBits      = modem.symbolsToBits(decidedSyms);
            nBits       = min(length(txRefBits), length(rxBits));
            nErrors     = sum(txRefBits(1:nBits) ~= rxBits(1:nBits));
            BER         = nErrors / nBits;
        end

        % ---- Constellation plots ------------------------------------
        function plotBeforeAfter(testCase, rxSym, frSym, titleStr, BER)
            figure('Name', titleStr, 'Position', [100 100 1200 500]);
            for p = 1:testCase.N_pol
                subplot(2, 2, (p-1)*2 + 1);
                plot(real(rxSym(:,p)), imag(rxSym(:,p)), '.', 'MarkerSize', 2);
                grid on; axis equal;
                title(sprintf('Before FR  \x2013  Pol %d', p));
                xlabel('In-Phase'); ylabel('Quadrature');

                subplot(2, 2, (p-1)*2 + 2);
                plot(real(frSym(:,p)), imag(frSym(:,p)), '.', 'MarkerSize', 2);
                grid on; axis equal;
                title(sprintf('After FR  \x2013  Pol %d', p));
                xlabel('In-Phase'); ylabel('Quadrature');
            end
            sgtitle(sprintf('QPSK: AWGN + %.1f MHz offset  |  %s  |  BER = %.2e', ...
                testCase.DeltaF_MHz, titleStr, BER));
        end

    end
end
