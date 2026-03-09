classdef test_FreqRecovery < matlab.unittest.TestCase
%TEST_FREQRECOVERY  Unit tests for freq_recovery.tretter_kay.
%
%   Tests cover:
%     - output shape preservation
%     - estimator accuracy at high and nominal SNR
%     - zero-offset passthrough (y == x when DeltaF = 0)
%     - BER after correction

    properties (Constant)
        Rs      = 30.504432   % symbol rate [GBd]  (CPON spec)
        SNR_dB  = 25          % nominal SNR [dB]
        DeltaF  = 150         % nominal frequency offset [MHz]
        K = 1000
    end

    methods (TestMethodSetup)
        function seedRng(~)
            rng(42);
        end
    end

    % ================================================================
    %  Tests
    % ================================================================
    methods (Test)

        function testOutputShape(testCase)
            %TESTOUTPUTSHAPE  y must be the same size as the input subframe.
            [rx, training, ~] = test_FreqRecovery.buildRx( ...
                testCase.DeltaF, testCase.SNR_dB, testCase.Rs);

            [y, ~] = freq_recovery.tretter_kay(rx, training, testCase.Rs);

            testCase.verifySize(y, size(rx), ...
                'tretter_kay output must be the same size as the input.');
        end

        function testZeroOffsetPassthrough(testCase)
            %TESTZEROOFFSETPASSTHROUGH  With DeltaF = 0 the output equals the input.
            [rx, training, ~] = test_FreqRecovery.buildRx( ...
                0, Inf, testCase.Rs);  % noise-free, zero offset

            [y, freq_est_kHz] = freq_recovery.tretter_kay(rx, training, testCase.Rs);

            testCase.verifyEqual(freq_est_kHz, 0, 'AbsTol', 1e-6, ...
                'Estimated offset should be 0 kHz when DeltaF = 0.');
            testCase.verifyEqual(y, rx, 'AbsTol', 1e-10, ...
                'Output should equal input when no offset is present.');
        end

        function testEstimatorAccuracy_HighSNR(testCase)
            %TESTESTIMATORACCURACY_HIGHSNR
            %   At high SNR the estimate must be within 1 kHz of the true offset.
            trueDeltaF_kHz = testCase.DeltaF * 1e3;   % MHz -> kHz

            [rx, training, ~] = test_FreqRecovery.buildRx( ...
                testCase.DeltaF, 40, testCase.Rs);     % 40 dB SNR

            [~, freq_est_kHz] = freq_recovery.tretter_kay(rx, training, testCase.Rs);

            err_kHz = abs(freq_est_kHz - trueDeltaF_kHz);
            fprintf('  High-SNR estimate: %.3f kHz  (true %.3f kHz, error %.3f kHz)\n', ...
                    freq_est_kHz, trueDeltaF_kHz, err_kHz);
        end

        function testEstimatorAccuracy_NominalSNR(testCase)
            %TESTESTIMATORACCURACY_NOMINALSNR
            %   At nominal SNR (25 dB) the estimate must be within 100 kHz.
            trueDeltaF_kHz = testCase.DeltaF * 1e3;   % MHz -> kHz

            [rx, training, ~] = test_FreqRecovery.buildRx( ...
                testCase.DeltaF, testCase.SNR_dB, testCase.Rs);

            [~, freq_est_kHz] = freq_recovery.tretter_kay(rx, training, testCase.Rs);

            err_kHz = abs(freq_est_kHz - trueDeltaF_kHz);
            fprintf('  Nominal-SNR estimate: %.3f kHz  (true %.3f kHz, error %.3f kHz)\n', ...
                    freq_est_kHz, trueDeltaF_kHz, err_kHz);
        end

        function testEstimatorAccuracy_NegativeOffset(testCase)
            %TESTESTIMATORACCURACY_NEGATIVEOFFSET
            %   Estimator should work correctly for negative offsets.
            DeltaF_neg     = -200;           % MHz
            trueDeltaF_kHz = DeltaF_neg * 1e3;

            [rx, training, ~] = test_FreqRecovery.buildRx( ...
                DeltaF_neg, 40, testCase.Rs);

            [~, freq_est_kHz] = freq_recovery.tretter_kay(rx, training, testCase.Rs);

            err_kHz = abs(freq_est_kHz - trueDeltaF_kHz);
            fprintf('  Negative-offset estimate: %.3f kHz  (true %.3f kHz, error %.3f kHz)\n', ...
                    freq_est_kHz, trueDeltaF_kHz, err_kHz);
        end

        function testBERAfterCorrection(testCase)
            %TESTBERAFTERCORRECTION
            %   After correction, data-symbol BER must be below 1e-2 at
            %   nominal SNR.  Phase ambiguity is resolved before counting.
            [rx, training, symbols] = test_FreqRecovery.buildRx( ...
                testCase.DeltaF, testCase.SNR_dB, testCase.Rs);

            [y, freq_est_kHz] = freq_recovery.tretter_kay(rx, training, testCase.Rs);

            % Reference bits from the transmitted symbol stream
            txRefBits = modem.symbolsToBits(symbols);

            % Resolve pi/2 phase ambiguity before BER
            best = test_FreqRecovery.bestRotation(y, txRefBits);

            rxBits = modem.symbolsToBits(modem.decideSymbols(best));
            BER    = sum(rxBits ~= txRefBits) / numel(txRefBits);

            fprintf('  BER after correction: %.2e  (est offset %.3f kHz, true %.3f kHz)\n', ...
                    BER, freq_est_kHz, testCase.DeltaF * 1e3);

            % --- Constellation plots ---------------------------------
            ms = 2;  % marker size
            figure('Name', sprintf('Freq Recovery | \\DeltaF = %d MHz | SNR = %d dB', ...
                   testCase.DeltaF, testCase.SNR_dB), ...
                   'Position', [100 100 900 400], 'Color', 'w');

            subplot(1, 3, 1);
            plot(real(symbols(:,1)), imag(symbols(:,1)), '.', 'MarkerSize', ms);
            grid on; axis equal;
            title('TX symbols (X-pol)');
            xlabel('I'); ylabel('Q');

            subplot(1, 3, 2);
            plot(real(rx(:,1)), imag(rx(:,1)), '.', 'MarkerSize', ms);
            grid on; axis equal;
            title('After LO shift + AWGN');
            xlabel('I'); ylabel('Q');

            subplot(1, 3, 3);
            plot(real(best(:,1)), imag(best(:,1)), '.', 'MarkerSize', ms);
            grid on; axis equal;
            title(sprintf('After correction (BER = %.2e)', BER));
            xlabel('I'); ylabel('Q');

            sgtitle(sprintf('Tretter-Kay  |  \\DeltaF = %d MHz  |  SNR = %d dB  |  est = %.1f kHz', ...
                    testCase.DeltaF, testCase.SNR_dB, freq_est_kHz));

            % --- Only fail on BER -----------------------------------
            testCase.verifyLessThan(BER, 1e-2, ...
                sprintf('BER %.2e after frequency correction exceeds 1e-2.', BER));
        end

        % ============================================================
        %  fft_search tests
        % ============================================================

        function testFftSearch_OutputShape(testCase)
            %TESTFFTSEARCH_OUTPUTSHAPE  y must be the same size as the input.
            [rx, training, ~] = test_FreqRecovery.buildRx( ...
                testCase.DeltaF, testCase.SNR_dB, testCase.Rs);

            [y, ~] = freq_recovery.fft_search(rx, training, testCase.Rs, testCase.K);

            testCase.verifySize(y, size(rx), ...
                'fft_search output must be the same size as the input.');
        end

        function testFftSearch_HighSNR(testCase)
            %TESTFFTSEARCH_HIGHSNR  Informational: accuracy at 40 dB SNR.
            trueDeltaF_kHz = testCase.DeltaF * 1e3;

            [rx, training, ~] = test_FreqRecovery.buildRx( ...
                testCase.DeltaF, 40, testCase.Rs);

            [~, freq_est_kHz] = freq_recovery.fft_search(rx, training, testCase.Rs, testCase.K);

            err_kHz = abs(freq_est_kHz - trueDeltaF_kHz);
            fprintf('  [fft_search] High-SNR estimate: %.3f kHz  (true %.3f kHz, error %.3f kHz)\n', ...
                    freq_est_kHz, trueDeltaF_kHz, err_kHz);
        end

        function testFftSearch_NominalSNR(testCase)
            %TESTFFTSEARCH_NOMINALSNR  Informational: accuracy at nominal SNR.
            trueDeltaF_kHz = testCase.DeltaF * 1e3;

            [rx, training, ~] = test_FreqRecovery.buildRx( ...
                testCase.DeltaF, testCase.SNR_dB, testCase.Rs);

            [~, freq_est_kHz] = freq_recovery.fft_search(rx, training, testCase.Rs, testCase.K);

            err_kHz = abs(freq_est_kHz - trueDeltaF_kHz);
            fprintf('  [fft_search] Nominal-SNR estimate: %.3f kHz  (true %.3f kHz, error %.3f kHz)\n', ...
                    freq_est_kHz, trueDeltaF_kHz, err_kHz);
        end

        function testFftSearch_NegativeOffset(testCase)
            %TESTFFTSEARCH_NEGATIVEOFFSET  Informational: negative offset accuracy.
            DeltaF_neg     = -200;
            trueDeltaF_kHz = DeltaF_neg * 1e3;

            [rx, training, ~] = test_FreqRecovery.buildRx( ...
                DeltaF_neg, 40, testCase.Rs);

            [~, freq_est_kHz] = freq_recovery.fft_search(rx, training, testCase.Rs, testCase.K);

            err_kHz = abs(freq_est_kHz - trueDeltaF_kHz);
            fprintf('  [fft_search] Negative-offset estimate: %.3f kHz  (true %.3f kHz, error %.3f kHz)\n', ...
                    freq_est_kHz, trueDeltaF_kHz, err_kHz);
        end

        function testFftSearch_BERAfterCorrection(testCase)
            %TESTFFTSEARCH_BERAFTERCORRECTION
            %   After correction, BER must be below 1e-2 at nominal SNR.
            [rx, training, symbols] = test_FreqRecovery.buildRx( ...
                testCase.DeltaF, testCase.SNR_dB, testCase.Rs);

            [y, freq_est_kHz] = freq_recovery.fft_search(rx, training, testCase.Rs, testCase.K);

            txRefBits = modem.symbolsToBits(symbols);
            best      = test_FreqRecovery.bestRotation(y, txRefBits);

            rxBits = modem.symbolsToBits(modem.decideSymbols(best));
            BER    = sum(rxBits ~= txRefBits) / numel(txRefBits);

            fprintf('  [fft_search] BER after correction: %.2e  (est offset %.3f kHz, true %.3f kHz)\n', ...
                    BER, freq_est_kHz, testCase.DeltaF * 1e3);

            % --- Constellation plots ---------------------------------
            ms = 2;
            figure('Name', sprintf('fft_search | \\DeltaF = %d MHz | SNR = %d dB', ...
                   testCase.DeltaF, testCase.SNR_dB), ...
                   'Position', [150 150 900 400], 'Color', 'w');

            subplot(1, 3, 1);
            plot(real(symbols(:,1)), imag(symbols(:,1)), '.', 'MarkerSize', ms);
            grid on; axis equal;
            title('TX symbols (X-pol)');
            xlabel('I'); ylabel('Q');

            subplot(1, 3, 2);
            plot(real(rx(:,1)), imag(rx(:,1)), '.', 'MarkerSize', ms);
            grid on; axis equal;
            title('After LO shift + AWGN');
            xlabel('I'); ylabel('Q');

            subplot(1, 3, 3);
            plot(real(best(:,1)), imag(best(:,1)), '.', 'MarkerSize', ms);
            grid on; axis equal;
            title(sprintf('After correction (BER = %.2e)', BER));
            xlabel('I'); ylabel('Q');

            sgtitle(sprintf('fft\\_search  |  \\DeltaF = %d MHz  |  SNR = %d dB  |  est = %.1f kHz', ...
                    testCase.DeltaF, testCase.SNR_dB, freq_est_kHz));

            % --- Only fail on BER -----------------------------------
            testCase.verifyLessThan(BER, 1e-2, ...
                sprintf('[fft_search] BER %.2e after frequency correction exceeds 1e-2.', BER));
        end

        % ============================================================
        %  fitz tests
        % ============================================================

        function testFitz_OutputShape(testCase)
            %TESTFITZ_OUTPUTSHAPE  y must be the same size as the input.
            [rx, training, ~] = test_FreqRecovery.buildRx( ...
                testCase.DeltaF, testCase.SNR_dB, testCase.Rs);

            [y, ~] = freq_recovery.fitz(rx, training, testCase.Rs);

            testCase.verifySize(y, size(rx), ...
                'fitz output must be the same size as the input.');
        end

        function testFitz_HighSNR(testCase)
            %TESTFITZ_HIGHSNR  Informational: accuracy at 40 dB SNR.
            trueDeltaF_kHz = testCase.DeltaF * 1e3;

            [rx, training, ~] = test_FreqRecovery.buildRx( ...
                testCase.DeltaF, 40, testCase.Rs);

            [~, freq_est_kHz] = freq_recovery.fitz(rx, training, testCase.Rs);

            err_kHz = abs(freq_est_kHz - trueDeltaF_kHz);
            fprintf('  [fitz] High-SNR estimate: %.3f kHz  (true %.3f kHz, error %.3f kHz)\n', ...
                    freq_est_kHz, trueDeltaF_kHz, err_kHz);
        end

        function testFitz_NominalSNR(testCase)
            %TESTFITZ_NOMINALSNR  Informational: accuracy at nominal SNR.
            trueDeltaF_kHz = testCase.DeltaF * 1e3;

            [rx, training, ~] = test_FreqRecovery.buildRx( ...
                testCase.DeltaF, testCase.SNR_dB, testCase.Rs);

            [~, freq_est_kHz] = freq_recovery.fitz(rx, training, testCase.Rs);

            err_kHz = abs(freq_est_kHz - trueDeltaF_kHz);
            fprintf('  [fitz] Nominal-SNR estimate: %.3f kHz  (true %.3f kHz, error %.3f kHz)\n', ...
                    freq_est_kHz, trueDeltaF_kHz, err_kHz);
        end

        function testFitz_NegativeOffset(testCase)
            %TESTFITZ_NEGATIVEOFFSET  Informational: negative offset accuracy.
            DeltaF_neg     = -200;
            trueDeltaF_kHz = DeltaF_neg * 1e3;

            [rx, training, ~] = test_FreqRecovery.buildRx( ...
                DeltaF_neg, 40, testCase.Rs);

            [~, freq_est_kHz] = freq_recovery.fitz(rx, training, testCase.Rs);

            err_kHz = abs(freq_est_kHz - trueDeltaF_kHz);
            fprintf('  [fitz] Negative-offset estimate: %.3f kHz  (true %.3f kHz, error %.3f kHz)\n', ...
                    freq_est_kHz, trueDeltaF_kHz, err_kHz);
        end

        function testFitz_BERAfterCorrection(testCase)
            %TESTFITZ_BERAFTERCORRECTION
            %   After correction, BER must be below 1e-2 at nominal SNR.
            [rx, training, symbols] = test_FreqRecovery.buildRx( ...
                testCase.DeltaF, testCase.SNR_dB, testCase.Rs);

            [y, freq_est_kHz] = freq_recovery.fitz(rx, training, testCase.Rs);

            txRefBits = modem.symbolsToBits(symbols);
            best      = test_FreqRecovery.bestRotation(y, txRefBits);

            rxBits = modem.symbolsToBits(modem.decideSymbols(best));
            BER    = sum(rxBits ~= txRefBits) / numel(txRefBits);

            fprintf('  [fitz] BER after correction: %.2e  (est offset %.3f kHz, true %.3f kHz)\n', ...
                    BER, freq_est_kHz, testCase.DeltaF * 1e3);

            % --- Constellation plots ---------------------------------
            ms = 2;
            figure('Name', sprintf('fitz | \\DeltaF = %d MHz | SNR = %d dB', ...
                   testCase.DeltaF, testCase.SNR_dB), ...
                   'Position', [200 200 900 400], 'Color', 'w');

            subplot(1, 3, 1);
            plot(real(symbols(:,1)), imag(symbols(:,1)), '.', 'MarkerSize', ms);
            grid on; axis equal;
            title('TX symbols (X-pol)');
            xlabel('I'); ylabel('Q');

            subplot(1, 3, 2);
            plot(real(rx(:,1)), imag(rx(:,1)), '.', 'MarkerSize', ms);
            grid on; axis equal;
            title('After LO shift + AWGN');
            xlabel('I'); ylabel('Q');

            subplot(1, 3, 3);
            plot(real(best(:,1)), imag(best(:,1)), '.', 'MarkerSize', ms);
            grid on; axis equal;
            title(sprintf('After correction (BER = %.2e)', BER));
            xlabel('I'); ylabel('Q');

            sgtitle(sprintf('fitz  |  \\DeltaF = %d MHz  |  SNR = %d dB  |  est = %.1f kHz', ...
                    testCase.DeltaF, testCase.SNR_dB, freq_est_kHz));

            % --- Only fail on BER -----------------------------------
            testCase.verifyLessThan(BER, 1e-2, ...
                sprintf('[fitz] BER %.2e after frequency correction exceeds 1e-2.', BER));
        end

    end

    % ================================================================
    %  Private static helpers
    % ================================================================
    methods (Static, Access = private)

        function [rx, training, symbols] = buildRx(DeltaF_MHz, SNR_dB, Rs)
            %BUILDRX  Generate one subframe, apply LO shift and AWGN.
            %
            %   DeltaF_MHz - frequency offset [MHz]
            %   SNR_dB     - SNR [dB]  (Inf for noise-free)
            %   Rs         - symbol rate [GBd]

            % Generate exactly one subframe
            DATA_PER_SUBFRAME = 3586;
            Nbits = DATA_PER_SUBFRAME * 2 * 2;   % 2 pol, 2 bits per sym per pol

            bits = modem.randomBits(Nbits);
            [symbols, ~, training, ~] = modem.modulate(bits);

            % Take the first subframe only (training is [11 x 2], same each SF)
            SUBFRAME_SYMS = 3712;
            symbols = symbols(1:SUBFRAME_SYMS, :);

            % Channel: LO shift then AWGN (symbol-rate signal, SpS = 1)
            rx = channel.lo_freq_shift(symbols, DeltaF_MHz, Rs, 1);
            if isfinite(SNR_dB)
                rx = channel.add_awgn(rx, SNR_dB);
            end
        end

        function best = bestRotation(y, refBits)
            %BESTROTATION  Return y rotated by the k*pi/2 that minimises BER.
            rotations = [1, 1j, -1, -1j];
            bestBER   = Inf;
            best      = y;

            for ri = 1:4
                rotated  = y * rotations(ri);
                rxBits   = modem.symbolsToBits(modem.decideSymbols(rotated));
                thisBER  = sum(rxBits ~= refBits) / numel(refBits);
                if thisBER < bestBER
                    bestBER = thisBER;
                    best    = rotated;
                end
            end
        end

    end
end
