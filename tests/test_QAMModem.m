classdef test_QAMModem < matlab.unittest.TestCase

    properties (TestParameter)
        M = {4, 16, 64, 256}
        N_pol = {1, 2}
    end

    methods (Test)
        function testModulateDemodulateRoundTrip(testCase, M, N_pol)
            % Verify that bits survive a modulate -> symbolsToBits round trip
            k = log2(M);
            Nbits = k * N_pol * 128;           % 128 symbols per pol

            bits = modem.randomBits(Nbits);
            symbols = modem.modulate(bits, M, N_pol);
            recoveredBits = modem.symbolsToBits(symbols, M);

            testCase.verifyEqual(recoveredBits, bits, ...
                sprintf('Round-trip failed for %d-QAM, %d pol(s).', M, N_pol));
        end

        function testDecideSymbolsRecoversBits(testCase, M, N_pol)
            % Verify that decideSymbols on clean symbols gives exact match
            k = log2(M);
            Nbits = k * N_pol * 64;

            bits = modem.randomBits(Nbits);
            symbols = modem.modulate(bits, M, N_pol);
            decided = modem.decideSymbols(symbols, M, N_pol);

            testCase.verifyEqual(decided, symbols, 'AbsTol', 1e-10, ...
                'decideSymbols should return the same symbols for clean input.');

            recoveredBits = modem.symbolsToBits(decided, M);
            testCase.verifyEqual(recoveredBits, bits, ...
                'Bits should survive modulate -> decide -> symbolsToBits.');
        end

        function testSymbolDimensions(testCase, M, N_pol)
            % Verify output symbol array has the expected shape
            k = log2(M);
            Ns = 50;
            Nbits = k * N_pol * Ns;

            bits = modem.randomBits(Nbits);
            symbols = modem.modulate(bits, M, N_pol);

            testCase.verifySize(symbols, [Ns, N_pol]);
        end

        function testUnitAveragePower(testCase, M, N_pol)
            % Verify constellation has approximately unit average power
            k = log2(M);
            Nbits = k * N_pol * 4096;

            bits = modem.randomBits(Nbits);
            symbols = modem.modulate(bits, M, N_pol);

            avgPower = mean(abs(symbols(:)).^2);
            testCase.verifyEqual(avgPower, 1, 'AbsTol', 0.05, ...
                'Average symbol power should be ~1.');
        end
    end

    methods (Test)
        function testInvalidMThrows(testCase)
            threw = false;
            try
                modem.modulate([0;1;0;1], 3, 1);
            catch
                threw = true;
            end
            testCase.verifyTrue(threw, ...
                'qam_modulate with non-power-of-2 M should throw an error.');
        end

        function testNotEnoughBitsThrows(testCase)
            threw = false;
            try
                modem.modulate([0;1], 16, 2);
            catch
                threw = true;
            end
            testCase.verifyTrue(threw, ...
                'modulate with too few bits should throw an error.');
        end

        function testRandomBitsLength(testCase)
            bits = modem.randomBits(200);
            testCase.verifyLength(bits, 200);
        end

        function testRandomBitsBinary(testCase)
            bits = modem.randomBits(500);
            testCase.verifyTrue(all(bits == 0 | bits == 1), ...
                'randomBits should only produce 0s and 1s.');
        end

        function testPulseShapingPlot(testCase)
            % Plot matched-filter output for rect and RRC pulse shaping
            M_  = 4;
            SpS = 2;
            Ns  = 64;
            k   = log2(M_);
            rolloff = 0.25;
            span    = 10;

            bits    = modem.randomBits(k * Ns);
            symbols = modem.modulate(bits, M_, 1);

            % Rectangular
            txRect = modem.rectPulse(symbols, SpS);
            rxRect = modem.matched_filter(txRect, SpS, 'rect');

            % RRC
            txRRC = modem.rrcPulse(symbols, SpS, rolloff, span);
            rxRRC = modem.matched_filter(txRRC, SpS, 'rrc', rolloff, span);

            t = (0:size(txRect,1)-1).' / SpS;

            figure('Name', 'Matched Filter Output', ...
                   'Position', [100 100 1200 800]);

            subplot(2,2,1);
            plot(t, real(txRect), t, real(rxRect), 'LineWidth', 0.8);
            hold on;
            stem((0:Ns-1).', real(symbols), 'k', 'MarkerSize', 3);
            hold off;
            xlabel('Symbol index'); ylabel('Re');
            title('Rectangular – real part');
            legend('Tx (shaped)', 'Rx (matched)', 'Symbols');
            xlim([0 20]); grid on;

            subplot(2,2,2);
            plot(t, imag(txRect), t, imag(rxRect), 'LineWidth', 0.8);
            hold on;
            stem((0:Ns-1).', imag(symbols), 'k', 'MarkerSize', 3);
            hold off;
            xlabel('Symbol index'); ylabel('Im');
            title('Rectangular – imag part');
            legend('Tx (shaped)', 'Rx (matched)', 'Symbols');
            xlim([0 20]); grid on;

            subplot(2,2,3);
            plot(t, real(txRRC), t, real(rxRRC), 'LineWidth', 0.8);
            hold on;
            stem((0:Ns-1).', real(symbols), 'k', 'MarkerSize', 3);
            hold off;
            xlabel('Symbol index'); ylabel('Re');
            title(sprintf('RRC (\\beta=%.2f) – real part', rolloff));
            legend('Tx (shaped)', 'Rx (matched)', 'Symbols');
            xlim([0 20]); grid on;

            subplot(2,2,4);
            plot(t, imag(txRRC), t, imag(rxRRC), 'LineWidth', 0.8);
            hold on;
            stem((0:Ns-1).', imag(symbols), 'k', 'MarkerSize', 3);
            hold off;
            xlabel('Symbol index'); ylabel('Im');
            title(sprintf('RRC (\\beta=%.2f) – imag part', rolloff));
            legend('Tx (shaped)', 'Rx (matched)', 'Symbols');
            xlim([0 20]); grid on;

            sgtitle('Pulse Shaping & Matched Filter Output (QPSK, 1 pol)');

            testCase.verifyTrue(true);
        end
    end
end
