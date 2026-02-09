classdef test_QAMModem < matlab.unittest.TestCase

    properties (TestParameter)
        M = {4, 16, 64, 256}
        N_pol = {1, 2}
    end

    methods (Test)
        function testModulateDemodulateRoundTrip(testCase, M, N_pol)
            % Verify that bits survive a modulate -> symbolsToBits round trip
            modem = QAMModem(M, N_pol);
            k = log2(M);
            Nbits = k * N_pol * 128;           % 128 symbols per pol

            bits = modem.randomBits(Nbits);
            symbols = modem.modulate(bits);
            recoveredBits = modem.symbolsToBits(symbols);

            testCase.verifyEqual(recoveredBits, bits, ...
                sprintf('Round-trip failed for %d-QAM, %d pol(s).', M, N_pol));
        end

        function testDecideSymbolsRecoversBits(testCase, M, N_pol)
            % Verify that decideSymbols on clean symbols gives exact match
            modem = QAMModem(M, N_pol);
            k = log2(M);
            Nbits = k * N_pol * 64;

            bits = modem.randomBits(Nbits);
            symbols = modem.modulate(bits);
            decided = modem.decideSymbols(symbols);

            testCase.verifyEqual(decided, symbols, 'AbsTol', 1e-10, ...
                'decideSymbols should return the same symbols for clean input.');

            recoveredBits = modem.symbolsToBits(decided);
            testCase.verifyEqual(recoveredBits, bits, ...
                'Bits should survive modulate -> decide -> symbolsToBits.');
        end

        function testSymbolDimensions(testCase, M, N_pol)
            % Verify output symbol array has the expected shape
            modem = QAMModem(M, N_pol);
            k = log2(M);
            Ns = 50;
            Nbits = k * N_pol * Ns;

            bits = modem.randomBits(Nbits);
            symbols = modem.modulate(bits);

            testCase.verifySize(symbols, [Ns, N_pol]);
        end

        function testUnitAveragePower(testCase, M, N_pol)
            % Verify constellation has approximately unit average power
            modem = QAMModem(M, N_pol);
            k = log2(M);
            Nbits = k * N_pol * 4096;

            bits = modem.randomBits(Nbits);
            symbols = modem.modulate(bits);

            avgPower = mean(abs(symbols(:)).^2);
            testCase.verifyEqual(avgPower, 1, 'AbsTol', 0.05, ...
                'Average symbol power should be ~1.');
        end
    end

    methods (Test)
        function testInvalidMThrows(testCase)
            threw = false;
            try
                QAMModem(3, 1);
            catch
                threw = true;
            end
            testCase.verifyTrue(threw, ...
                'QAMModem with non-power-of-2 M should throw an error.');
        end

        function testNotEnoughBitsThrows(testCase)
            modem = QAMModem(16, 2);
            threw = false;
            try
                modem.modulate([0;1]);
            catch
                threw = true;
            end
            testCase.verifyTrue(threw, ...
                'modulate with too few bits should throw an error.');
        end

        function testRandomBitsLength(testCase)
            modem = QAMModem(16, 1);
            bits = modem.randomBits(200);
            testCase.verifyLength(bits, 200);
        end

        function testRandomBitsBinary(testCase)
            modem = QAMModem(16, 1);
            bits = modem.randomBits(500);
            testCase.verifyTrue(all(bits == 0 | bits == 1), ...
                'randomBits should only produce 0s and 1s.');
        end
    end
end
