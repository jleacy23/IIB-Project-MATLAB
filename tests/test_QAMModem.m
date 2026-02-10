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

            bits = qam_randomBits(Nbits);
            symbols = qam_modulate(bits, M, N_pol);
            recoveredBits = qam_symbolsToBits(symbols, M);

            testCase.verifyEqual(recoveredBits, bits, ...
                sprintf('Round-trip failed for %d-QAM, %d pol(s).', M, N_pol));
        end

        function testDecideSymbolsRecoversBits(testCase, M, N_pol)
            % Verify that decideSymbols on clean symbols gives exact match
            k = log2(M);
            Nbits = k * N_pol * 64;

            bits = qam_randomBits(Nbits);
            symbols = qam_modulate(bits, M, N_pol);
            decided = qam_decideSymbols(symbols, M, N_pol);

            testCase.verifyEqual(decided, symbols, 'AbsTol', 1e-10, ...
                'decideSymbols should return the same symbols for clean input.');

            recoveredBits = qam_symbolsToBits(decided, M);
            testCase.verifyEqual(recoveredBits, bits, ...
                'Bits should survive modulate -> decide -> symbolsToBits.');
        end

        function testSymbolDimensions(testCase, M, N_pol)
            % Verify output symbol array has the expected shape
            k = log2(M);
            Ns = 50;
            Nbits = k * N_pol * Ns;

            bits = qam_randomBits(Nbits);
            symbols = qam_modulate(bits, M, N_pol);

            testCase.verifySize(symbols, [Ns, N_pol]);
        end

        function testUnitAveragePower(testCase, M, N_pol)
            % Verify constellation has approximately unit average power
            k = log2(M);
            Nbits = k * N_pol * 4096;

            bits = qam_randomBits(Nbits);
            symbols = qam_modulate(bits, M, N_pol);

            avgPower = mean(abs(symbols(:)).^2);
            testCase.verifyEqual(avgPower, 1, 'AbsTol', 0.05, ...
                'Average symbol power should be ~1.');
        end
    end

    methods (Test)
        function testInvalidMThrows(testCase)
            threw = false;
            try
                qam_modulate([0;1;0;1], 3, 1);
            catch
                threw = true;
            end
            testCase.verifyTrue(threw, ...
                'qam_modulate with non-power-of-2 M should throw an error.');
        end

        function testNotEnoughBitsThrows(testCase)
            threw = false;
            try
                qam_modulate([0;1], 16, 2);
            catch
                threw = true;
            end
            testCase.verifyTrue(threw, ...
                'modulate with too few bits should throw an error.');
        end

        function testRandomBitsLength(testCase)
            bits = qam_randomBits(200);
            testCase.verifyLength(bits, 200);
        end

        function testRandomBitsBinary(testCase)
            bits = qam_randomBits(500);
            testCase.verifyTrue(all(bits == 0 | bits == 1), ...
                'randomBits should only produce 0s and 1s.');
        end
    end
end
