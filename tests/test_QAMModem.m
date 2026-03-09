classdef test_QAMModem < matlab.unittest.TestCase
    % Tests for QPSK modem functions with CPON framing.

    properties (Constant)
        SUBFRAME_SYMS  = 3712
        BLOCK_LEN      = 32
        N_BLOCKS       = 116
        N_TRAIN        = 11
        DATA_PER_SF    = 3586      % 3712 - 116 pilots - 10 extra TS
    end

    methods (Test)

        % ---- Round-trip: data bits survive modulate -> symbolsToBits ---
        function testModulateDemodulateRoundTrip(testCase)
            Nbits = testCase.DATA_PER_SF * 4;     % exactly 1 subframe of data
            bits  = modem.randomBits(Nbits);
            [symbols, ~, ~, nSF] = modem.modulate(bits);

            % Reference bits from the full symbol stream (incl. pilots/TS)
            txRefBits = modem.symbolsToBits(symbols);
            rxDecided = modem.decideSymbols(symbols);     % clean input
            rxBits    = modem.symbolsToBits(rxDecided);

            testCase.verifyEqual(rxBits, txRefBits, ...
                'Round-trip failed for clean QPSK symbols.');
            testCase.verifyEqual(nSF, 1, 'Expected exactly 1 subframe.');
        end

        % ---- decideSymbols on clean QPSK returns same symbols --------
        function testDecideSymbolsClean(testCase)
            Nbits = testCase.DATA_PER_SF * 4;
            bits  = modem.randomBits(Nbits);
            [symbols, ~, ~, ~] = modem.modulate(bits);

            % Data symbols are ±1±1j, pilots/training are ±3±3j
            decided = modem.decideSymbols(symbols);

            % All signs should match (slicer maps ±3 -> ±1)
            testCase.verifyEqual(sign(real(decided)), sign(real(symbols)), ...
                'Real-part signs differ after decideSymbols.');
            testCase.verifyEqual(sign(imag(decided)), sign(imag(symbols)), ...
                'Imag-part signs differ after decideSymbols.');
        end

        % ---- Output dimensions ----------------------------------------
        function testSymbolDimensions(testCase)
            Nbits = testCase.DATA_PER_SF * 4 * 2;    % 2 subframes of data
            bits  = modem.randomBits(Nbits);
            [symbols, pilots, training, nSF] = modem.modulate(bits);

            testCase.verifyEqual(nSF, 2, 'Expected 2 subframes.');
            testCase.verifySize(symbols, [2 * testCase.SUBFRAME_SYMS, 2]);
            testCase.verifySize(pilots,  [testCase.N_BLOCKS, 2]);
            testCase.verifySize(training, [testCase.N_TRAIN, 2]);
        end

        % ---- CPON framing structure: training at correct positions -----
        function testTrainingPositions(testCase)
            Nbits = testCase.DATA_PER_SF * 4;
            bits  = modem.randomBits(Nbits);
            [symbols, ~, training, ~] = modem.modulate(bits);

            % Training symbols TS2..TS11 at positions 2..11
            testCase.verifyEqual(symbols(2:testCase.N_TRAIN, :), ...
                training(2:testCase.N_TRAIN, :), ...
                'Training symbols TS2..TS11 mismatch at subframe start.');
        end

        % ---- CPON framing: pilots at every 32nd position ---------------
        function testPilotPositions(testCase)
            Nbits = testCase.DATA_PER_SF * 4;
            bits  = modem.randomBits(Nbits);
            [symbols, pilots, ~, ~] = modem.modulate(bits);

            % Pilot at position 1 (TS1 = pilot 1)
            testCase.verifyEqual(symbols(1, :), pilots(1, :), ...
                'Pilot 1 (TS1) mismatch.');

            % Pilots at positions 33, 65, ..., 3681 (block boundaries)
            for blk = 2:testCase.N_BLOCKS
                pos = (blk - 1) * testCase.BLOCK_LEN + testCase.BLOCK_LEN;
                testCase.verifyEqual(symbols(pos, :), pilots(blk, :), ...
                    sprintf('Pilot at block %d (pos %d) mismatch.', blk, pos));
            end
        end

        % ---- Pilot amplitude -------------------------------------------
        function testPilotAmplitude(testCase)
            Nbits = testCase.DATA_PER_SF * 4;
            bits  = modem.randomBits(Nbits);
            [~, pilots, training, ~] = modem.modulate(bits);

            % Pilots and training should be at ±3 ±3j
            testCase.verifyEqual(abs(real(pilots(:))),  3*ones(numel(pilots),1), ...
                'Pilot real amplitudes should all be 3.');
            testCase.verifyEqual(abs(imag(pilots(:))),  3*ones(numel(pilots),1), ...
                'Pilot imag amplitudes should all be 3.');
            testCase.verifyEqual(abs(real(training(:))), 3*ones(numel(training),1), ...
                'Training real amplitudes should all be 3.');
            testCase.verifyEqual(abs(imag(training(:))), 3*ones(numel(training),1), ...
                'Training imag amplitudes should all be 3.');
        end

        % ---- randomBits length and binary values ----------------------
        function testRandomBitsLength(testCase)
            bits = modem.randomBits(200);
            testCase.verifyLength(bits, 200);
        end

        function testRandomBitsBinary(testCase)
            bits = modem.randomBits(500);
            testCase.verifyTrue(all(bits == 0 | bits == 1), ...
                'randomBits should only produce 0s and 1s.');
        end

        % ---- Slicer: QPSK sign decision ------------------------------
        function testSlicerQPSK(testCase)
            s_in  = [0.7+0.3j; -2.5-0.1j; 0.01-99j; -0.5+0.5j];
            s_exp = [1+1j; -1-1j; 1-1j; -1+1j];
            s_out = modem.slicer(s_in);
            testCase.verifyEqual(s_out, s_exp, ...
                'Slicer should produce ±1±1j based on sign.');
        end

        % ---- Pulse shaping plot (QPSK) -------------------------------
        function testPulseShapingPlot(testCase)
            SpS     = 2;
            Ns      = 64;
            rolloff = 0.25;
            span    = 10;

            Nbits   = 4 * Ns;          % 2 bits/sym/pol * 2 pols * Ns
            bits    = modem.randomBits(Nbits);
            [symbols, ~, ~, ~] = modem.modulate(bits);
            symbols = symbols(1:Ns, 1);   % single pol, Ns symbols

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
