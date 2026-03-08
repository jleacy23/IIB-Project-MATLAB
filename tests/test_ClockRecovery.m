classdef test_ClockRecovery < matlab.unittest.TestCase
    %TEST_CLOCKRECOVERY  Verify clk_recovery corrects constant timing
    %  offsets and sampling-frequency offsets on RRC-shaped signals.
    %
    %  Approach: symbols are pulse-shaped at a high oversampling factor
    %  (SpS_hi = 16) so that arbitrary sub-sample timing shifts can be
    %  applied by selecting different sample phases.  The signal is then
    %  decimated to 2 Sa/symbol before being fed to clk_recovery.

    properties (Constant)
        M       = 4              % QPSK
        N_pol   = 1              % single polarisation (clk_recovery works per-pol)
        Ns      = 4096           % symbols
        SpS     = 2              % target samples per symbol
        SpS_hi  = 16             % high-resolution oversampling

        % RRC parameters
        Rolloff = 0.25
        Span    = 10

        % DPLL loop-filter constants (empirically tuned for Nyquist TED)
        ki = 1e-4
        kp = 5e-3

        % Pass/fail
        BER_THRESHOLD = 1e-2
    end

    methods (TestMethodSetup)
        function seedRng(~)
            rng(42);
        end
    end

    % ================================================================
    methods (Test)

        % -------- Constant timing offset ----------------------------
        function testConstantTimingOffset(testCase)
            SpS_hi_ = testCase.SpS_hi;
            SpS_    = testCase.SpS;
            Ns_     = testCase.Ns;
            M_      = testCase.M;

            % --- Tx: generate symbols & pulse-shape at high SpS ---
            k      = log2(M_);
            Nbits  = k * testCase.N_pol * Ns_;
            txBits = modem.randomBits(Nbits);
            symbols = modem.modulate(txBits, M_, testCase.N_pol);

            txHi = modem.rrcPulse(symbols, SpS_hi_, testCase.Rolloff, testCase.Span);

            % --- Matched filter at high SpS ---
            rxHi = modem.matched_filter(txHi, SpS_hi_, 'rrc', ...
                testCase.Rolloff, testCase.Span);

            % --- Introduce a constant timing offset ---
            %  Shift by 3 high-res samples = 3/16 of a symbol period
            timingOffset = 3;   % samples at SpS_hi
            rxShifted = circshift(rxHi, timingOffset);

            % --- Decimate to 2 Sa/symbol ---
            decFactor = SpS_hi_ / SpS_;
            rx2 = rxShifted(1:decFactor:end, :);

            % --- Clock recovery ---
            crOut = clk_recovery.recovery(rx2, 'Nyquist', Ns_, ...
                testCase.ki, testCase.kp);

            % --- Downsample to symbol rate & demodulate ---
            crSym = crOut(1:SpS_:end);

            % Resolve phase ambiguity
            [BER, crSym] = bestRotationBER(testCase, crSym, txBits, M_);

            fprintf('Constant offset BER = %.2e\n', BER);

            % --- Plot ---
            rxSym = rx2(1:SpS_:end);
            plotBeforeAfter(testCase, rxSym, crSym, ...
                'Constant Timing Offset', M_, BER);

            % --- Verify ---
            testCase.verifyTrue(all(isfinite(crSym(:))), ...
                'Clock-recovery output contains NaN/Inf.');
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD, ...
                sprintf('Constant-offset BER %.2e exceeds threshold.', BER));
        end

        % -------- Sampling-frequency offset (drifting timing) -------
        function testSamplingFrequencyOffset(testCase)
            SpS_hi_ = testCase.SpS_hi;
            SpS_    = testCase.SpS;
            Ns_     = testCase.Ns;
            M_      = testCase.M;

            % --- Tx ---
            k      = log2(M_);
            Nbits  = k * testCase.N_pol * Ns_;
            txBits = modem.randomBits(Nbits);
            symbols = modem.modulate(txBits, M_, testCase.N_pol);

            txHi = modem.rrcPulse(symbols, SpS_hi_, testCase.Rolloff, testCase.Span);

            % --- Matched filter at high SpS ---
            rxHi = modem.matched_filter(txHi, SpS_hi_, 'rrc', ...
                testCase.Rolloff, testCase.Span);

            % --- Introduce a sampling-frequency offset ---
            %  Resample as if the receiver clock runs at a slightly
            %  different rate.  A 100 ppm offset is typical.
            ppm = 100;                           % parts per million
            Nhi = size(rxHi, 1);
            tOrig    = (0:Nhi-1).';              % original sample instants
            tSkewed  = tOrig * (1 + ppm*1e-6);   % skewed sample instants
            rxSkewed = interp1(tOrig, rxHi, tSkewed, 'spline', 0);

            % --- Decimate to 2 Sa/symbol ---
            decFactor = SpS_hi_ / SpS_;
            rx2 = rxSkewed(1:decFactor:end, :);

            % --- Clock recovery ---
            crOut = clk_recovery.recovery(rx2, 'Nyquist', Ns_, ...
                testCase.ki, testCase.kp);

            % --- Downsample to symbol rate & demodulate ---
            crSym = crOut(1:SpS_:end);

            % Resolve phase ambiguity
            [BER, crSym] = bestRotationBER(testCase, crSym, txBits, M_);

            fprintf('SFO (%d ppm) BER = %.2e\n', ppm, BER);

            % --- Plot ---
            rxSym = rx2(1:SpS_:end);
            plotBeforeAfter(testCase, rxSym, crSym, ...
                sprintf('SFO %d ppm', ppm), M_, BER);

            % --- Verify ---
            testCase.verifyTrue(all(isfinite(crSym(:))), ...
                'Clock-recovery output contains NaN/Inf.');
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD, ...
                sprintf('SFO BER %.2e exceeds threshold.', BER));
        end

    end

    % ================================================================
    %  Helpers
    % ================================================================
    methods (Access = private)

        function [BER, bestSym] = bestRotationBER(testCase, crSym, txBits, M)
            %BESTROTATIONBER  Try all four pi/2 rotations, return lowest BER.
            bestBER = Inf;
            bestSym = crSym;
            N_pol_ = testCase.N_pol;
            for kk = 0:3
                rotated = crSym .* exp(-1j * kk * pi/2);
                decided = modem.decideSymbols(rotated, M, N_pol_);
                rxBits  = modem.symbolsToBits(decided, M);
                nBits   = min(length(rxBits), length(txBits));
                nErr    = sum(txBits(1:nBits) ~= rxBits(1:nBits));
                thisBER = nErr / nBits;
                if thisBER < bestBER
                    bestBER = thisBER;
                    bestSym = rotated;
                end
            end
            BER = bestBER;
        end

        function plotBeforeAfter(~, rxSym, crSym, titleStr, M, BER)
            figure('Name', titleStr, 'Position', [100 100 900 400]);

            subplot(1, 2, 1);
            plot(real(rxSym), imag(rxSym), '.', 'MarkerSize', 2);
            grid on; axis equal;
            title('Before Clock Recovery');
            xlabel('I'); ylabel('Q');

            subplot(1, 2, 2);
            plot(real(crSym), imag(crSym), '.', 'MarkerSize', 2);
            grid on; axis equal;
            title('After Clock Recovery');
            xlabel('I'); ylabel('Q');

            sgtitle(sprintf('%d-QAM  |  %s  |  BER = %.2e', M, titleStr, BER));
        end
    end
end
