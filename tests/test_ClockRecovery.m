classdef test_ClockRecovery < matlab.unittest.TestCase
    %TEST_CLOCKRECOVERY  Verify clk_recovery.recovery_godard corrects
    %  constant timing offsets and sampling-frequency offsets on
    %  RRC-shaped signals, in both feedforward and feedback modes.
    %
    %  Approach: symbols are pulse-shaped at a high oversampling factor
    %  (SpS_hi = 16) so that arbitrary sub-sample timing shifts can be
    %  applied by selecting different sample phases.  The signal is then
    %  decimated to 2 Sa/symbol before being fed to recovery_godard.

    properties (Constant)
        N_pol   = 1              % single polarisation (per-pol operation)
        Ns      = 2^18           % symbols
        SpS     = 2              % target samples per symbol
        SpS_hi  = 16             % high-resolution oversampling

        % RRC parameters
        Rolloff = 0.25
        Span    = 10

        % Modified Godard estimator parameters
        N_fft = 256

        % PI loop-filter gains for feedback-mode Godard (tuned for the
        % unnormalised imag(S) error signal used by recovery_godard).
        ki_godard = 1e-5
        kp_godard = 1e-4

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

        % -------- Modified Godard feedforward: constant timing -------
        function testGodardFF_ConstantTimingOffset(testCase)
            runGodardScenario(testCase, 'feedforward', 'constant');
        end

        % -------- Modified Godard feedforward: SFO -------------------
        function testGodardFF_SamplingFrequencyOffset(testCase)
            runGodardScenario(testCase, 'feedforward', 'sfo');
        end

        % -------- Modified Godard feedback: constant timing ----------
        function testGodardFB_ConstantTimingOffset(testCase)
            runGodardScenario(testCase, 'feedback', 'constant');
        end

        % -------- Modified Godard feedback: SFO ----------------------
        function testGodardFB_SamplingFrequencyOffset(testCase)
            runGodardScenario(testCase, 'feedback', 'sfo');
        end

    end

    % ================================================================
    %  Helpers
    % ================================================================
    methods (Access = private)

        function runGodardScenario(testCase, mode, impairment)
            SpS_hi_ = testCase.SpS_hi;
            SpS_    = testCase.SpS;
            Ns_     = testCase.Ns;

            % --- Tx: generate symbols & pulse-shape at high SpS ---
            Nbits  = 2 * Ns_;              % QPSK: 2 bits/sym, single pol
            txBits = modem.randomBits(Nbits);
            symbols = modem.modulate(txBits);
            symbols = symbols(:, 1);        % single pol

            txHi = modem.rrcPulse(symbols, SpS_hi_, testCase.Rolloff, testCase.Span);

            % --- Matched filter at high SpS ---
            rxHi = modem.matched_filter(txHi, SpS_hi_, 'rrc', ...
                testCase.Rolloff, testCase.Span);

            % --- Impairment (via channel.apply_timing_error) ---
            switch impairment
                case 'constant'
                    tau0  = 3 / SpS_hi_;    % 3/16 of a symbol period
                    rxImp = channel.apply_timing_error(rxHi, 0, tau0, SpS_hi_);
                    impLabel = sprintf('Constant Offset tau0 = %.3f T', tau0);
                case 'sfo'
                    ppm   = 100;
                    rxImp = channel.apply_timing_error(rxHi, ppm, 0, SpS_hi_);
                    impLabel = sprintf('SFO %d ppm', ppm);
                otherwise
                    error('Unknown impairment: %s', impairment);
            end

            % --- Decimate to 2 Sa/symbol ---
            decFactor = SpS_hi_ / SpS_;
            rx2 = rxImp(1:decFactor:end, :);

            % --- Clock recovery (Modified Godard) ---
            switch mode
                case 'feedforward'
                    crOut = clk_recovery.recovery_godard(rx2, Ns_, ...
                        testCase.N_fft, testCase.Rolloff, 'feedforward');
                    modeLabel = 'FF';
                case 'feedback'
                    crOut = clk_recovery.recovery_godard(rx2, Ns_, ...
                        testCase.N_fft, testCase.Rolloff, 'feedback', ...
                        testCase.ki_godard, testCase.kp_godard);
                    modeLabel = 'FB';
                otherwise
                    error('Unknown mode: %s', mode);
            end

            % --- Downsample to symbol rate ---
            crSym = crOut(1:SpS_:end);

            % Discard half a block at each end (FF: extrapolation; FB:
            % loop-filter transient) for the BER calculation.
            skipSym  = ceil(testCase.N_fft / (2 * SpS_));
            crSymBER = crSym(skipSym+1:end-skipSym);

            % Resolve phase ambiguity (BER on settled portion only)
            txRefBits = modem.symbolsToBits(symbols);
            skipBits  = skipSym * 2;        % QPSK: 2 bits per symbol
            [BER, crSymBER] = bestRotationBER(testCase, crSymBER, ...
                txRefBits(skipBits+1:end-skipBits));

            fprintf('Godard %s %s BER = %.2e\n', modeLabel, impLabel, BER);

            % --- Plot ---
            rxSym = rx2(1:SpS_:end);
            plotBeforeAfter(testCase, rxSym, crSymBER, ...
                sprintf('Godard %s  |  %s', modeLabel, impLabel), BER);

            % --- Verify ---
            testCase.verifyTrue(all(isfinite(crSym(:))), ...
                'Clock-recovery output contains NaN/Inf.');
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD, ...
                sprintf('Godard %s %s BER %.2e exceeds threshold.', ...
                modeLabel, impLabel, BER));
        end

        function [BER, bestSym] = bestRotationBER(~, crSym, txRefBits)
            %BESTROTATIONBER  Try all four pi/2 rotations, return lowest BER.
            bestBER = Inf;
            bestSym = crSym;
            for kk = 0:3
                rotated = crSym .* exp(-1j * kk * pi/2);
                decided = modem.decideSymbols(rotated);
                rxBits  = modem.symbolsToBits(decided);
                nBits   = min(length(rxBits), length(txRefBits));
                nErr    = sum(txRefBits(1:nBits) ~= rxBits(1:nBits));
                thisBER = nErr / nBits;
                if thisBER < bestBER
                    bestBER = thisBER;
                    bestSym = rotated;
                end
            end
            BER = bestBER;
        end

        function plotBeforeAfter(~, rxSym, crSym, titleStr, BER)
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

            sgtitle(sprintf('QPSK  |  %s  |  BER = %.2e', titleStr, BER));
        end
    end
end
