classdef results_FreqRecovery < matlab.unittest.TestCase
    %RESULTS_FREQRECOVERY  Examine the effect of pulse shaping on the
    %  4th-power frequency-offset spectrum used by cr_freq_recovery.
    %
    %  For each pulse type (rect, RRC) the test:
    %    1. QAM-modulates random symbols and pulse-shapes them.
    %    2. Applies AWGN + a known frequency offset via channel_lo_freq_shift.
    %    3. Applies matching matched filtering + downsampling.
    %    4. Raises the received symbols to the 4th power and takes the FFT.
    %    5. Finds the peak frequency and compares with the true offset.
    %    6. Plots the 4th-power spectra for rect vs RRC side by side.

    properties (Constant)
        M       = 4              % QPSK
        N_pol   = 2
        Ns      = 2^16           % symbols per polarisation
        SpS     = 2              % samples per symbol
        Rs      = 32             % symbol rate [GBd]
        SNR_dB  = 25             % [dB]
        DeltaF  = 800            % frequency offset [MHz]

        % RRC parameters
        Rolloff = 0.1
        Span    = 10             % filter span in symbols
    end

    methods (TestMethodSetup)
        function seedRng(~)
            rng(42);
        end
    end

    methods (Test)

        function testPulseShapingEffect(testCase)
            M_      = testCase.M;
            Ns_     = testCase.Ns;
            N_pol_  = testCase.N_pol;
            SpS_    = testCase.SpS;
            Rs_     = testCase.Rs;          % GBd
            DeltaF_ = testCase.DeltaF;      % MHz
            SNR_    = testCase.SNR_dB;

            % --- Tx: generate symbols ---
            k      = log2(M_);
            Nbits  = k * N_pol_ * Ns_;
            txBits = modem.randomBits(Nbits);
            symbols = modem.modulate(txBits, M_, N_pol_);

            % =============================================================
            %  Path A: no pulse shaping (symbol-rate)
            % =============================================================
            rxNone = channel.add_awgn(symbols, SNR_);
            rxNone = channel.lo_freq_shift(rxNone, DeltaF_, Rs_, 1);
            symNone = rxNone;   % already at symbol rate

            % =============================================================
            %  Path B: rectangular pulse shaping
            % =============================================================
            txRect = modem.rectPulse(symbols, SpS_);
            rxRect = channel.add_awgn(txRect, SNR_);
            rxRect = channel.lo_freq_shift(rxRect, DeltaF_, Rs_, SpS_);

            % Matched filter + downsample
            rxRectFilt = modem.matched_filter(rxRect, SpS_, 'rect');
            symRect    = rxRectFilt(1:SpS_:end, :);

            % =============================================================
            %  Path C: RRC pulse shaping
            % =============================================================
            txRRC = modem.rrcPulse(symbols, SpS_, testCase.Rolloff, testCase.Span);
            rxRRC = channel.add_awgn(txRRC, SNR_);
            rxRRC = channel.lo_freq_shift(rxRRC, DeltaF_, Rs_, SpS_);

            % Matched filter + downsample
            rxRRCFilt = modem.matched_filter(rxRRC, SpS_, 'rrc', ...
                testCase.Rolloff, testCase.Span);
            symRRC    = rxRRCFilt(1:SpS_:end, :);

            % =============================================================
            %  4th-power spectrum analysis (pol 1 only)
            % =============================================================
            NsOut = size(symRect, 1);
            rs_Hz = Rs_ * 1e9;                                   % Hz
            fAxis = (-1/2 + 1/NsOut : 1/NsOut : 1/2) * rs_Hz;   % Hz

            specNone = fftshift(abs(fft(symNone(:,1).^4)));
            specRect = fftshift(abs(fft(symRect(:,1).^4)));
            specRRC  = fftshift(abs(fft(symRRC(:,1).^4)));

            % Restrict to positive frequencies
            posIdx = fAxis >= 0;
            snP = specNone; snP(~posIdx) = 0;
            srP = specRect; srP(~posIdx) = 0;
            ssP = specRRC;  ssP(~posIdx) = 0;

            [~, idxNone] = max(snP);
            [~, idxRect] = max(srP);
            [~, idxRRC]  = max(ssP);

            estNone_Hz = fAxis(idxNone) / 4;   % undo 4th power
            estRect_Hz = fAxis(idxRect) / 4;
            estRRC_Hz  = fAxis(idxRRC)  / 4;

            trueF_Hz = DeltaF_ * 1e6;          % MHz -> Hz
            freqRes  = rs_Hz / NsOut;           % spectral resolution

            fprintf('  True offset       = %.2f MHz\n', DeltaF_);
            fprintf('  None estimate     = %.2f MHz  (error %.2f kHz)\n', ...
                     estNone_Hz/1e6, abs(estNone_Hz - trueF_Hz)/1e3);
            fprintf('  Rect estimate     = %.2f MHz  (error %.2f kHz)\n', ...
                     estRect_Hz/1e6, abs(estRect_Hz - trueF_Hz)/1e3);
            fprintf('  RRC  estimate     = %.2f MHz  (error %.2f kHz)\n', ...
                     estRRC_Hz/1e6,  abs(estRRC_Hz  - trueF_Hz)/1e3);
            fprintf('  Spectral res      = %.2f kHz\n', freqRes/1e3);

            % --- Verify all estimates are within one spectral bin ---
            testCase.verifyLessThan(abs(estNone_Hz - trueF_Hz), freqRes, ...
                sprintf('None peak error %.2f Hz exceeds resolution %.2f Hz.', ...
                         abs(estNone_Hz - trueF_Hz), freqRes));
            testCase.verifyLessThan(abs(estRect_Hz - trueF_Hz), freqRes, ...
                sprintf('Rect peak error %.2f Hz exceeds resolution %.2f Hz.', ...
                         abs(estRect_Hz - trueF_Hz), freqRes));
            testCase.verifyLessThan(abs(estRRC_Hz - trueF_Hz), freqRes, ...
                sprintf('RRC peak error %.2f Hz exceeds resolution %.2f Hz.', ...
                         abs(estRRC_Hz - trueF_Hz), freqRes));

            % =============================================================
            %  Plot: 4th-power spectra comparison
            % =============================================================
            fAxis_MHz = fAxis / 1e6;

            figure('Name', 'Pulse Shaping Effect on 4th-Power Spectrum', ...
                   'Position', [100 100 1600 500]);

            subplot(1, 3, 1);
            plot(fAxis_MHz, 10*log10(specNone / max(specNone)), 'LineWidth', 0.8);
            hold on;
            xline(4*DeltaF_, 'r--', 'LineWidth', 1.2);
            xline(-4*DeltaF_, 'r--', 'LineWidth', 1.2);
            hold off;
            xlabel('Frequency [MHz]'); ylabel('Normalised PSD [dB]');
            title('No pulse shaping');
            xlim([-Rs_*1e3/2, Rs_*1e3/2]);
            ylim([-60 0]);
            grid on;
            legend('|FFT\{x^4\}|', '4\Deltaf (true)', 'Location', 'south');

            subplot(1, 3, 2);
            plot(fAxis_MHz, 10*log10(specRect / max(specRect)), 'LineWidth', 0.8);
            hold on;
            xline(4*DeltaF_, 'r--', 'LineWidth', 1.2);
            xline(-4*DeltaF_, 'r--', 'LineWidth', 1.2);
            hold off;
            xlabel('Frequency [MHz]'); ylabel('Normalised PSD [dB]');
            title('Rectangular pulse shaping');
            xlim([-Rs_*1e3/2, Rs_*1e3/2]);
            ylim([-60 0]);
            grid on;
            legend('|FFT\{x^4\}|', '4\Deltaf (true)', 'Location', 'south');

            subplot(1, 3, 3);
            plot(fAxis_MHz, 10*log10(specRRC / max(specRRC)), 'LineWidth', 0.8);
            hold on;
            xline(4*DeltaF_, 'r--', 'LineWidth', 1.2);
            xline(-4*DeltaF_, 'r--', 'LineWidth', 1.2);
            hold off;
            xlabel('Frequency [MHz]'); ylabel('Normalised PSD [dB]');
            title(sprintf('RRC pulse shaping (\\beta = %.2f)', testCase.Rolloff));
            xlim([-Rs_*1e3/2, Rs_*1e3/2]);
            ylim([-60 0]);
            grid on;
            legend('|FFT\{x^4\}|', '4\Deltaf (true)', 'Location', 'south');

            sgtitle(sprintf('4th-Power Spectrum  |  %d-QAM  |  \\Deltaf = %d MHz  |  SNR = %d dB', ...
                             M_, DeltaF_, SNR_));
        end

    end
end
