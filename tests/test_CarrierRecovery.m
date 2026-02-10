classdef test_CarrierRecovery < matlab.unittest.TestCase
    % Tests for carrier recovery (Viterbi-Viterbi).
    % Applies AWGN + phase noise, runs VV carrier recovery, plots
    % before/after constellations and checks BER.

    properties (Constant)
        N_pol   = 2
        Ns      = 8192          % symbols per polarisation
        SpS     = 1             % symbol-rate processing (no pulse shaping)

        % System
        Rs      = 32            % [GBd]
        SNR_dB  = 20            % [dB]
        Linewidth = 100e4       % laser linewidth [Hz]

        % Channel (unused impairments set to benign values)
        L       = 80            % fibre length [km]
        D       = 0             % no CD
        CWL     = 1550          % [nm]
        DGDSpec = 0             % no PMD
        N_pmd   = 1
        LW      = 100e4         % phase-noise linewidth [Hz]

        % Carrier recovery
        NTaps   = 15

        % Pass / fail
        BER_THRESHOLD = 5e-2
    end

    methods (TestMethodSetup)
        function seedRng(~)
            rng(42);
        end
    end

    % ================================================================
    methods (Test)

        % -------- 4-QAM -----------------------------
        function testQPSK(testCase)
            M = 4;
            [rxSym, crSym, BER] = runScenario(testCase, M);
            plotBeforeAfter(testCase, rxSym, crSym, ...
                '4-QAM  |  VV', M, BER);

            testCase.verifyTrue(all(isfinite(crSym(:))), ...
                'Carrier-recovery output contains NaN/Inf.');
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD, ...
                sprintf('4-QAM BER %.2e exceeds threshold.', BER));
        end
    end

    % ================================================================
    %  Helpers
    % ================================================================
    methods (Access = private)

        function [rxSym, crSym, BER] = runScenario(testCase, M)
            rng(42);

            % --- Tx ---
            k      = log2(M);
            Nbits  = k * testCase.N_pol * testCase.Ns;
            txBits = qam_randomBits(Nbits);
            symbols = qam_modulate(txBits, M, testCase.N_pol);

            % --- Channel: AWGN + phase noise ---
            rxSym = channel_add_awgn(symbols, testCase.SNR_dB);
            rxSym = channel_add_phase_noise(rxSym, testCase.Rs, testCase.LW);

            % --- Carrier Recovery (Viterbi-Viterbi) ---
            symEnergy = mean(abs(symbols(:)).^2);
            VVFilter  = cr_genVVFilter(testCase.Linewidth, testCase.Rs, ...
                testCase.SNR_dB, symEnergy, testCase.N_pol, testCase.NTaps);
            crSym = cr_viterbiViterbi(rxSym, testCase.N_pol, testCase.NTaps, ...
                VVFilter);

            % --- Demodulate & BER ---
            decidedSyms = qam_decideSymbols(crSym, M, testCase.N_pol);
            rxBits      = qam_symbolsToBits(decidedSyms, M);

            nErrors = sum(txBits ~= rxBits);
            BER     = nErrors / length(txBits);
            fprintf('%d-QAM BER = %.2e  (%d / %d)\n', ...
                     M, BER, nErrors, length(txBits));
        end

        function plotBeforeAfter(testCase, rxSym, crSym, titleStr, M, BER)
            figure('Name', titleStr, 'Position', [100 100 1200 500]);
            for p = 1:testCase.N_pol
                % Before carrier recovery
                subplot(2, 2, (p-1)*2 + 1);
                plot(real(rxSym(:,p)), imag(rxSym(:,p)), '.', ...
                     'MarkerSize', 2);
                grid on; axis equal;
                title(sprintf('Before CR  –  Pol %d', p));
                xlabel('In-Phase'); ylabel('Quadrature');

                % After carrier recovery
                subplot(2, 2, (p-1)*2 + 2);
                plot(real(crSym(:,p)), imag(crSym(:,p)), '.', ...
                     'MarkerSize', 2);
                grid on; axis equal;
                title(sprintf('After CR  –  Pol %d', p));
                xlabel('In-Phase'); ylabel('Quadrature');
            end
            sgtitle(sprintf('%d-QAM: AWGN + Phase Noise  |  %s  |  BER = %.2e', ...
                             M, titleStr, BER));
        end
    end
end
