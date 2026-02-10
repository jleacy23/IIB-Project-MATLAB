classdef test_CDEqualizer < matlab.unittest.TestCase
    % Tests for CD equalizer: CD-only BER check and CD+AWGN visual test.

    properties (Constant)
        % Modulation
        M       = 16
        N_pol   = 2
        Ns      = 4096          % symbols per polarisation
        SpS     = 2             % samples per symbol

        % System
        Rs      = 32            % symbol rate [GBd]
        L       = 80            % fibre length [km]
        D       = 17            % dispersion coeff [ps/nm/km]
        CWL     = 1550          % central wavelength [nm]

        % Channel (unused impairments set to benign values)
        DGDSpec = 0             % no PMD
        N_pmd   = 1
        LW      = 0             % no phase noise

        % CD Equalizer
        NFFT    = 512

        % Thresholds
        BER_CD_ONLY = 1e-3
    end

    methods (TestMethodSetup)
        function seedRng(~)
            rng(42);
        end
    end

    % ----------------------------------------------------------------
    %  Test 1 — CD only: BER must be below threshold
    % ----------------------------------------------------------------
    methods (Test)
        function testCDOnlyBER(testCase)
            % --- Tx ---
            k     = log2(testCase.M);
            Nbits = k * testCase.N_pol * testCase.Ns;
            txBits   = qam_randomBits(Nbits);
            symbols  = qam_modulate(txBits, testCase.M, testCase.N_pol);
            txSig    = qam_rectPulse(symbols, testCase.SpS);

            % --- Channel (CD only, high SNR, no phase noise/PMD) ---
            rxSig = channel_add_chromatic_dispersion(txSig, ...
                testCase.L, testCase.SpS, testCase.Rs, testCase.D, testCase.CWL);

            % --- CD Equalizer ---
            eqSig = cdeq_equalize(rxSig, testCase.D, testCase.L, ...
                testCase.CWL, testCase.Rs, testCase.N_pol, ...
                testCase.SpS, testCase.NFFT);

            % --- Downsample & recover bits ---
            eqSymbols    = eqSig(1:testCase.SpS:end, :);
            decidedSyms  = qam_decideSymbols(eqSymbols, testCase.M, testCase.N_pol);
            rxBits       = qam_symbolsToBits(decidedSyms, testCase.M);

            % --- BER ---
            nErrors = sum(txBits ~= rxBits);
            BER     = nErrors / length(txBits);
            fprintf('CD-only BER = %.2e  (%d errors / %d bits)\n', ...
                     BER, nErrors, length(txBits));

            testCase.verifyLessThan(BER, testCase.BER_CD_ONLY, ...
                sprintf('CD-only BER %.2e exceeds threshold %.2e.', ...
                         BER, testCase.BER_CD_ONLY));
        end

        % ------------------------------------------------------------
        %  Test 2 — CD + AWGN: visual constellation comparison
        % ------------------------------------------------------------
        function testCDPlusAWGNConstellation(testCase)
            SNR_dB = 25;

            % --- Tx ---
            k     = log2(testCase.M);
            Nbits = k * testCase.N_pol * testCase.Ns;
            txBits   = qam_randomBits(Nbits);
            symbols  = qam_modulate(txBits, testCase.M, testCase.N_pol);
            txSig    = qam_rectPulse(symbols, testCase.SpS);

            % --- Channel (CD + AWGN) ---
            rxSig = channel_add_chromatic_dispersion(txSig, ...
                testCase.L, testCase.SpS, testCase.Rs, testCase.D, testCase.CWL);
            rxSig = channel_add_awgn(rxSig, SNR_dB);

            % --- CD Equalizer ---
            eqSig = cdeq_equalize(rxSig, testCase.D, testCase.L, ...
                testCase.CWL, testCase.Rs, testCase.N_pol, ...
                testCase.SpS, testCase.NFFT);

            % --- Downsample ---
            rxSymbols = rxSig(1:testCase.SpS:end, :);
            eqSymbols = eqSig(1:testCase.SpS:end, :);

            % --- BER after equalization ---
            decidedSyms = qam_decideSymbols(eqSymbols, testCase.M, testCase.N_pol);
            rxBits      = qam_symbolsToBits(decidedSyms, testCase.M);
            nErrors = sum(txBits ~= rxBits);
            BER     = nErrors / length(txBits);
            fprintf('CD+AWGN BER = %.2e  (%d errors / %d bits)\n', ...
                     BER, nErrors, length(txBits));

            % --- Plot before / after CD equalization ---
            figure('Name','CD+AWGN Equalization Test', ...
                   'Position',[100 100 1200 500]);

            for p = 1:testCase.N_pol
                % Before equalization
                subplot(2, 2, (p-1)*2 + 1);
                plot(real(rxSymbols(:,p)), imag(rxSymbols(:,p)), '.', ...
                     'MarkerSize', 2);
                grid on; axis equal;
                title(sprintf('Before CD EQ  –  Pol %d', p));
                xlabel('In-Phase'); ylabel('Quadrature');

                % After equalization
                subplot(2, 2, (p-1)*2 + 2);
                plot(real(eqSymbols(:,p)), imag(eqSymbols(:,p)), '.', ...
                     'MarkerSize', 2);
                grid on; axis equal;
                title(sprintf('After CD EQ  –  Pol %d', p));
                xlabel('In-Phase'); ylabel('Quadrature');
            end
            sgtitle(sprintf('16-QAM: CD + AWGN (%d dB)  |  BER = %.2e', ...
                             SNR_dB, BER));

            % Sanity checks
            testCase.verifyTrue(all(isfinite(eqSig(:))), ...
                'Equalizer output contains NaN/Inf.');
            testCase.verifySize(eqSig(1:testCase.SpS:end, :), ...
                size(symbols), ...
                'Equalized symbol array shape mismatch.');
        end

    end
end
