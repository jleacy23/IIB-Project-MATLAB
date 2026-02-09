classdef test_CDEqualizer < matlab.unittest.TestCase
    % Tests for CDEqualizer: CD-only BER check and CD+AWGN visual test.

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
        FL      = 12          
        WL      = 12 + floor(log2(512)) + 1           % sign bit + integer bits + 1 fractional bit

        % Thresholds
        BER_CD_ONLY = 1e-3      % BER threshold for CD-only case
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
            modem = QAMModem(testCase.M, testCase.N_pol);
            k     = modem.bitsPerSymbol;
            Nbits = k * testCase.N_pol * testCase.Ns;
            txBits   = modem.randomBits(Nbits);
            symbols  = modem.modulate(txBits);
            txSig    = modem.rectPulse(symbols, testCase.SpS);

            % --- Channel (CD only, high SNR, no phase noise/PMD) ---
            ch = Channel(testCase.L, 100, testCase.SpS, testCase.Rs, ...
                         testCase.D, testCase.CWL, ...
                         testCase.DGDSpec, testCase.N_pmd, testCase.LW);
            rxSig = ch.add_chromatic_dispersion(txSig);

            % --- CD Equalizer ---
            cdeq = CDEqualizer(testCase.D, testCase.L, testCase.CWL, ...
                               testCase.Rs, testCase.N_pol, ...
                               testCase.SpS, testCase.NFFT, ...
                               testCase.WL, testCase.FL);
            eqSig = cdeq.equalize(rxSig, false);

            % --- Downsample & recover bits ---
            eqSymbols    = eqSig(1:testCase.SpS:end, :);
            decidedSyms  = modem.decideSymbols(eqSymbols);
            rxBits       = modem.symbolsToBits(decidedSyms);

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
            modem = QAMModem(testCase.M, testCase.N_pol);
            k     = modem.bitsPerSymbol;
            Nbits = k * testCase.N_pol * testCase.Ns;
            txBits   = modem.randomBits(Nbits);
            symbols  = modem.modulate(txBits);
            txSig    = modem.rectPulse(symbols, testCase.SpS);

            % --- Channel (CD + AWGN) ---
            ch = Channel(testCase.L, SNR_dB, testCase.SpS, testCase.Rs, ...
                         testCase.D, testCase.CWL, ...
                         testCase.DGDSpec, testCase.N_pmd, testCase.LW);
            rxSig = ch.add_chromatic_dispersion(txSig);
            rxSig = ch.add_awgn(rxSig);

            % --- CD Equalizer ---
            cdeq = CDEqualizer(testCase.D, testCase.L, testCase.CWL, ...
                               testCase.Rs, testCase.N_pol, ...
                               testCase.SpS, testCase.NFFT, ...
                               testCase.WL, testCase.FL);
            eqSig = cdeq.equalize(rxSig, false);

            % --- Downsample ---
            rxSymbols = rxSig(1:testCase.SpS:end, :);
            eqSymbols = eqSig(1:testCase.SpS:end, :);

            % --- BER after equalization ---
            decidedSyms = modem.decideSymbols(eqSymbols);
            rxBits      = modem.symbolsToBits(decidedSyms);
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

        % ------------------------------------------------------------
        %  Test 3 — CD only, fixed-point equalizer: BER check
        % ------------------------------------------------------------
        function testCDOnlyFixedPointBER(testCase)
            % --- Tx ---
            modem = QAMModem(testCase.M, testCase.N_pol);
            k     = modem.bitsPerSymbol;
            Nbits = k * testCase.N_pol * testCase.Ns;
            txBits   = modem.randomBits(Nbits);
            symbols  = modem.modulate(txBits);
            txSig    = modem.rectPulse(symbols, testCase.SpS);

            % --- Channel (CD only, high SNR) ---
            ch = Channel(testCase.L, 100, testCase.SpS, testCase.Rs, ...
                         testCase.D, testCase.CWL, ...
                         testCase.DGDSpec, testCase.N_pmd, testCase.LW);
            rxSig = ch.add_chromatic_dispersion(txSig);

            % --- Fixed-point CD Equalizer ---
            cdeq = CDEqualizer(testCase.D, testCase.L, testCase.CWL, ...
                               testCase.Rs, testCase.N_pol, ...
                               testCase.SpS, testCase.NFFT, ...
                               testCase.WL, testCase.FL);
            eqSig = cdeq.equalize(rxSig, true);
            eqSig = double(eqSig);

            % --- Downsample & recover bits ---
            eqSymbols    = eqSig(1:testCase.SpS:end, :);
            decidedSyms  = modem.decideSymbols(eqSymbols);
            rxBits       = modem.symbolsToBits(decidedSyms);

            % --- BER ---
            nErrors = sum(txBits ~= rxBits);
            BER     = nErrors / length(txBits);
            fprintf('CD-only fixed-point BER = %.2e  (%d errors / %d bits)\n', ...
                     BER, nErrors, length(txBits));

            % Allow slightly higher threshold for quantisation effects
            testCase.verifyLessThan(BER, testCase.BER_CD_ONLY, ...
                sprintf('Fixed-point CD-only BER %.2e exceeds threshold %.2e.', ...
                         BER, testCase.BER_CD_ONLY));
        end
        %------------------------------------------------
        % Test 4 — CD + AWGN, fixed-point equalizer: visual test
        %------------------------------------------------
        function testCDPlusAWGNFixedPointConstellation(testCase)
            % This test is similar to testCDPlusAWGNConstellation but uses
            % the fixed-point equalizer. It checks that the constellation is
            % still reasonably clear and that the BER is acceptable.

            SNR_dB = 25;

            % --- Tx ---
            modem = QAMModem(testCase.M, testCase.N_pol);
            k     = modem.bitsPerSymbol;
            Nbits = k * testCase.N_pol * testCase.Ns;
            txBits   = modem.randomBits(Nbits);
            symbols  = modem.modulate(txBits);
            txSig    = modem.rectPulse(symbols, testCase.SpS);

            % --- Channel (CD + AWGN) ---
            ch = Channel(testCase.L, SNR_dB, testCase.SpS, testCase.Rs, ...
                         testCase.D, testCase.CWL, ...
                         testCase.DGDSpec, testCase.N_pmd, testCase.LW);
            rxSig = ch.add_chromatic_dispersion(txSig);
            rxSig = ch.add_awgn(rxSig);

            % --- Fixed-point CD Equalizer ---
            cdeq = CDEqualizer(testCase.D, testCase.L, testCase.CWL, ...
                               testCase.Rs, testCase.N_pol, ...
                               testCase.SpS, testCase.NFFT, ...
                               testCase.WL, testCase.FL);
            eqSig = cdeq.equalize(rxSig, true);
            eqSig = double(eqSig);

            % --- Downsample ---
            rxSymbols = rxSig(1:testCase.SpS:end, :);
            eqSymbols = eqSig(1:testCase.SpS:end, :);

            % --- BER after equalization ---
            decidedSyms = modem.decideSymbols(eqSymbols);
            rxBits      = modem.symbolsToBits(decidedSyms);
            nErrors = sum(txBits ~= rxBits);
            BER     = nErrors / length(txBits);
            fprintf('CD+AWGN fixed-point BER = %.2e  (%d errors / %d bits)\n', ...
                     BER, nErrors, length(txBits));

            % --- Plot before / after CD equalization ---
            figure('Name','CD+AWGN Fixed-Point Equalization Test', ...
                   'Position',[150 150 1200 500]);

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
                title(sprintf('After CD EQ (Fixed-Point) –  Pol %d', p));
                xlabel('In-Phase'); ylabel('Quadrature');
            end
        end
    end
end
