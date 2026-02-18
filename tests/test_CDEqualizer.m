classdef test_CDEqualizer < matlab.unittest.TestCase
    % Tests for CD equalizer: CD-only BER check and CD+AWGN visual test.

    properties (Constant)
        % Modulation
        M       = 16
        N_pol   = 2
        Ns      = 1024          % symbols per polarisation
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

        % ============================================================
        %  Fixed-point tests
        % ============================================================

        % -------- CD only: fxp32 BER check ----------------------------
        function testCDOnlyBER_Fxp32(testCase)
            T = cdeq_equalize_fxp_types('fixed16');

            % --- Tx ---
            k     = log2(testCase.M);
            Nbits = k * testCase.N_pol * testCase.Ns;
            txBits   = qam_randomBits(Nbits);
            symbols  = qam_modulate(txBits, testCase.M, testCase.N_pol);
            txSig    = qam_rectPulse(symbols, testCase.SpS);

            % --- Channel (CD only) ---
            rxSig = channel_add_chromatic_dispersion(txSig, ...
                testCase.L, testCase.SpS, testCase.Rs, testCase.D, testCase.CWL);

            % --- Cast to fi ---
            rxSig_fi = cast(rxSig, 'like', T.x);

            % --- CD Equalizer (fxp MATLAB) ---
            eqSig = cdeq_equalize_fxp(rxSig_fi, testCase.D, testCase.L, ...
                testCase.CWL, testCase.Rs, testCase.N_pol, testCase.SpS, ...
                testCase.NFFT, false, T)

            % --- BER ---
            eqSymbols   = double(eqSig(1:testCase.SpS:end, :));
            decidedSyms = qam_decideSymbols(eqSymbols, testCase.M, testCase.N_pol);
            rxBits      = qam_symbolsToBits(decidedSyms, testCase.M);
            nErrors = sum(txBits ~= rxBits);
            BER     = nErrors / length(txBits);
            fprintf('CD-only FXP32 BER = %.2e  (%d errors / %d bits)\n', ...
                     BER, nErrors, length(txBits));

            testCase.verifyLessThan(BER, testCase.BER_CD_ONLY, ...
                sprintf('FXP32 CD-only BER %.2e exceeds threshold %.2e.', ...
                         BER, testCase.BER_CD_ONLY));
        end

        % -------- CD only: fxp32 vs float NRMSE ----------------------
        function testCDOnly_Fxp32_vs_Float(testCase)
            T = cdeq_equalize_fxp_types('fixed16');

            % --- Tx ---
            k     = log2(testCase.M);
            Nbits = k * testCase.N_pol * testCase.Ns;
            txBits   = qam_randomBits(Nbits);
            symbols  = qam_modulate(txBits, testCase.M, testCase.N_pol);
            txSig    = qam_rectPulse(symbols, testCase.SpS);

            % --- Channel (CD only) ---
            rxSig = channel_add_chromatic_dispersion(txSig, ...
                testCase.L, testCase.SpS, testCase.Rs, testCase.D, testCase.CWL);

            % --- Float reference ---
            eqRef = cdeq_equalize(rxSig, testCase.D, testCase.L, ...
                testCase.CWL, testCase.Rs, testCase.N_pol, testCase.SpS, ...
                testCase.NFFT);

            % --- FXP ---
            rxSig_fi = cast(rxSig, 'like', T.x);
            eqFxp    = cdeq_equalize_fxp(rxSig_fi, testCase.D, testCase.L, ...
                testCase.CWL, testCase.Rs, testCase.N_pol, testCase.SpS, ...
                testCase.NFFT, false, T);

            % --- NRMSE ---
            nrmse = norm(double(eqFxp) - eqRef) / norm(eqRef);
            fprintf('CD-only FXP32 vs float NRMSE = %.4e\n', nrmse);

            testCase.verifyLessThan(nrmse, 0.05, ...
                sprintf('FXP32 NRMSE %.4e exceeds 5%% threshold.', nrmse));
        end

        % -------- CD only: fxp MEX bit-exact with MATLAB fxp ----------
        function testCDOnly_Fxp32_MexMatch(testCase)
            testCase.assumeTrue(exist('cdeq_equalize_fxp_mex', 'file') == 3, ...
                'cdeq_equalize_fxp_mex not found — run build_cdeq_equalize_fxp_mex first.');

            T = cdeq_equalize_fxp_types('fixed16');

            % --- Tx ---
            k     = log2(testCase.M);
            Nbits = k * testCase.N_pol * testCase.Ns;
            txBits   = qam_randomBits(Nbits);
            symbols  = qam_modulate(txBits, testCase.M, testCase.N_pol);
            txSig    = qam_rectPulse(symbols, testCase.SpS);

            % --- Channel (CD only) ---
            rxSig = channel_add_chromatic_dispersion(txSig, ...
                testCase.L, testCase.SpS, testCase.Rs, testCase.D, testCase.CWL);

            rxSig_fi = cast(rxSig, 'like', T.x);

            % --- MATLAB fxp ---
            eqML = cdeq_equalize_fxp(rxSig_fi, testCase.D, testCase.L, ...
                testCase.CWL, testCase.Rs, testCase.N_pol, testCase.SpS, ...
                testCase.NFFT, false, T);

            % --- MEX fxp ---
            eqMEX = cdeq_equalize_fxp_mex(rxSig_fi, testCase.D, testCase.L, ...
                testCase.CWL, testCase.Rs, testCase.N_pol, testCase.SpS, ...
                testCase.NFFT, false, T);

            % --- Verify bit-exact ---
            testCase.verifyEqual(double(eqMEX), double(eqML), ...
                'MEX output must be bit-exact with MATLAB fxp output.');

            if isa(eqML, 'embedded.fi') && isa(eqMEX, 'embedded.fi')
                testCase.verifyEqual(eqMEX.WordLength, eqML.WordLength, ...
                    'MEX WordLength differs from MATLAB.');
                testCase.verifyEqual(eqMEX.FractionLength, eqML.FractionLength, ...
                    'MEX FractionLength differs from MATLAB.');
            end
        end

        % -------- CD + AWGN: fxp32 visual comparison -----------------
        function testCDPlusAWGN_Fxp32(testCase)
            SNR_dB = 25;
            T = cdeq_equalize_fxp_types('fixed16');

            % --- Tx ---
            k     = log2(testCase.M);
            Nbits = k * testCase.N_pol * testCase.Ns;
            txBits   = qam_randomBits(Nbits);
            symbols  = qam_modulate(txBits, testCase.M, testCase.N_pol);
            txSig    = qam_rectPulse(symbols, testCase.SpS);

            % --- Channel ---
            rxSig = channel_add_chromatic_dispersion(txSig, ...
                testCase.L, testCase.SpS, testCase.Rs, testCase.D, testCase.CWL);
            rxSig = channel_add_awgn(rxSig, SNR_dB);

            % --- Float reference ---
            eqRef = cdeq_equalize(rxSig, testCase.D, testCase.L, ...
                testCase.CWL, testCase.Rs, testCase.N_pol, testCase.SpS, ...
                testCase.NFFT);

            % --- FXP ---
            rxSig_fi = cast(rxSig, 'like', T.x);
            eqFxp    = cdeq_equalize_fxp(rxSig_fi, testCase.D, testCase.L, ...
                testCase.CWL, testCase.Rs, testCase.N_pol, testCase.SpS, ...
                testCase.NFFT, false, T);

            eqRefSym = eqRef(1:testCase.SpS:end, :);
            eqFxpSym = double(eqFxp(1:testCase.SpS:end, :));
            rxSym    = rxSig(1:testCase.SpS:end, :);

            % --- Plot ---
            figure('Name', 'CD+AWGN FXP32 Comparison', ...
                   'Position', [100 100 1400 700]);
            for p = 1:testCase.N_pol
                subplot(2, 3, (p-1)*3 + 1);
                plot(real(rxSym(:,p)), imag(rxSym(:,p)), '.', 'MarkerSize', 2);
                grid on; axis equal;
                title(sprintf('Before CD EQ – Pol %d', p));

                subplot(2, 3, (p-1)*3 + 2);
                plot(real(eqRefSym(:,p)), imag(eqRefSym(:,p)), '.', 'MarkerSize', 2);
                grid on; axis equal;
                title(sprintf('Float CD EQ – Pol %d', p));

                subplot(2, 3, (p-1)*3 + 3);
                plot(real(eqFxpSym(:,p)), imag(eqFxpSym(:,p)), '.', 'MarkerSize', 2);
                grid on; axis equal;
                title(sprintf('FXP32 CD EQ – Pol %d', p));
            end
            sgtitle(sprintf('16-QAM: CD + AWGN (%d dB)  |  Float vs FXP32', SNR_dB));

            % Sanity checks
            testCase.verifyTrue(all(isfinite(double(eqFxp(:)))), ...
                'FXP32 output contains NaN/Inf.');
        end

    end
end
