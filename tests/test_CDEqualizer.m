classdef test_CDEqualizer < matlab.unittest.TestCase
    % Tests for CD equalizer: CD-only BER check and CD+AWGN visual test.

    properties (Constant)
        % Modulation
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

        % Pulse shaping (Nyquist / raised-cosine)
        Rolloff = 0.25
        Span    = 10

        % Fixed-point configuration
        FxpConfig = 'fixed16'

        % Thresholds
        BER_CD_ONLY = 1e-3
    end

    methods (TestClassSetup)
        function buildFxpMex(testCase)
            % Build the CD fixed-point MEX binaries once so the fxp tests
            % run against the compiled (fast) versions rather than the
            % interpreted fi datapath.  Mirrors test_AdaptiveEqualizer.
            thisDir  = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(thisDir);
            addpath(genpath(fullfile(repoRoot, 'src')));
            addpath(fullfile(repoRoot, 'build'));

            P = struct();
            P.FxpConfig_CD = testCase.FxpConfig;
            P.N_pol        = testCase.N_pol;
            P.D            = testCase.D;
            P.L            = testCase.L;
            P.CWL          = testCase.CWL;
            P.Rs           = testCase.Rs;
            P.SpS          = testCase.SpS;
            P.NFFT         = testCase.NFFT;
            P.po2Twiddle   = false;

            cfg = coder.config('mex');
            build_cd_eq_equalize_fxp_mex(P, cfg);       % freq-domain
            build_cd_eq_equalize_td_fxp_mex(P, cfg);    % time-domain
        end
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
            Nbits    = 4 * testCase.Ns;
            txBits   = modem.randomBits(Nbits);
            symbols  = modem.modulate(txBits);
            txSig    = modem.nyquistPulse(symbols, testCase.SpS, ...
                testCase.Rolloff, testCase.Span);

            % --- Channel (CD only, high SNR, no phase noise/PMD) ---
            rxSig = channel.add_chromatic_dispersion(txSig, ...
                testCase.L, testCase.SpS, testCase.Rs, testCase.D, testCase.CWL);

            % --- CD Equalizer ---
            eqSig = cd_eq.equalize(rxSig, testCase.D, testCase.L, ...
                testCase.CWL, testCase.Rs, testCase.N_pol, ...
                testCase.SpS, testCase.NFFT);

            % --- Downsample & recover bits ---
            eqSymbols    = eqSig(1:testCase.SpS:end, :);
            decidedSyms  = modem.decideSymbols(eqSymbols);
            rxBits       = modem.symbolsToBits(decidedSyms);

            % --- BER ---
            txRefBits = modem.symbolsToBits(symbols);
            nBits   = min(length(txRefBits), length(rxBits));
            nErrors = sum(txRefBits(1:nBits) ~= rxBits(1:nBits));
            BER     = nErrors / nBits;
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
            Nbits    = 4 * testCase.Ns;
            txBits   = modem.randomBits(Nbits);
            symbols  = modem.modulate(txBits);
            txSig    = modem.nyquistPulse(symbols, testCase.SpS, ...
                testCase.Rolloff, testCase.Span);

            % --- Channel (CD + AWGN) ---
            rxSig = channel.add_chromatic_dispersion(txSig, ...
                testCase.L, testCase.SpS, testCase.Rs, testCase.D, testCase.CWL);
            rxSig = channel.add_awgn(rxSig, SNR_dB);

            % --- CD Equalizer ---
            eqSig = cd_eq.equalize(rxSig, testCase.D, testCase.L, ...
                testCase.CWL, testCase.Rs, testCase.N_pol, ...
                testCase.SpS, testCase.NFFT);

            % --- Downsample ---
            rxSymbols = rxSig(1:testCase.SpS:end, :);
            eqSymbols = eqSig(1:testCase.SpS:end, :);

            % --- BER after equalization ---
            decidedSyms = modem.decideSymbols(eqSymbols);
            rxBits      = modem.symbolsToBits(decidedSyms);
            txRefBits   = modem.symbolsToBits(symbols);
            nBits   = min(length(txRefBits), length(rxBits));
            nErrors = sum(txRefBits(1:nBits) ~= rxBits(1:nBits));
            BER     = nErrors / nBits;
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
            sgtitle(sprintf('QPSK: CD + AWGN (%d dB)  |  BER = %.2e', ...
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
            T = cd_eq.equalize_fxp_types('fixed16');

            % --- Tx ---
            Nbits    = 4 * testCase.Ns;
            txBits   = modem.randomBits(Nbits);
            symbols  = modem.modulate(txBits);
            txSig    = modem.nyquistPulse(symbols, testCase.SpS, ...
                testCase.Rolloff, testCase.Span);

            % --- Channel (CD only) ---
            rxSig = channel.add_chromatic_dispersion(txSig, ...
                testCase.L, testCase.SpS, testCase.Rs, testCase.D, testCase.CWL);

            % --- Cast to fi ---
            rxSig_fi = cast(rxSig, 'like', T.x);

            % --- CD Equalizer (fxp MEX) ---
            eqSig = cd_eq.equalize_fxp_mex(rxSig_fi, testCase.D, testCase.L, ...
                testCase.CWL, testCase.Rs, testCase.N_pol, testCase.SpS, ...
                testCase.NFFT, false, T);

            % --- BER ---
            eqSymbols   = double(eqSig(1:testCase.SpS:end, :));
            decidedSyms = modem.decideSymbols(eqSymbols);
            rxBits      = modem.symbolsToBits(decidedSyms);
            txRefBits   = modem.symbolsToBits(symbols);
            nBits   = min(length(txRefBits), length(rxBits));
            nErrors = sum(txRefBits(1:nBits) ~= rxBits(1:nBits));
            BER     = nErrors / nBits;
            fprintf('CD-only FXP32 BER = %.2e  (%d errors / %d bits)\n', ...
                     BER, nErrors, length(txBits));

            testCase.verifyLessThan(BER, testCase.BER_CD_ONLY, ...
                sprintf('FXP32 CD-only BER %.2e exceeds threshold %.2e.', ...
                         BER, testCase.BER_CD_ONLY));
        end

        % -------- CD only: fxp32 vs float NRMSE ----------------------
        function testCDOnly_Fxp32_vs_Float(testCase)
            T = cd_eq.equalize_fxp_types('fixed16');

            % --- Tx ---
            Nbits    = 4 * testCase.Ns;
            txBits   = modem.randomBits(Nbits);
            symbols  = modem.modulate(txBits);
            txSig    = modem.nyquistPulse(symbols, testCase.SpS, ...
                testCase.Rolloff, testCase.Span);

            % --- Channel (CD only) ---
            rxSig = channel.add_chromatic_dispersion(txSig, ...
                testCase.L, testCase.SpS, testCase.Rs, testCase.D, testCase.CWL);

            % --- Float reference ---
            eqRef = cd_eq.equalize(rxSig, testCase.D, testCase.L, ...
                testCase.CWL, testCase.Rs, testCase.N_pol, testCase.SpS, ...
                testCase.NFFT);

            % --- FXP (MEX) ---
            rxSig_fi = cast(rxSig, 'like', T.x);
            eqFxp    = cd_eq.equalize_fxp_mex(rxSig_fi, testCase.D, testCase.L, ...
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
            testCase.assumeTrue(exist('cd_eq.equalize_fxp_mex', 'file') == 3, ...
                'cd_eq.equalize_fxp_mex not found — run build_cd_eq_equalize_fxp_mex first.');

            T = cd_eq.equalize_fxp_types('fixed16');

            % --- Tx ---
            Nbits    = 4 * testCase.Ns;
            txBits   = modem.randomBits(Nbits);
            symbols  = modem.modulate(txBits);
            txSig    = modem.nyquistPulse(symbols, testCase.SpS, ...
                testCase.Rolloff, testCase.Span);

            % --- Channel (CD only) ---
            rxSig = channel.add_chromatic_dispersion(txSig, ...
                testCase.L, testCase.SpS, testCase.Rs, testCase.D, testCase.CWL);

            rxSig_fi = cast(rxSig, 'like', T.x);

            % --- MATLAB fxp ---
            eqML = cd_eq.equalize_fxp(rxSig_fi, testCase.D, testCase.L, ...
                testCase.CWL, testCase.Rs, testCase.N_pol, testCase.SpS, ...
                testCase.NFFT, false, T);

            % --- MEX fxp ---
            eqMEX = cd_eq.equalize_fxp_mex(rxSig_fi, testCase.D, testCase.L, ...
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
            T = cd_eq.equalize_fxp_types('fixed16');

            % --- Tx ---
            Nbits    = 4 * testCase.Ns;
            txBits   = modem.randomBits(Nbits);
            symbols  = modem.modulate(txBits);
            txSig    = modem.nyquistPulse(symbols, testCase.SpS, ...
                testCase.Rolloff, testCase.Span);

            % --- Channel ---
            rxSig = channel.add_chromatic_dispersion(txSig, ...
                testCase.L, testCase.SpS, testCase.Rs, testCase.D, testCase.CWL);
            rxSig = channel.add_awgn(rxSig, SNR_dB);

            % --- Float reference ---
            eqRef = cd_eq.equalize(rxSig, testCase.D, testCase.L, ...
                testCase.CWL, testCase.Rs, testCase.N_pol, testCase.SpS, ...
                testCase.NFFT);

            % --- FXP (MEX) ---
            rxSig_fi = cast(rxSig, 'like', T.x);
            eqFxp    = cd_eq.equalize_fxp_mex(rxSig_fi, testCase.D, testCase.L, ...
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
            sgtitle(sprintf('QPSK: CD + AWGN (%d dB)  |  Float vs FXP32', SNR_dB));

            % Sanity checks
            testCase.verifyTrue(all(isfinite(double(eqFxp(:)))), ...
                'FXP32 output contains NaN/Inf.');
        end

        % ============================================================
        %  Time-domain (FIR) fixed-point tests
        % ============================================================

        % -------- CD only: time-domain fxp BER check -----------------
        function testCDOnlyBER_TD_Fxp(testCase)
            T = cd_eq.equalize_td_fxp_types('fixed16');

            % --- Tx ---
            Nbits    = 4 * testCase.Ns;
            txBits   = modem.randomBits(Nbits);
            symbols  = modem.modulate(txBits);
            txSig    = modem.nyquistPulse(symbols, testCase.SpS, ...
                testCase.Rolloff, testCase.Span);

            % --- Channel (CD only) ---
            rxSig = channel.add_chromatic_dispersion(txSig, ...
                testCase.L, testCase.SpS, testCase.Rs, testCase.D, testCase.CWL);

            % --- Cast to fi ---
            rxSig_fi = cast(rxSig, 'like', T.x);

            % --- Time-domain CD Equalizer (fxp MEX) ---
            NTap  = cd_eq.computeOverlap(testCase.D, testCase.L, testCase.CWL, ...
                testCase.Rs, testCase.SpS, testCase.NFFT);
            eqSig = cd_eq.equalize_td_fxp_mex(rxSig_fi, testCase.D, testCase.L, ...
                testCase.CWL, testCase.Rs, testCase.N_pol, testCase.SpS, NTap, T);

            % --- BER ---
            eqSymbols   = double(eqSig(1:testCase.SpS:end, :));
            decidedSyms = modem.decideSymbols(eqSymbols);
            rxBits      = modem.symbolsToBits(decidedSyms);
            txRefBits   = modem.symbolsToBits(symbols);
            nBits   = min(length(txRefBits), length(rxBits));
            nErrors = sum(txRefBits(1:nBits) ~= rxBits(1:nBits));
            BER     = nErrors / nBits;
            fprintf('CD-only TD-FXP BER = %.2e  (%d errors / %d bits)\n', ...
                     BER, nErrors, length(txBits));

            testCase.verifyLessThan(BER, testCase.BER_CD_ONLY, ...
                sprintf('TD-FXP CD-only BER %.2e exceeds threshold %.2e.', ...
                         BER, testCase.BER_CD_ONLY));
        end

        % -------- CD only: time-domain fxp vs time-domain float ------
        %  Isolates quantization error against the same algorithm's float
        %  reference (equalize_td), so the threshold can be tight.
        function testCDOnly_TD_Fxp_vs_Float(testCase)
            T = cd_eq.equalize_td_fxp_types('fixed16');

            % --- Tx ---
            Nbits    = 4 * testCase.Ns;
            txBits   = modem.randomBits(Nbits);
            symbols  = modem.modulate(txBits);
            txSig    = modem.nyquistPulse(symbols, testCase.SpS, ...
                testCase.Rolloff, testCase.Span);

            % --- Channel (CD only) ---
            rxSig = channel.add_chromatic_dispersion(txSig, ...
                testCase.L, testCase.SpS, testCase.Rs, testCase.D, testCase.CWL);

            % --- Time-domain float reference (same algorithm) ---
            NTap  = cd_eq.computeOverlap(testCase.D, testCase.L, testCase.CWL, ...
                testCase.Rs, testCase.SpS, testCase.NFFT);
            eqRef = cd_eq.equalize_td(rxSig, testCase.D, testCase.L, ...
                testCase.CWL, testCase.Rs, testCase.N_pol, testCase.SpS, NTap);

            % --- Time-domain FXP (MEX) ---
            rxSig_fi = cast(rxSig, 'like', T.x);
            eqFxp    = cd_eq.equalize_td_fxp_mex(rxSig_fi, testCase.D, testCase.L, ...
                testCase.CWL, testCase.Rs, testCase.N_pol, testCase.SpS, NTap, T);

            % --- NRMSE (quantization error only) ---
            nrmse = norm(double(eqFxp) - eqRef) / norm(eqRef);
            fprintf('CD-only TD-FXP vs TD-float NRMSE = %.4e\n', nrmse);

            testCase.verifyLessThan(nrmse, 0.05, ...
                sprintf('TD-FXP NRMSE %.4e exceeds 5%% threshold.', nrmse));
        end

        % -------- CD only: time-domain float vs freq-domain float ----
        %  Cross-checks that the two algorithms (FIR vs overlap-save)
        %  realise the same all-pass response.
        function testCDOnly_TD_vs_FD_Float(testCase)
            % --- Tx ---
            Nbits    = 4 * testCase.Ns;
            txBits   = modem.randomBits(Nbits);
            symbols  = modem.modulate(txBits);
            txSig    = modem.nyquistPulse(symbols, testCase.SpS, ...
                testCase.Rolloff, testCase.Span);

            % --- Channel (CD only) ---
            rxSig = channel.add_chromatic_dispersion(txSig, ...
                testCase.L, testCase.SpS, testCase.Rs, testCase.D, testCase.CWL);

            % --- Freq-domain float (overlap-save) ---
            eqFD = cd_eq.equalize(rxSig, testCase.D, testCase.L, ...
                testCase.CWL, testCase.Rs, testCase.N_pol, testCase.SpS, ...
                testCase.NFFT);

            % --- Time-domain float (FIR) ---
            NTap = cd_eq.computeOverlap(testCase.D, testCase.L, testCase.CWL, ...
                testCase.Rs, testCase.SpS, testCase.NFFT);
            eqTD = cd_eq.equalize_td(rxSig, testCase.D, testCase.L, ...
                testCase.CWL, testCase.Rs, testCase.N_pol, testCase.SpS, NTap);

            % --- NRMSE between the two algorithms ---
            nrmse = norm(eqTD - eqFD) / norm(eqFD);
            fprintf('CD-only TD-float vs FD-float NRMSE = %.4e\n', nrmse);

            testCase.verifyLessThan(nrmse, 0.10, ...
                sprintf('TD vs FD float NRMSE %.4e exceeds 10%% threshold.', nrmse));
        end

        % -------- CD only: time-domain fxp MEX bit-exact -------------
        function testCDOnly_TD_MexMatch(testCase)
            testCase.assumeTrue(exist('cd_eq.equalize_td_fxp_mex', 'file') == 3, ...
                'cd_eq.equalize_td_fxp_mex not found — run build_cd_eq_equalize_td_fxp_mex first.');

            T = cd_eq.equalize_td_fxp_types('fixed16');

            % --- Tx ---
            Nbits    = 4 * testCase.Ns;
            txBits   = modem.randomBits(Nbits);
            symbols  = modem.modulate(txBits);
            txSig    = modem.nyquistPulse(symbols, testCase.SpS, ...
                testCase.Rolloff, testCase.Span);

            % --- Channel (CD only) ---
            rxSig = channel.add_chromatic_dispersion(txSig, ...
                testCase.L, testCase.SpS, testCase.Rs, testCase.D, testCase.CWL);

            rxSig_fi = cast(rxSig, 'like', T.x);

            % --- MATLAB fxp ---
            NTap  = cd_eq.computeOverlap(testCase.D, testCase.L, testCase.CWL, ...
                testCase.Rs, testCase.SpS, testCase.NFFT);
            eqML  = cd_eq.equalize_td_fxp(rxSig_fi, testCase.D, testCase.L, ...
                testCase.CWL, testCase.Rs, testCase.N_pol, testCase.SpS, NTap, T);

            % --- MEX fxp ---
            eqMEX = cd_eq.equalize_td_fxp_mex(rxSig_fi, testCase.D, testCase.L, ...
                testCase.CWL, testCase.Rs, testCase.N_pol, testCase.SpS, NTap, T);

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

        % -------- CD only: time-domain vs freq-domain fxp BER --------
        function testCDOnly_TD_vs_FD_Fxp_BER(testCase)
            T_td = cd_eq.equalize_td_fxp_types('fixed16');
            T_fd = cd_eq.equalize_fxp_types('fixed16');

            % --- Tx ---
            Nbits    = 4 * testCase.Ns;
            txBits   = modem.randomBits(Nbits);
            symbols  = modem.modulate(txBits);
            txSig    = modem.nyquistPulse(symbols, testCase.SpS, ...
                testCase.Rolloff, testCase.Span);

            % --- Channel (CD only) ---
            rxSig = channel.add_chromatic_dispersion(txSig, ...
                testCase.L, testCase.SpS, testCase.Rs, testCase.D, testCase.CWL);

            % --- Freq-domain fxp (MEX) ---
            rxSig_fd = cast(rxSig, 'like', T_fd.x);
            eqFD = cd_eq.equalize_fxp_mex(rxSig_fd, testCase.D, testCase.L, ...
                testCase.CWL, testCase.Rs, testCase.N_pol, testCase.SpS, ...
                testCase.NFFT, false, T_fd);

            % --- Time-domain fxp (MEX) ---
            rxSig_td = cast(rxSig, 'like', T_td.x);
            NTap = cd_eq.computeOverlap(testCase.D, testCase.L, testCase.CWL, ...
                testCase.Rs, testCase.SpS, testCase.NFFT);
            eqTD = cd_eq.equalize_td_fxp_mex(rxSig_td, testCase.D, testCase.L, ...
                testCase.CWL, testCase.Rs, testCase.N_pol, testCase.SpS, NTap, T_td);

            % --- BER for each ---
            txRefBits = modem.symbolsToBits(symbols);

            berFD = localBER(eqFD, testCase.SpS, txRefBits);
            berTD = localBER(eqTD, testCase.SpS, txRefBits);
            fprintf('CD-only FXP BER:  freq-domain = %.2e | time-domain = %.2e\n', ...
                     berFD, berTD);

            testCase.verifyLessThan(berTD, testCase.BER_CD_ONLY, ...
                sprintf('TD-FXP BER %.2e exceeds threshold %.2e.', ...
                         berTD, testCase.BER_CD_ONLY));
        end

    end
end

% --------------------------------------------------------------------
%  Local helper: BER of an equalized signal against reference bits
% --------------------------------------------------------------------
function BER = localBER(eqSig, SpS, txRefBits)
    eqSymbols   = double(eqSig(1:SpS:end, :));
    decidedSyms = modem.decideSymbols(eqSymbols);
    rxBits      = modem.symbolsToBits(decidedSyms);
    nBits   = min(length(txRefBits), length(rxBits));
    nErrors = sum(txRefBits(1:nBits) ~= rxBits(1:nBits));
    BER     = nErrors / nBits;
end
