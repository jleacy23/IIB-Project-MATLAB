classdef test_Channel < matlab.unittest.TestCase
    % Visual and shape tests for the Channel functions.
    % Generates 16-QAM constellations and passes them through various
    % impairment combinations, plotting the result for visual inspection.

    properties (Constant)
        % Simulation parameters
        M       = 16
        N_pol   = 2
        Ns      = 4096          % symbols per polarization
        SpS     = 2             % samples per symbol
        Rs      = 32            % symbol rate [GBd]
        L       = 80            % fibre length [km]
        SNR_dB  = 25            % SNR [dB]
        D       = 17            % dispersion coeff [ps/nm/km]
        CWL     = 1550          % central wavelength [nm]
        DGDSpec = 0.1           % PMD coeff [ps/sqrt(km)]
        N_pmd   = 10            % PMD stages
        LW      = 100e3         % laser linewidth [Hz]
    end

    methods (TestMethodSetup)
        function seedRng(~)
            rng(42);            % reproducible noise draws
        end
    end

    % ----- helper ---------------------------------------------------
    methods (Access = private)
        function [txSig, symbols] = generateTx(testCase)
            k       = log2(testCase.M);
            Nbits   = k * testCase.N_pol * testCase.Ns;
            bits    = qam_randomBits(Nbits);
            symbols = qam_modulate(bits, testCase.M, testCase.N_pol);
            txSig   = qam_rectPulse(symbols, testCase.SpS);
        end

        function sym = downsample(testCase, sig)
            % Extract symbol-rate samples from the oversampled signal
            sym = sig(1:testCase.SpS:end, :);
        end
    end

    % ====== Shape-check tests =======================================
    methods (Test)
        function testAWGNOutputShape(testCase)
            [txSig, ~] = generateTx(testCase);
            rxSig = channel_add_awgn(txSig, testCase.SNR_dB);
            testCase.verifySize(rxSig, size(txSig), ...
                'AWGN output size must match input size.');
        end

        function testPhaseNoiseOutputShape(testCase)
            [txSig, ~] = generateTx(testCase);
            rxSig = channel_add_phase_noise(txSig, testCase.Rs, testCase.LW);
            testCase.verifySize(rxSig, size(txSig), ...
                'Phase-noise output size must match input size.');
        end

        function testCDOutputShape(testCase)
            [txSig, ~] = generateTx(testCase);
            rxSig = channel_add_chromatic_dispersion(txSig, ...
                testCase.L, testCase.SpS, testCase.Rs, testCase.D, testCase.CWL);
            testCase.verifySize(rxSig, size(txSig), ...
                'CD output size must match input size.');
        end

        function testPMDOutputShape(testCase)
            [txSig, ~] = generateTx(testCase);
            rxSig = channel_add_pmd(txSig, testCase.L, testCase.SpS, ...
                testCase.Rs, testCase.DGDSpec, testCase.N_pmd);
            testCase.verifySize(rxSig, size(txSig), ...
                'PMD output size must match input size.');
        end

        function testAllImpairOutputShape(testCase)
            [txSig, ~] = generateTx(testCase);
            rxSig = channel_add_awgn(txSig, testCase.SNR_dB);
            rxSig = channel_add_phase_noise(rxSig, testCase.Rs, testCase.LW);
            rxSig = channel_add_chromatic_dispersion(rxSig, ...
                testCase.L, testCase.SpS, testCase.Rs, testCase.D, testCase.CWL);
            rxSig = channel_add_pmd(rxSig, testCase.L, testCase.SpS, ...
                testCase.Rs, testCase.DGDSpec, testCase.N_pmd);
            testCase.verifySize(rxSig, size(txSig), ...
                'Combined-impairment output size must match input size.');
        end
    end

    % ====== Visual constellation tests ==============================
    methods (Test)
        function testVisualConstellations(testCase)
            [txSig, symbols] = generateTx(testCase);

            % --- 1. AWGN only ---
            rng(42);
            rx_awgn = channel_add_awgn(txSig, testCase.SNR_dB);
            sym_awgn = downsample(testCase, rx_awgn);

            % --- 2. AWGN + phase noise ---
            rng(42);
            rx_pn = channel_add_awgn(txSig, testCase.SNR_dB);
            rx_pn = channel_add_phase_noise(rx_pn, testCase.Rs, testCase.LW);
            sym_pn = downsample(testCase, rx_pn);

            % --- 3. AWGN + chromatic dispersion ---
            rng(42);
            rx_cd = channel_add_awgn(txSig, testCase.SNR_dB);
            rx_cd = channel_add_chromatic_dispersion(rx_cd, ...
                testCase.L, testCase.SpS, testCase.Rs, testCase.D, testCase.CWL);
            sym_cd = downsample(testCase, rx_cd);

            % --- 4. AWGN + PMD ---
            rng(42);
            rx_pmd = channel_add_awgn(txSig, testCase.SNR_dB);
            rx_pmd = channel_add_pmd(rx_pmd, testCase.L, testCase.SpS, ...
                testCase.Rs, testCase.DGDSpec, testCase.N_pmd);
            sym_pmd = downsample(testCase, rx_pmd);

            % --- 5. All combined ---
            rng(42);
            rx_all = channel_add_awgn(txSig, testCase.SNR_dB);
            rx_all = channel_add_phase_noise(rx_all, testCase.Rs, testCase.LW);
            rx_all = channel_add_chromatic_dispersion(rx_all, ...
                testCase.L, testCase.SpS, testCase.Rs, testCase.D, testCase.CWL);
            rx_all = channel_add_pmd(rx_all, testCase.L, testCase.SpS, ...
                testCase.Rs, testCase.DGDSpec, testCase.N_pmd);
            sym_all = downsample(testCase, rx_all);

            % --- Plot ---
            labels = {'AWGN', 'AWGN + Phase Noise', ...
                      'AWGN + CD', 'AWGN + PMD', 'All Impairments'};
            data   = {sym_awgn, sym_pn, sym_cd, sym_pmd, sym_all};

            figure('Name','Channel Visual Test','Position',[100 100 1200 900]);
            for k = 1:numel(data)
                for p = 1:testCase.N_pol
                    idx = (k-1)*testCase.N_pol + p;
                    subplot(5, testCase.N_pol, idx);
                    plot(real(data{k}(:,p)), imag(data{k}(:,p)), '.', ...
                         'MarkerSize', 2);
                    grid on; axis equal;
                    title(sprintf('%s  –  Pol %d', labels{k}, p));
                    xlabel('In-Phase'); ylabel('Quadrature');
                end
            end
            sgtitle('16-QAM Constellations under Channel Impairments');

            % Verify no NaNs or Infs in any output
            for k = 1:numel(data)
                testCase.verifyTrue(all(isfinite(data{k}(:))), ...
                    sprintf('%s produced NaN/Inf values.', labels{k}));
            end
        end
    end
end
