classdef test_AdaptiveEqualizer < matlab.unittest.TestCase
    % Visual tests for adaptive equalizer: AWGN+PMD then adaptive EQ.

    properties (Constant)
        N_pol   = 2
        Ns      = 8192          % symbols per polarisation
        SpS     = 2

        % System
        Rs      = 32            % [GBd]
        L       = 80            % [km]
        SNR_dB  = 25            % [dB]
        D       = 0             % no CD for this test
        CWL     = 1550          % [nm]
        DGDSpec = 0.5           % PMD coeff [ps/sqrt(km)]
        N_pmd   = 10
        LW      = 0             % no phase noise

        % Adaptive EQ common settings
        NTaps   = 15
        Mu      = 1e-3
        N1      = 2000          % single-spike re-init iteration
        NOut    = 500           % discard transient
    end

    methods (TestMethodSetup)
        function seedRng(~)
            rng(42);
        end
    end

    % ================================================================
    methods (Test)

        % -------- 4-QAM CMA -----------------------
        function testQPSK_CMA(testCase)
            M = 4;
            [rxSym, eqSym] = runScenario(testCase, M, 'CMA', ...
                testCase.NTaps, testCase.Mu, true, testCase.N1, [], ...
                testCase.NOut);
            plotBeforeAfter(testCase, rxSym, eqSym, ...
                '4-QAM  |  CMA', M, 'AWGN + PMD');

            testCase.verifyTrue(all(isfinite(eqSym(:))), ...
                'Equalizer output contains NaN/Inf.');
        end

        % -------- 16-QAM CMA+RDE ------------------
        function test16QAM_CMARERDE(testCase)
            M = 16;
            [rxSym, eqSym] = runScenario(testCase, M, 'CMA+RDE', ...
                testCase.NTaps, testCase.Mu, true, testCase.N1, 4000, ...
                testCase.NOut);
            plotBeforeAfter(testCase, rxSym, eqSym, ...
                '16-QAM  |  CMA+RDE', M, 'AWGN + PMD');

            testCase.verifyTrue(all(isfinite(eqSym(:))), ...
                'Equalizer output contains NaN/Inf.');
        end

        % -------- 16-QAM CMA+RDE + phase noise -----
        function test16QAM_CMARERDE_PhaseNoise(testCase)
            M = 16;
            [rxSym, eqSym] = runScenario(testCase, M, 'CMA+RDE', ...
                testCase.NTaps, testCase.Mu, true, testCase.N1, 4000, ...
                testCase.NOut, true);
            plotBeforeAfter(testCase, rxSym, eqSym, ...
                '16-QAM  |  CMA+RDE  |  PN', M, 'AWGN + PMD + Phase Noise');

            testCase.verifyTrue(all(isfinite(eqSym(:))), ...
                'Equalizer output contains NaN/Inf.');
        end

    end

    % ================================================================
    %  Helpers
    % ================================================================
    methods (Access = private)

        function [rxSym, eqSym] = runScenario(testCase, M, Eq, ...
                NTaps, Mu, SingleSpike, N1, N2, NOut, addPhaseNoise)
            if nargin < 10
                addPhaseNoise = false;
            end
            rng(42);

            % --- Tx ---
            k      = log2(M);
            Nbits  = k * testCase.N_pol * testCase.Ns;
            bits   = qam_randomBits(Nbits);
            symbols = qam_modulate(bits, M, testCase.N_pol);
            txSig  = qam_rectPulse(symbols, testCase.SpS);

            % --- Channel: AWGN + PMD (+ optional phase noise) ---
            rxSig = channel_add_awgn(txSig, testCase.SNR_dB);
            if addPhaseNoise
                linewidth = 100e4;   % 1 MHz
                rxSig = channel_add_phase_noise(rxSig, testCase.Rs, linewidth);
            end
            rxSig = channel_add_pmd(rxSig, testCase.L, testCase.SpS, ...
                testCase.Rs, testCase.DGDSpec, testCase.N_pmd);


            % --- Downsample before EQ for reference constellation ---
            rxSym = rxSig(1:testCase.SpS:end, :);

            % --- Adaptive Equalizer ---
            eqSig = adeq_equalize(rxSig, testCase.SpS, Eq, NTaps, Mu, ...
                SingleSpike, N1, N2, NOut);

            eqSym = eqSig;   % already at symbol rate after equalize
        end

        function plotBeforeAfter(testCase, rxSym, eqSym, titleStr, M, channelStr)
            if nargin < 6
                channelStr = 'AWGN + PMD';
            end
            figure('Name', titleStr, 'Position', [100 100 1200 500]);
            for p = 1:testCase.N_pol
                % Before
                subplot(2, 2, (p-1)*2 + 1);
                plot(real(rxSym(:,p)), imag(rxSym(:,p)), '.', ...
                     'MarkerSize', 2);
                grid on; axis equal;
                title(sprintf('Before Adaptive EQ  –  Pol %d', p));
                xlabel('In-Phase'); ylabel('Quadrature');

                % After
                subplot(2, 2, (p-1)*2 + 2);
                plot(real(eqSym(:,p)), imag(eqSym(:,p)), '.', ...
                     'MarkerSize', 2);
                grid on; axis equal;
                title(sprintf('After Adaptive EQ  –  Pol %d', p));
                xlabel('In-Phase'); ylabel('Quadrature');
            end
            sgtitle(sprintf('%d-QAM: %s  |  %s', M, channelStr, titleStr));
        end
    end
end
