classdef test_AdaptiveEqualizer < matlab.unittest.TestCase
    % Visual tests for AdaptiveEqualizer: AWGN+PMD then adaptive EQ.
    % Four scenarios are plotted:
    %   1. 4-QAM,  CMA only,       floating point
    %   2. 16-QAM, CMA+RDE,        floating point
    %   3. 4-QAM,  CMA only,       fixed point
    %   4. 16-QAM, CMA+RDE,        fixed point

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

        % Fixed-point settings
        WL      = 16
        FL      = 12
    end

    methods (TestMethodSetup)
        function seedRng(~)
            rng(42);
        end
    end

    % ================================================================
    methods (Test)

        % -------- 4-QAM CMA (floating point) -----------------------
        function testQPSK_CMA_Float(testCase)
            M = 4;
            paramDE.Eq          = 'CMA';
            paramDE.NTaps       = testCase.NTaps;
            paramDE.Mu          = testCase.Mu;
            paramDE.SingleSpike = true;
            paramDE.N1          = testCase.N1;
            paramDE.NOut        = testCase.NOut;

            [rxSym, eqSym] = runScenario(testCase, M, paramDE, false);
            plotBeforeAfter(testCase, rxSym, eqSym, ...
                '4-QAM  |  CMA  |  Float', M, 'AWGN + PMD');

            testCase.verifyTrue(all(isfinite(eqSym(:))), ...
                'Equalizer output contains NaN/Inf.');
        end

        % -------- 16-QAM CMA+RDE (floating point) ------------------
        function test16QAM_CMARERDE_Float(testCase)
            M = 16;
            paramDE.Eq          = 'CMA+RDE';
            paramDE.NTaps       = testCase.NTaps;
            paramDE.Mu          = testCase.Mu;
            paramDE.SingleSpike = true;
            paramDE.N1          = testCase.N1;
            paramDE.N2          = 4000;     % switch CMA->RDE
            paramDE.NOut        = testCase.NOut;

            [rxSym, eqSym] = runScenario(testCase, M, paramDE, false);
            plotBeforeAfter(testCase, rxSym, eqSym, ...
                '16-QAM  |  CMA+RDE  |  Float', M, 'AWGN + PMD');

            testCase.verifyTrue(all(isfinite(eqSym(:))), ...
                'Equalizer output contains NaN/Inf.');
        end

        % -------- 4-QAM CMA (fixed point) --------------------------
        function testQPSK_CMA_FixedPoint(testCase)
            M = 4;
            paramDE.Eq          = 'CMA';
            paramDE.NTaps       = testCase.NTaps;
            paramDE.Mu          = testCase.Mu;
            paramDE.SingleSpike = true;
            paramDE.N1          = testCase.N1;
            paramDE.NOut        = testCase.NOut;

            [rxSym, eqSym] = runScenario(testCase, M, paramDE, true);
            plotBeforeAfter(testCase, rxSym, eqSym, ...
                '4-QAM  |  CMA  |  Fixed Point', M, 'AWGN + PMD');

            testCase.verifyTrue(all(isfinite(eqSym(:))), ...
                'Equalizer output contains NaN/Inf.');
        end

        % -------- 16-QAM CMA+RDE (fixed point) ---------------------
        function test16QAM_CMARERDE_FixedPoint(testCase)
            M = 16;
            paramDE.Eq          = 'CMA+RDE';
            paramDE.NTaps       = testCase.NTaps;
            paramDE.Mu          = testCase.Mu;
            paramDE.SingleSpike = true;
            paramDE.N1          = testCase.N1;
            paramDE.N2          = 4000;
            paramDE.NOut        = testCase.NOut;

            [rxSym, eqSym] = runScenario(testCase, M, paramDE, true);
            plotBeforeAfter(testCase, rxSym, eqSym, ...
                '16-QAM  |  CMA+RDE  |  Fixed Point', M, 'AWGN + PMD');

            testCase.verifyTrue(all(isfinite(eqSym(:))), ...
                'Equalizer output contains NaN/Inf.');
        end

        % -------- 16-QAM CMA+RDE + phase noise (floating point) -----
        function test16QAM_CMARERDE_PhaseNoise_Float(testCase)
            M = 16;
            paramDE.Eq          = 'CMA+RDE';
            paramDE.NTaps       = testCase.NTaps;
            paramDE.Mu          = testCase.Mu;
            paramDE.SingleSpike = true;
            paramDE.N1          = testCase.N1;
            paramDE.N2          = 4000;
            paramDE.NOut        = testCase.NOut;

            [rxSym, eqSym] = runScenario(testCase, M, paramDE, false, true);
            plotBeforeAfter(testCase, rxSym, eqSym, ...
                '16-QAM  |  CMA+RDE  |  Float + PN', M, 'AWGN + PMD + Phase Noise');

            testCase.verifyTrue(all(isfinite(eqSym(:))), ...
                'Equalizer output contains NaN/Inf.');
        end

        % -------- 16-QAM CMA+RDE + phase noise (fixed point) --------
        function test16QAM_CMARERDE_PhaseNoise_FixedPoint(testCase)
            M = 16;
            paramDE.Eq          = 'CMA+RDE';
            paramDE.NTaps       = testCase.NTaps;
            paramDE.Mu          = testCase.Mu;
            paramDE.SingleSpike = true;
            paramDE.N1          = testCase.N1;
            paramDE.N2          = 4000;
            paramDE.NOut        = testCase.NOut;

            [rxSym, eqSym] = runScenario(testCase, M, paramDE, true, true);
            plotBeforeAfter(testCase, rxSym, eqSym, ...
                '16-QAM  |  CMA+RDE  |  FxP + PN', M, 'AWGN + PMD + Phase Noise');

            testCase.verifyTrue(all(isfinite(eqSym(:))), ...
                'Equalizer output contains NaN/Inf.');
        end
    end

    % ================================================================
    %  Helpers
    % ================================================================
    methods (Access = private)

        function [rxSym, eqSym] = runScenario(testCase, M, paramDE, useFixedPoint, addPhaseNoise)
            if nargin < 5
                addPhaseNoise = false;
            end
            rng(42);

            % --- Tx ---
            modem  = QAMModem(M, testCase.N_pol);
            k      = modem.bitsPerSymbol;
            Nbits  = k * testCase.N_pol * testCase.Ns;
            bits   = modem.randomBits(Nbits);
            symbols = modem.modulate(bits);
            txSig  = modem.rectPulse(symbols, testCase.SpS);

            % --- Channel: AWGN + PMD (+ optional phase noise) ---
            linewidth = 0;
            if addPhaseNoise
                linewidth = 100e3;   % 100 kHz
            end
            ch = Channel(testCase.L, testCase.SNR_dB, testCase.SpS, ...
                         testCase.Rs, testCase.D, testCase.CWL, ...
                         testCase.DGDSpec, testCase.N_pmd, linewidth);
            rxSig = ch.add_awgn(txSig);
            rxSig = ch.add_pmd(rxSig);
            if addPhaseNoise
                rxSig = ch.add_phase_noise(rxSig);
            end

            % --- Downsample before EQ for reference constellation ---
            rxSym = rxSig(1:testCase.SpS:end, :);

            % --- Adaptive Equalizer ---
            aeq   = AdaptiveEqualizer(paramDE, testCase.WL, testCase.FL);
            eqSig = aeq.equalize(rxSig, testCase.SpS, useFixedPoint);
            eqSig = double(eqSig);

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
