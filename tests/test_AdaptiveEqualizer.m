classdef test_AdaptiveEqualizer < matlab.unittest.TestCase
    % Visual tests for adaptive equalizer: AWGN+PMD then adaptive EQ.

    properties (Constant)
        N_pol   = 2
        Ns      = 2^17          % symbols per polarisation
        SpS     = 2

        % System
        Rs      = 32            % [GBd]
        L       = 80            % [km]
        SNR_dB  = 25            % [dB]
        D       = 0             % no CD for this test
        CWL     = 1550          % [nm]
        DGDSpec = 0.5           % PMD coeff [ps/sqrt(km)]
        N_pmd   = 5
        LW      = 0             % no phase noise

        % Pulse shaping
        Rolloff = 0.25
        Span    = 10

        % Adaptive EQ common settings
        NTaps   = 15
        Mu      = 1e-3
        N1      = 1000          % single-spike re-init iteration
        NOut    = 2000           % discard transient

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

        % -------- QPSK CMA -----------------------
        function testQPSK_CMA(testCase)
            [rxSym, eqSym, txBits, symbols] = runScenario(testCase, ...
                testCase.NTaps, testCase.Mu, true, testCase.N1, ...
                testCase.NOut);
            plotBeforeAfter(testCase, rxSym, eqSym, ...
                'QPSK  |  CMA', 'AWGN + PMD');

            testCase.verifyTrue(all(isfinite(eqSym(:))), ...
                'Equalizer output contains NaN/Inf.');

            % --- BER (resolve per-pol rotation + swap ambiguity) ---
            refSyms = symbols(testCase.NOut+1:end, :);
            % check refSyms and eqSym are the same size
            testCase.verifyEqual(size(refSyms), size(eqSym), ...
                'Reference symbols and equalizer output sizes differ.');
            totalErrors = 0;
            totalBits   = 0;
            for p = 1:testCase.N_pol
                refBitsPol = modem.symbolsToBits(refSyms(:,p));
                bestPolBER = Inf;
                % try both EQ outputs (equalizer may swap pols)
                for q = 1:testCase.N_pol
                    for kk = 0:31
                        rotated = eqSym(:,q) .* exp(-1j * kk * pi/16);
                        decSym  = modem.decideSymbols(rotated);
                        decBits = modem.symbolsToBits(decSym);
                        polBER  = sum(refBitsPol ~= decBits) / numel(refBitsPol);
                        if polBER < bestPolBER, bestPolBER = polBER; end
                    end
                end
                totalErrors = totalErrors + bestPolBER * numel(refBitsPol);
                totalBits   = totalBits + numel(refBitsPol);
            end
            bestBER = totalErrors / totalBits;
            fprintf('QPSK CMA BER = %.2e\n', bestBER);
            testCase.verifyLessThan(bestBER, testCase.BER_THRESHOLD, ...
                sprintf('QPSK CMA BER %.2e exceeds threshold.', bestBER));
        end

        % -------- QPSK CMA + phase noise -----
        function testQPSK_CMA_PhaseNoise(testCase)
            [rxSym, eqSym] = runScenario(testCase, ...
                testCase.NTaps, testCase.Mu, true, testCase.N1, ...
                testCase.NOut, true);
            plotBeforeAfter(testCase, rxSym, eqSym, ...
                'QPSK  |  CMA  |  PN', 'AWGN + PMD + Phase Noise');

            testCase.verifyTrue(all(isfinite(eqSym(:))), ...
                'Equalizer output contains NaN/Inf.');
        end

        % ============================================================
        %  Fixed-point tests: MATLAB fxp vs MEX fxp (same types)
        % ============================================================

        % % -------- QPSK CMA (fixed-point MATLAB vs MEX) ------
        % function testQPSK_CMA_Fxp(testCase)
        %     [rxSym, eqML, eqMEX] = runScenarioFxpBoth(testCase, ...
        %         testCase.NTaps, testCase.Mu, true, testCase.N1, ...
        %         testCase.NOut, 'fixed32');

        %     verifyFxpOutput(testCase, eqMEX, double(eqMEX), 'QPSK CMA FXP32-MEX');
        %     verifyFxpMexMatchesMatlab(testCase, eqML, eqMEX, 'QPSK CMA FXP32');

        %     plotFxpComparison(testCase, rxSym, eqML, eqMEX, ...
        %         'QPSK  |  CMA  |  FXP32');
        % end

        % % -------- QPSK CMA + phase noise (fixed-point MATLAB vs MEX) ---
        % function testQPSK_CMA_PhaseNoise_Fxp(testCase)
        %     [rxSym, eqML, eqMEX] = runScenarioFxpBoth(testCase, ...
        %         testCase.NTaps, testCase.Mu, true, testCase.N1, ...
        %         testCase.NOut, 'fixed32', true);

        %     verifyFxpOutput(testCase, eqMEX, double(eqMEX), 'QPSK CMA PN FXP32-MEX');
        %     verifyFxpMexMatchesMatlab(testCase, eqML, eqMEX, 'QPSK CMA PN FXP32');

        %     plotFxpComparison(testCase, rxSym, eqML, eqMEX, ...
        %         'QPSK  |  CMA  |  PN  |  FXP32');
        % end

    end

    % ================================================================
    %  Helpers
    % ================================================================
    methods (Access = private)

        function [rxSym, eqSym, bits, symbols] = runScenario(testCase, ...
                NTaps, Mu, SingleSpike, N1, NOut, addPhaseNoise, SignOnly)
            if nargin < 7
                addPhaseNoise = false;
            end
            if nargin < 8
                SignOnly = false;
            end
            rng(42);

            % --- Tx ---
            Nbits  = 4 * testCase.Ns;
            bits   = modem.randomBits(Nbits);
            symbols = modem.modulate(bits);
            % txSig  = modem.rrcPulse(symbols, testCase.SpS, ...
            %     testCase.Rolloff, testCase.Span);
            % duplicate for SpS > 1
            txSig = repelem(symbols, testCase.SpS, 1);

            % --- Channel: AWGN + PMD (+ optional phase noise) ---
            rxSig = channel.add_awgn(txSig, testCase.SNR_dB);
            if addPhaseNoise
                linewidth = 100e4;   % 1 MHz
                rxSig = channel.add_phase_noise(rxSig, testCase.Rs, linewidth);
            end
            rxSig = channel.add_pmd(rxSig, testCase.L, testCase.SpS, ...
                testCase.Rs, testCase.DGDSpec, testCase.N_pmd);

            % --- Matched filter ---
            % rxSig = modem.matched_filter(rxSig, testCase.SpS, 'rrc', ...
            %     testCase.Rolloff, testCase.Span);

            rxSym = rxSig(1:testCase.SpS:end, :);

            % --- Adaptive Equalizer ---
            eqSig = adaptive_eq.equalize(rxSig, testCase.SpS, NTaps, Mu, ...
                SingleSpike, N1, NOut, SignOnly);

            eqSym = eqSig;   % already at symbol rate after equalize
        end

        function [rxSym, eqML, eqMEX] = runScenarioFxpBoth(testCase, ...
                NTaps, Mu, SingleSpike, N1, NOut, fxpConfig, addPhaseNoise, SignOnly)
            %RUNSCENARIOFXPBOTH  Run both MATLAB fxp and MEX fxp on
            %  identical fi input, returning both outputs for comparison.
            if nargin < 8
                addPhaseNoise = false;
            end
            if nargin < 9
                SignOnly = false;
            end
            rng(42);

            % --- Tx ---
            Nbits  = 4 * testCase.Ns;
            bits   = modem.randomBits(Nbits);
            symbols = modem.modulate(bits);
            txSig  = modem.rrcPulse(symbols, testCase.SpS, ...
                testCase.Rolloff, testCase.Span);

            % --- Channel: AWGN + PMD (+ optional phase noise) ---
            rxSig = channel.add_awgn(txSig, testCase.SNR_dB);
            if addPhaseNoise
                linewidth = 100e4;   % 1 MHz
                rxSig = channel.add_phase_noise(rxSig, testCase.Rs, linewidth);
            end
            rxSig = channel.add_pmd(rxSig, testCase.L, testCase.SpS, ...
                testCase.Rs, testCase.DGDSpec, testCase.N_pmd);

            % --- Matched filter ---
            rxSig = modem.matched_filter(rxSig, testCase.SpS, 'rrc', ...
                testCase.Rolloff, testCase.Span);

            % --- Downsample before EQ for reference constellation ---
            rxSym = rxSig(1:testCase.SpS:end, :);

            % --- Cast input to fixed-point ---
            T = adaptive_eq.equalize_fxp_types(fxpConfig);
            rxSig_fi = cast(rxSig, 'like', T.x);

            % --- Run MATLAB fixed-point function ---
            eqML = adaptive_eq.equalize_fxp(rxSig_fi, ...
                double(testCase.SpS), double(NTaps), double(Mu), ...
                SingleSpike, double(N1), double(NOut), logical(SignOnly), T);

            % --- Run fixed-point MEX ---
            eqMEX = adaptive_eq.equalize_fxp_mex(rxSig_fi, ...
                double(testCase.SpS), double(NTaps), double(Mu), ...
                SingleSpike, double(N1), double(NOut), logical(SignOnly), T);
        end

        function verifyFxpOutput(testCase, eqSym, eqSym_dbl, tag)
            % --- Shape ---
            testCase.verifyEqual(size(eqSym, 2), 2, ...
                sprintf('%s: output must have 2 columns (polarisations).', tag));
            testCase.verifyGreaterThan(size(eqSym, 1), 0, ...
                sprintf('%s: output must have at least 1 row.', tag));

            % --- Type ---
            testCase.verifyTrue(isa(eqSym, 'embedded.fi'), ...
                sprintf('%s: output must be fi (embedded.fi) type.', tag));

            % --- Complexity ---
            testCase.verifyTrue(~isreal(eqSym_dbl), ...
                sprintf('%s: output must be complex-valued.', tag));

            % --- Finiteness ---
            testCase.verifyTrue(all(isfinite(eqSym_dbl(:))), ...
                sprintf('%s: output contains NaN/Inf.', tag));

            % --- Non-trivial output (not all zeros) ---
            testCase.verifyGreaterThan(max(abs(eqSym_dbl(:))), 0, ...
                sprintf('%s: output is all zeros.', tag));

            % --- Reasonable amplitude (signal normalised ~1) ---
            maxAmp = max(abs(eqSym_dbl(:)));
            testCase.verifyLessThan(maxAmp, 10, ...
                sprintf('%s: output amplitude suspiciously large (%.2f).', tag, maxAmp));
        end

        function verifyFxpMexMatchesMatlab(testCase, eqML, eqMEX, tag)
            % Verify MEX output is bit-exact with MATLAB fxp output
            testCase.verifyEqual(size(eqMEX), size(eqML), ...
                sprintf('%s: MEX and MATLAB output sizes differ.', tag));

            mlDbl  = double(eqML);
            mexDbl = double(eqMEX);

            testCase.verifyEqual(mexDbl, mlDbl, ...
                sprintf('%s: MEX output differs from MATLAB fxp output (should be bit-exact).', tag));

            % Also verify matching fi properties
            if isa(eqML, 'embedded.fi') && isa(eqMEX, 'embedded.fi')
                testCase.verifyEqual(eqMEX.WordLength, eqML.WordLength, ...
                    sprintf('%s: MEX WordLength differs from MATLAB.', tag));
                testCase.verifyEqual(eqMEX.FractionLength, eqML.FractionLength, ...
                    sprintf('%s: MEX FractionLength differs from MATLAB.', tag));
            end
        end

        function plotFxpComparison(testCase, rxSym, eqML, eqMEX, titleStr)
            mlDbl  = double(eqML);
            mexDbl = double(eqMEX);

            figure('Name', titleStr, 'Position', [100 100 1400 700]);
            for p = 1:testCase.N_pol
                % Before EQ
                subplot(2, 3, (p-1)*3 + 1);
                plot(real(rxSym(:,p)), imag(rxSym(:,p)), '.', 'MarkerSize', 2);
                grid on; axis equal;
                title(sprintf('Before EQ – Pol %d', p));
                xlabel('I'); ylabel('Q');

                % MATLAB fxp
                subplot(2, 3, (p-1)*3 + 2);
                plot(real(mlDbl(:,p)), imag(mlDbl(:,p)), '.', 'MarkerSize', 2);
                grid on; axis equal;
                title(sprintf('MATLAB FXP – Pol %d', p));
                xlabel('I'); ylabel('Q');

                % MEX fxp
                subplot(2, 3, (p-1)*3 + 3);
                plot(real(mexDbl(:,p)), imag(mexDbl(:,p)), '.', 'MarkerSize', 2);
                grid on; axis equal;
                title(sprintf('MEX FXP – Pol %d', p));
                xlabel('I'); ylabel('Q');
            end
            sgtitle(sprintf('QPSK  |  %s  |  MATLAB vs MEX', titleStr));
        end

        function plotBeforeAfter(testCase, rxSym, eqSym, titleStr, channelStr)
            if nargin < 5
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
            sgtitle(sprintf('QPSK: %s  |  %s', channelStr, titleStr));
        end
    end
end
