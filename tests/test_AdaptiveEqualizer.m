classdef test_AdaptiveEqualizer < matlab.unittest.TestCase
    % Visual test for the parallel-lane adaptive equalizer:
    % AWGN + PMD channel, then parallel CMA, then best BER over a
    % phase-rotation sweep with polarisation-ambiguity resolution.

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
        NOut    = 2000          % discard transient
        PLanes  = 32             % number of parallel lanes

        % Fixed-point
        FxpConfig  = 'fixed32'   % equalize_fxp_types configuration
        UpdateStep = 1           % output samples between weight updates

        % Pass / fail
        BER_THRESHOLD = 5e-2
    end

    methods (TestClassSetup)
        function buildFxpMex(testCase)
            % Build the fixed-point MEX once before the fxp test runs so
            % it picks up the current equalize_fxp (incl. PLanes arg).
            thisDir  = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(thisDir);
            addpath(genpath(fullfile(repoRoot, 'src')));
            addpath(fullfile(repoRoot, 'build'));

            P = struct();
            P.FxpConfig_AEQ   = testCase.FxpConfig;
            P.SpS             = testCase.SpS;
            P.AEQ_NTaps       = testCase.NTaps;
            P.AEQ_Mu          = testCase.Mu;
            P.AEQ_SingleSpike = true;
            P.AEQ_N1          = testCase.N1;
            P.AEQ_NOut        = testCase.NOut;
            P.AEQ_SignOnly    = false;
            P.AEQ_UpdateStep  = testCase.UpdateStep;
            P.AEQ_PLanes      = testCase.PLanes;

            cfg = coder.config('mex');
            build_adaptive_eq_equalize_fxp_mex(P, cfg);
        end
    end

    methods (TestMethodSetup)
        function seedRng(~)
            rng(42);
        end
    end

    % ================================================================
    methods (Test)

        % -------- QPSK parallel CMA: AWGN + PMD -----------------------
        function testQPSK_ParallelCMA(testCase)
            [rxSym, eqSym, symbols] = runScenario(testCase, ...
                testCase.NTaps, testCase.Mu, true, testCase.N1, ...
                testCase.NOut, testCase.PLanes);

            plotBeforeAfter(testCase, rxSym, eqSym, ...
                sprintf('QPSK  |  Parallel CMA  |  %d lanes', testCase.PLanes), ...
                'AWGN + PMD');

            testCase.verifyTrue(all(isfinite(eqSym(:))), ...
                'Equalizer output contains NaN/Inf.');

            % --- BER (resolve per-pol rotation + swap ambiguity) ---
            refSyms = symbols(testCase.NOut+1:end, :);
            testCase.verifyEqual(size(refSyms), size(eqSym), ...
                'Reference symbols and equalizer output sizes differ.');

            bestBER = computeBestBER(testCase, refSyms, eqSym);
            fprintf('QPSK parallel CMA (%d lanes) BER = %.2e\n', ...
                testCase.PLanes, bestBER);
            testCase.verifyLessThan(bestBER, testCase.BER_THRESHOLD, ...
                sprintf('QPSK parallel CMA BER %.2e exceeds threshold.', ...
                bestBER));
        end

        % -------- QPSK parallel CMA (fixed-point): AWGN + PMD ---------
        function testQPSK_ParallelCMA_Fxp(testCase)
            [rxSym, eqSym, symbols] = runScenarioFxp(testCase, ...
                testCase.NTaps, testCase.Mu, true, testCase.N1, ...
                testCase.NOut, testCase.PLanes, testCase.FxpConfig);

            eqSymD = double(eqSym);
            plotBeforeAfter(testCase, rxSym, eqSymD, ...
                sprintf('QPSK  |  Parallel CMA  |  %s  |  %d lanes', ...
                testCase.FxpConfig, testCase.PLanes), 'AWGN + PMD');

            testCase.verifyTrue(isa(eqSym, 'embedded.fi'), ...
                'Fixed-point equalizer output must be a fi object.');
            testCase.verifyTrue(all(isfinite(eqSymD(:))), ...
                'Fixed-point equalizer output contains NaN/Inf.');

            % --- BER (resolve per-pol rotation + swap ambiguity) ---
            refSyms = symbols(testCase.NOut+1:end, :);
            testCase.verifyEqual(size(refSyms), size(eqSymD), ...
                'Reference symbols and equalizer output sizes differ.');

            bestBER = computeBestBER(testCase, refSyms, eqSymD);
            fprintf('QPSK parallel CMA FXP %s (%d lanes) BER = %.2e\n', ...
                testCase.FxpConfig, testCase.PLanes, bestBER);
            testCase.verifyLessThan(bestBER, testCase.BER_THRESHOLD, ...
                sprintf('QPSK parallel CMA FXP BER %.2e exceeds threshold.', ...
                bestBER));
        end

    end

    % ================================================================
    %  Helpers
    % ================================================================
    methods (Access = private)

        function [rxSig, rxSym, symbols] = genRx(testCase)
            % Shared Tx + channel (AWGN + PMD) generation.
            rng(42);

            % --- Tx ---
            % Plain QPSK symbols at +/-1 +/-1j (no pilots/training, which
            % would disrupt the blind equaliser).
            symbols = (2*randi([0 1], testCase.Ns, testCase.N_pol) - 1) ...
                + 1j*(2*randi([0 1], testCase.Ns, testCase.N_pol) - 1);
            % duplicate for SpS > 1
            txSig = repelem(symbols, testCase.SpS, 1);

            % --- Channel: AWGN + PMD ---
            rxSig = channel.add_awgn(txSig, testCase.SNR_dB);
            rxSig = channel.add_pmd(rxSig, testCase.L, testCase.SpS, ...
                testCase.Rs, testCase.DGDSpec, testCase.N_pmd);

            rxSym = rxSig(1:testCase.SpS:end, :);
        end

        function [rxSym, eqSym, symbols] = runScenario(testCase, ...
                NTaps, Mu, SingleSpike, N1, NOut, PLanes, SignOnly)
            if nargin < 8
                SignOnly = false;
            end

            [rxSig, rxSym, symbols] = genRx(testCase);

            % --- Parallel adaptive equalizer (floating-point) ---
            eqSig = adaptive_eq.equalize(rxSig, testCase.SpS, NTaps, Mu, ...
                SingleSpike, N1, NOut, SignOnly, PLanes);

            eqSym = eqSig;   % already at symbol rate after equalize
        end

        function [rxSym, eqSym, symbols] = runScenarioFxp(testCase, ...
                NTaps, Mu, SingleSpike, N1, NOut, PLanes, FxpConfig, SignOnly)
            if nargin < 9
                SignOnly = false;
            end

            [rxSig, rxSym, symbols] = genRx(testCase);

            % --- Cast input to fixed-point ---
            T = adaptive_eq.equalize_fxp_types(FxpConfig);
            rxSig_fi = cast(rxSig, 'like', T.x);

            % --- Parallel adaptive equalizer (fixed-point MEX) ---
            eqSym = adaptive_eq.equalize_fxp_mex(rxSig_fi, ...
                double(testCase.SpS), double(NTaps), double(Mu), ...
                SingleSpike, double(N1), double(NOut), logical(SignOnly), ...
                double(testCase.UpdateStep), T, double(PLanes));
        end

        function bestBER = computeBestBER(testCase, refSyms, eqSym)
            % Resolve per-pol phase rotation (multiples of pi/16) and
            % the polarisation swap ambiguity, return the lowest BER.
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
        end

        function plotBeforeAfter(testCase, rxSym, eqSym, titleStr, channelStr)
            if nargin < 5
                channelStr = 'AWGN + PMD';
            end
            % Discard the same transient samples as the equalizer so
            % only the kept (non-discarded) points are plotted.
            rxSym = rxSym(testCase.NOut+1:end, :);

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
