classdef test_ClockRecovery < matlab.unittest.TestCase
    %TEST_CLOCKRECOVERY  Verify clk_recovery.recovery_godard corrects
    %  constant timing offsets and sampling-frequency offsets on
    %  RRC-shaped signals (feedback Modified Godard).
    %
    %  Approach: symbols are pulse-shaped at a high oversampling factor
    %  (SpS_hi = 16) so that arbitrary sub-sample timing shifts can be
    %  applied by selecting different sample phases.  The signal is then
    %  decimated to 2 Sa/symbol before being fed to recovery_godard.
    %
    %  Each scenario performs a 2D sweep over the PI loop-filter gains
    %  (ki, kp) and reports the combination with the lowest per-quadrant
    %  cluster variance.  This is a phase-ambiguity-invariant measure of
    %  constellation tightness and is sensitive to timing residual even
    %  in regimes where BER has saturated near zero (or near random).

    properties (Constant)
        N_pol   = 1              % single polarisation (per-pol operation)
        Ns      = 2^18           % symbols
        SpS     = 2              % target samples per symbol
        SpS_hi  = 16             % high-resolution oversampling

        % RRC parameters
        Rolloff = 0.25
        Span    = 10

        % Modified Godard estimator parameters
        N_fft = 128

        % DPLL clock recovery parameters
        NLanes_DPLL = 2          % parallel lanes per block in recovery()

        % PI loop-filter gain sweep ranges (log-spaced)
        ki_sweep = logspace(-7, -3, 9)
        kp_sweep = logspace(-6, -2, 9)

        % Pass/fail (BER computed at the optimum for reporting only)
        BER_THRESHOLD = 1e-2
    end

    methods (TestMethodSetup)
        function seedRng(~)
            rng(42);
        end
    end

    % ================================================================
    methods (Test)

        % -------- Modified Godard: constant timing -------------------
        function testGodard_ConstantTimingOffset(testCase)
            runGodardScenario(testCase, 'constant');
        end

        % -------- Modified Godard: SFO -------------------------------
        function testGodard_SamplingFrequencyOffset(testCase)
            runGodardScenario(testCase, 'sfo');
        end

        % -------- DPLL: constant timing ------------------------------
        function testDPLL_ConstantTimingOffset(testCase)
            runDPLLScenario(testCase, 'constant');
        end

        % -------- DPLL: SFO ------------------------------------------
        function testDPLL_SamplingFrequencyOffset(testCase)
            runDPLLScenario(testCase, 'sfo');
        end

    end

    % ================================================================
    %  Helpers
    % ================================================================
    methods (Access = private)

        function runGodardScenario(testCase, impairment)
            SpS_hi_ = testCase.SpS_hi;
            SpS_    = testCase.SpS;
            Ns_     = testCase.Ns;

            % --- Tx: generate symbols & pulse-shape at high SpS ---
            Nbits  = 2 * Ns_;              % QPSK: 2 bits/sym, single pol
            txBits = modem.randomBits(Nbits);
            symbols = modem.modulate(txBits);
            symbols = symbols(:, 1);        % single pol

            txHi = modem.rrcPulse(symbols, SpS_hi_, testCase.Rolloff, testCase.Span);

            % --- Matched filter at high SpS ---
            rxHi = modem.matched_filter(txHi, SpS_hi_, 'rrc', ...
                testCase.Rolloff, testCase.Span);

            % --- Impairment (via channel.apply_timing_error) ---
            switch impairment
                case 'constant'
                    tau0  = 3 / SpS_hi_;    % 3/16 of a symbol period
                    rxImp = channel.apply_timing_error(rxHi, 0, tau0, SpS_hi_);
                    impLabel = sprintf('Constant Offset tau0 = %.3f T', tau0);
                case 'sfo'
                    ppm   = 40;
                    rxImp = channel.apply_timing_error(rxHi, ppm, 0, SpS_hi_);
                    impLabel = sprintf('SFO %d ppm', ppm);
                otherwise
                    error('Unknown impairment: %s', impairment);
            end

            % --- Decimate to 2 Sa/symbol ---
            decFactor = SpS_hi_ / SpS_;
            rx2 = rxImp(1:decFactor:end, :);

            txRefBits = modem.symbolsToBits(symbols);
            skipSym   = ceil(testCase.N_fft / (2 * SpS_));
            skipBits  = skipSym * 2;        % QPSK: 2 bits per symbol

            % --- 2D sweep over (ki, kp) ---
            ki_vec  = testCase.ki_sweep;
            kp_vec  = testCase.kp_sweep;
            var_map = nan(numel(ki_vec), numel(kp_vec));

            bestVar   = Inf;
            best_ki   = NaN;
            best_kp   = NaN;
            bestCrSym = [];

            for ii = 1:numel(ki_vec)
                for jj = 1:numel(kp_vec)
                    ki = ki_vec(ii);
                    kp = kp_vec(jj);

                    crOut = clk_recovery.recovery_godard(rx2, Ns_, ...
                        testCase.N_fft, testCase.Rolloff, ki, kp);

                    crSym  = crOut(1:SpS_:end);
                    crEval = crSym(skipSym+1:end-skipSym);

                    if ~all(isfinite(crEval))
                        var_map(ii, jj) = NaN;
                        continue;
                    end

                    % Normalise amplitude so the metric is scale-free
                    crEvalN = crEval / sqrt(mean(abs(crEval).^2));

                    qVar = quadrantVariance(crEvalN);
                    var_map(ii, jj) = qVar;

                    if qVar < bestVar
                        bestVar   = qVar;
                        best_ki   = ki;
                        best_kp   = kp;
                        bestCrSym = crEvalN;
                    end
                end
            end

            % --- BER at the optimum (reporting only) ---
            bestBER = bestRotationBER(testCase, bestCrSym, ...
                txRefBits(skipBits+1:end-skipBits));

            fprintf(['Godard %s  min quad-var = %.3e  (BER = %.2e)  ' ...
                'at ki = %.3g, kp = %.3g\n'], impLabel, bestVar, ...
                bestBER, best_ki, best_kp);

            % --- Plot: variance heat-map + best-case constellations ---
            plotSweep(testCase, ki_vec, kp_vec, var_map, ...
                best_ki, best_kp, impLabel);

            rxSym = rx2(1:SpS_:end);
            plotBeforeAfter(testCase, rxSym, bestCrSym, ...
                sprintf('Godard  |  %s  |  best ki=%.2g kp=%.2g', ...
                impLabel, best_ki, best_kp), bestVar);

            % --- Verify ---
            testCase.verifyTrue(isfinite(bestVar), ...
                'No finite quadrant-variance configuration found.');
            testCase.verifyLessThan(bestBER, testCase.BER_THRESHOLD, ...
                sprintf(['Godard %s best BER %.2e (ki=%.3g, kp=%.3g) ' ...
                'exceeds threshold.'], impLabel, bestBER, best_ki, best_kp));
        end

        function runDPLLScenario(testCase, impairment)
            SpS_hi_ = testCase.SpS_hi;
            SpS_    = testCase.SpS;
            Ns_     = testCase.Ns;
            NLanes_ = testCase.NLanes_DPLL;

            % --- Tx: generate symbols & pulse-shape at high SpS ---
            Nbits  = 2 * Ns_;
            txBits = modem.randomBits(Nbits);
            symbols = modem.modulate(txBits);
            symbols = symbols(:, 1);

            txHi = modem.rrcPulse(symbols, SpS_hi_, testCase.Rolloff, testCase.Span);

            rxHi = modem.matched_filter(txHi, SpS_hi_, 'rrc', ...
                testCase.Rolloff, testCase.Span);

            switch impairment
                case 'constant'
                    tau0  = 3 / SpS_hi_;
                    rxImp = channel.apply_timing_error(rxHi, 0, tau0, SpS_hi_);
                    impLabel = sprintf('Constant Offset tau0 = %.3f T', tau0);
                case 'sfo'
                    ppm   = 40;
                    rxImp = channel.apply_timing_error(rxHi, ppm, 0, SpS_hi_);
                    impLabel = sprintf('SFO %d ppm', ppm);
                otherwise
                    error('Unknown impairment: %s', impairment);
            end

            decFactor = SpS_hi_ / SpS_;
            rx2 = rxImp(1:decFactor:end, :);

            txRefBits = modem.symbolsToBits(symbols);
            skipSym   = max(ceil(Ns_ * 0.1), 100);  % PLL transient
            skipBits  = skipSym * 2;

            % --- 2D sweep over (ki, kp) ---
            ki_vec  = testCase.ki_sweep;
            kp_vec  = testCase.kp_sweep;
            var_map = nan(numel(ki_vec), numel(kp_vec));

            bestVar   = Inf;
            best_ki   = NaN;
            best_kp   = NaN;
            bestCrSym = [];

            for ii = 1:numel(ki_vec)
                for jj = 1:numel(kp_vec)
                    ki = ki_vec(ii);
                    kp = kp_vec(jj);

                    crOut = clk_recovery.recovery(rx2, 'Nyquist', Ns_, ...
                        ki, kp, NLanes_);

                    % Try both symbol phases; keep the tighter one
                    crSymA = crOut(1:SpS_:end);
                    crSymB = crOut(2:SpS_:end);

                    [qVar, crEvalN] = pickBestPhase(crSymA, crSymB, skipSym);
                    var_map(ii, jj) = qVar;

                    if isfinite(qVar) && qVar < bestVar
                        bestVar   = qVar;
                        best_ki   = ki;
                        best_kp   = kp;
                        bestCrSym = crEvalN;
                    end
                end
            end

            bestBER = bestRotationBER(testCase, bestCrSym, ...
                txRefBits(skipBits+1:end-skipBits));

            fprintf(['DPLL %s  min quad-var = %.3e  (BER = %.2e)  ' ...
                'at ki = %.3g, kp = %.3g  (NLanes = %d)\n'], impLabel, ...
                bestVar, bestBER, best_ki, best_kp, NLanes_);

            plotSweep(testCase, ki_vec, kp_vec, var_map, ...
                best_ki, best_kp, sprintf('DPLL | %s', impLabel));

            rxSym = rx2(1:SpS_:end);
            plotBeforeAfter(testCase, rxSym, bestCrSym, ...
                sprintf('DPLL  |  %s  |  best ki=%.2g kp=%.2g', ...
                impLabel, best_ki, best_kp), bestVar);

            testCase.verifyTrue(isfinite(bestVar), ...
                'No finite quadrant-variance configuration found.');
            testCase.verifyLessThan(bestBER, testCase.BER_THRESHOLD, ...
                sprintf(['DPLL %s best BER %.2e (ki=%.3g, kp=%.3g) ' ...
                'exceeds threshold.'], impLabel, bestBER, best_ki, best_kp));
        end

        function BER = bestRotationBER(~, crSym, txRefBits)
            %BESTROTATIONBER  Try all four pi/2 rotations, return lowest BER.
            BER = Inf;
            for kk = 0:3
                rotated = crSym .* exp(-1j * kk * pi/2);
                decided = modem.decideSymbols(rotated);
                rxBits  = modem.symbolsToBits(decided);
                nBits   = min(length(rxBits), length(txRefBits));
                nErr    = sum(txRefBits(1:nBits) ~= rxBits(1:nBits));
                thisBER = nErr / nBits;
                if thisBER < BER
                    BER = thisBER;
                end
            end
        end

        function plotSweep(~, ki_vec, kp_vec, var_map, best_ki, best_kp, impLabel)
            figure('Name', sprintf('Godard PI sweep | %s', impLabel), ...
                'Position', [100 100 700 500]);

            logVar = log10(max(var_map, 1e-6));
            imagesc(log10(kp_vec), log10(ki_vec), logVar);
            set(gca, 'YDir', 'normal');
            xlabel('log_{10}(k_p)');
            ylabel('log_{10}(k_i)');
            cb = colorbar;
            ylabel(cb, 'log_{10}(mean quadrant variance)');
            hold on;
            plot(log10(best_kp), log10(best_ki), 'rx', ...
                'MarkerSize', 12, 'LineWidth', 2);
            hold off;
            title(sprintf('PI gain sweep  |  %s', impLabel));
        end

        function plotBeforeAfter(~, rxSym, crSym, titleStr, qVar)
            figure('Name', titleStr, 'Position', [100 100 900 400]);

            subplot(1, 2, 1);
            plot(real(rxSym), imag(rxSym), '.', 'MarkerSize', 2);
            grid on; axis equal;
            title('Before Clock Recovery');
            xlabel('I'); ylabel('Q');

            subplot(1, 2, 2);
            plot(real(crSym), imag(crSym), '.', 'MarkerSize', 2);
            grid on; axis equal;
            title('After Clock Recovery');
            xlabel('I'); ylabel('Q');

            sgtitle(sprintf('QPSK  |  %s  |  quad-var = %.3e', titleStr, qVar));
        end
    end
end


% --------------------------------------------------------------------
function v = quadrantVariance(sym)
%QUADRANTVARIANCE  Mean cluster variance about per-quadrant centroids.
%   Bins each sample by the sign of its real/imag parts (the four QPSK
%   quadrants), computes the centroid of each non-empty bin, and returns
%   the mean of E[|x - centroid|^2] across the bins.  This is invariant
%   to ±π/2 phase ambiguity and grows monotonically as timing residual
%   smears the clusters.

    sym = sym(:);
    qIdx = (real(sym) > 0) + 2 * (imag(sym) > 0);   % 0..3

    accum = 0;
    cnt   = 0;
    for q = 0:3
        mask = (qIdx == q);
        n    = sum(mask);
        if n < 2
            continue;
        end
        cluster = sym(mask);
        centroid = mean(cluster);
        accum    = accum + sum(abs(cluster - centroid).^2);
        cnt      = cnt + n;
    end

    if cnt == 0
        v = NaN;
    else
        v = accum / cnt;
    end
end


function [v, crEvalN] = pickBestPhase(crSymA, crSymB, skipSym)
%PICKBESTPHASE  Evaluate quadrant variance on both candidate symbol
%   phases (Sa/symbol = 2 -> two possible strobe phases) and return the
%   one with the tighter clusters together with its normalised samples.

    function [v_, x_] = evalPhase(crSym)
        if length(crSym) <= 2*skipSym
            v_ = NaN; x_ = [];
            return;
        end
        crEval = crSym(skipSym+1:end-skipSym);
        if ~all(isfinite(crEval))
            v_ = NaN; x_ = [];
            return;
        end
        x_ = crEval / sqrt(mean(abs(crEval).^2));
        v_ = quadrantVariance(x_);
    end

    [vA, xA] = evalPhase(crSymA);
    [vB, xB] = evalPhase(crSymB);

    if ~isfinite(vA) && ~isfinite(vB)
        v = NaN; crEvalN = [];
    elseif ~isfinite(vB) || (isfinite(vA) && vA <= vB)
        v = vA; crEvalN = xA;
    else
        v = vB; crEvalN = xB;
    end
end
