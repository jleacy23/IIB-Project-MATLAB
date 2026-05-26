classdef test_EqClkCombined < matlab.unittest.TestCase
    %TEST_EQCLKCOMBINED  Verify the four combined equalisation +
    %   clock-recovery blocks in +eq_clk on an RRC DP-QPSK signal
    %   impaired with a small amount of chromatic dispersion, PMD,
    %   constant timing offset and AWGN.  Each block:
    %     1. combined_adaptive_gardner          - CMA + Gardner DPLL
    %     2. combined_cd_td_gardner_adaptive    - CD/MF FIR + Gardner + CMA
    %     3. combined_cd_fd_godard_adaptive     - CD/MF (overlap-save) +
    %                                             Godard (shared FFT) + CMA
    %     4. combined_cd_fd_gardner_adaptive    - CD/MF (overlap-save) +
    %                                             Gardner + CMA
    %   is exercised in its own test, plotted before/after, and the
    %   recovered BER (per-pol phase + pol-swap + small lag search) is
    %   verified against a loose threshold.

    properties (Constant)
        N_pol  = 2
        Ns     = 2^14         % symbols per polarisation
        SpS    = 2

        % --- System ---
        Rs      = 32          % [GBd]
        L_km    = 80          % small CD
        D       = 17          % [ps/(nm*km)]
        CWL     = 1550        % [nm]
        DGDSpec = 0.1        % small PMD [ps/sqrt(km)]
        N_pmd   = 1
        SNR_dB  = 25
        tau0    = 0       % constant timing offset [symbol periods]
        SFO_ppm = 40           % sample-frequency offset [ppm] (0 = off)

        % --- Pulse shaping ---
        Rolloff = 0.25
        Span    = 10

        % --- Adaptive EQ ---
        NTapsAdapt     = 3
        Mu             = 1e-3
        SingleSpike    = true
        N1             = 100       % iteration to reinit y-pol weights
        NOutAdapt      = 200       % samples to discard after equalisation
        SignOnly       = true
        PLanes         = 32       % adaptive_eq parallel lanes
        Mode           = 0       % 0 = CMA, 1 = pilot-aided
        BlockLen       = []      % [] -> defaults to PLanes
        SubframeBlocks = 0       % 0 disables CMA subframe skip
        NLanes         = 32       % clk_recovery parallel lanes

        % --- CD time-domain FIR length (block 2) ---
        NTapCD     = 31

        % --- CD frequency-domain overlap-save (blocks 3, 4) ---
        NFFT       = 256
        NOverlap   = 64

        % --- DPLL / Godard loop-filter gain sweep (log-spaced) ---
        ki_sweep_gardner = logspace(-7, -4, 4)
        kp_sweep_gardner = logspace(-6, -3, 4)
        ki_sweep_godard  = logspace(-6, -3, 4)
        kp_sweep_godard  = logspace(-5, -2, 4)

        % --- Pass / fail ---
        BER_THRESHOLD = 5e-2
    end

    methods (TestMethodSetup)
        function seedRng(~)
            rng(42);
        end
    end

    % ================================================================
    methods (Test)

        function test_block1_adaptive_gardner(testCase)
            label = 'Block 1: Gardner -> Adaptive CMA';
            [rxSig, symbols] = genRx(testCase);
            runFn = @(ki, kp) eq_clk.combined_adaptive_gardner( ...
                rxSig, testCase.SpS, ki, kp, ...
                testCase.Ns, testCase.NLanes, true, adaptOpts(testCase));
            [y, ki_b, kp_b, qVar_b] = runSweep(testCase, runFn, ...
                testCase.ki_sweep_gardner, testCase.kp_sweep_gardner, label);
            evaluate(testCase, rxSig, y, symbols, ...
                sprintf('%s  (best ki=%.2g kp=%.2g qVar=%.2e)', ...
                label, ki_b, kp_b, qVar_b));
        end

        function test_block2_cd_td_gardner_adaptive(testCase)
            label = 'Block 2: CD/MF FIR -> Gardner -> Adaptive CMA';
            [rxSig, symbols] = genRx(testCase);
            runFn = @(ki, kp) eq_clk.combined_cd_td_gardner_adaptive( ...
                rxSig, testCase.SpS, testCase.NTapCD, ...
                testCase.D, testCase.L_km, testCase.CWL, testCase.Rs, ...
                testCase.Rolloff, testCase.Span, ki, kp, ...
                testCase.Ns, testCase.NLanes, adaptOpts(testCase));
            [y, ki_b, kp_b, qVar_b] = runSweep(testCase, runFn, ...
                testCase.ki_sweep_gardner, testCase.kp_sweep_gardner, label);
            evaluate(testCase, rxSig, y, symbols, ...
                sprintf('%s  (best ki=%.2g kp=%.2g qVar=%.2e)', ...
                label, ki_b, kp_b, qVar_b));
        end

        function test_block3_cd_fd_godard_adaptive(testCase)
            label = 'Block 3: CD/MF (overlap-save) + Godard (shared FFT) -> Adaptive CMA';
            [rxSig, symbols] = genRx(testCase);
            runFn = @(ki, kp) eq_clk.combined_cd_fd_godard_adaptive( ...
                rxSig, testCase.SpS, testCase.NFFT, testCase.NOverlap, ...
                testCase.D, testCase.L_km, testCase.CWL, testCase.Rs, ...
                testCase.Rolloff, ki, kp, ...
                testCase.Ns, adaptOpts(testCase));
            [y, ki_b, kp_b, qVar_b] = runSweep(testCase, runFn, ...
                testCase.ki_sweep_godard, testCase.kp_sweep_godard, label);
            evaluate(testCase, rxSig, y, symbols, ...
                sprintf('%s  (best ki=%.2g kp=%.2g qVar=%.2e)', ...
                label, ki_b, kp_b, qVar_b));
        end

        function test_block4_cd_fd_gardner_adaptive(testCase)
            label = 'Block 4: CD/MF (overlap-save) -> Gardner -> Adaptive CMA';
            [rxSig, symbols] = genRx(testCase);
            runFn = @(ki, kp) eq_clk.combined_cd_fd_gardner_adaptive( ...
                rxSig, testCase.SpS, testCase.NFFT, testCase.NOverlap, ...
                testCase.D, testCase.L_km, testCase.CWL, testCase.Rs, ...
                testCase.Rolloff, ki, kp, ...
                testCase.Ns, testCase.NLanes, adaptOpts(testCase));
            [y, ki_b, kp_b, qVar_b] = runSweep(testCase, runFn, ...
                testCase.ki_sweep_gardner, testCase.kp_sweep_gardner, label);
            evaluate(testCase, rxSig, y, symbols, ...
                sprintf('%s  (best ki=%.2g kp=%.2g qVar=%.2e)', ...
                label, ki_b, kp_b, qVar_b));
        end

    end

    % ================================================================
    %  Helpers
    % ================================================================
    methods (Access = private)

        function [bestY, bestKi, bestKp, bestQVar] = runSweep(testCase, ...
                runFn, ki_vec, kp_vec, label)
            % Sweep (ki, kp), evaluate per-quadrant variance of the
            % equaliser output (per polarisation, then averaged), and
            % return the configuration with the lowest variance.
            bestY    = [];
            bestKi   = NaN;
            bestKp   = NaN;
            bestQVar = Inf;

            fprintf('\n=== %s: (ki, kp) sweep [%d x %d] ===\n', ...
                label, numel(ki_vec), numel(kp_vec));
            for ii = 1:numel(ki_vec)
                for jj = 1:numel(kp_vec)
                    ki = ki_vec(ii);
                    kp = kp_vec(jj);
                    try
                        y = runFn(ki, kp);
                    catch ME
                        fprintf('  ki=%.2e kp=%.2e -> ERROR: %s\n', ...
                            ki, kp, ME.message);
                        continue;
                    end
                    if ~all(isfinite(y(:)))
                        fprintf('  ki=%.2e kp=%.2e -> diverged\n', ki, kp);
                        continue;
                    end
                    yEval = trimWarmup(testCase, y);
                    qVar  = quadrantVariance(yEval);
                    fprintf('  ki=%.2e kp=%.2e -> quad-var = %.3e\n', ...
                        ki, kp, qVar);
                    if isfinite(qVar) && qVar < bestQVar
                        bestQVar = qVar;
                        bestY    = y;
                        bestKi   = ki;
                        bestKp   = kp;
                    end
                end
            end
            fprintf('  --> best (ki, kp) = (%.2e, %.2e)  min quad-var = %.3e\n', ...
                bestKi, bestKp, bestQVar);
        end

        function yEval = trimWarmup(testCase, y)
            % Drop the same NOutAdapt-aligned transient used by evaluate(),
            % keeping the sweep metric on the same converged window.
            if size(y, 1) > testCase.NOutAdapt
                yEval = y(testCase.NOutAdapt + 1 : end, :);
            else
                yEval = y;
            end
        end

        function opts = adaptOpts(testCase)
            opts = struct( ...
                'NTaps',          testCase.NTapsAdapt, ...
                'Mu',             testCase.Mu, ...
                'SingleSpike',    testCase.SingleSpike, ...
                'N1',             testCase.N1, ...
                'NOut',           testCase.NOutAdapt, ...
                'SignOnly',       testCase.SignOnly, ...
                'PLanes',         testCase.PLanes, ...
                'Mode',           testCase.Mode, ...
                'Pilots',         [], ...
                'BlockLen',       testCase.BlockLen, ...
                'SubframeBlocks', testCase.SubframeBlocks);
        end

        function [rxSig, symbols] = genRx(testCase)
            % --- Tx: plain DP-QPSK at +/-1 +/-1j, RRC pulse-shaped ----
            symbols = (2*randi([0 1], testCase.Ns, testCase.N_pol) - 1) ...
                + 1j*(2*randi([0 1], testCase.Ns, testCase.N_pol) - 1);

            txSig = modem.rrcPulse(symbols, testCase.SpS, ...
                testCase.Rolloff, testCase.Span);

            % --- Channel: small CD + PMD + timing offset + AWGN -------
            rxSig = channel.add_chromatic_dispersion(txSig, testCase.L_km, ...
                testCase.SpS, testCase.Rs, testCase.D, testCase.CWL);
            rxSig = channel.add_pmd(rxSig, testCase.L_km, testCase.SpS, ...
                testCase.Rs, testCase.DGDSpec, testCase.N_pmd);
            rxSig = channel.apply_timing_error(rxSig, testCase.SFO_ppm, ...
                testCase.tau0, testCase.SpS);
            rxSig = channel.add_awgn(rxSig, testCase.SNR_dB);
        end

        function evaluate(testCase, rxSig, eqSym, symbols, label)
            testCase.verifyTrue(all(isfinite(eqSym(:))), ...
                sprintf('%s: output contains NaN/Inf.', label));

            % Transient discard is owned by NOutAdapt (-> adaptive_eq.NOut);
            % align reference symbols to the post-NOut equaliser output.
            Nout    = size(eqSym, 1);
            eqEval  = eqSym;
            refSyms = symbols(testCase.NOutAdapt + 1 : ...
                              min(testCase.NOutAdapt + Nout, ...
                                  size(symbols, 1)), :);
            Nuse    = min(size(eqEval, 1), size(refSyms, 1));
            eqEval  = eqEval(1:Nuse, :);
            refSyms = refSyms(1:Nuse, :);

            BER = computeBestBER(testCase, refSyms, eqEval);
            fprintf('%s\n   BER = %.3e   (used %d symbols)\n', ...
                label, BER, Nuse);

            % --- Plot before/after constellations -------------------
            rxSym = rxSig(1:testCase.SpS:end, :);
            plotBeforeAfter(testCase, rxSym, eqEval, label, BER);

            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD, ...
                sprintf('%s BER %.2e exceeds threshold %.2e.', ...
                label, BER, testCase.BER_THRESHOLD));
        end

        function bestBER = computeBestBER(testCase, refSyms, eqSym)
            % Resolve per-pol phase rotation (16 angles), pol swap, and a
            % small integer-sample lag (+/- LAG_MAX) ambiguity, return the
            % aggregate lowest BER over both polarisations.
            LAG_MAX = 5;
            totalErrors = 0;
            totalBits   = 0;
            for p = 1:testCase.N_pol
                refBitsPol = modem.symbolsToBits(refSyms(:,p));
                bestPolBER = Inf;
                for q = 1:testCase.N_pol
                    for lag = -LAG_MAX:LAG_MAX
                        if lag >= 0
                            ref = refBitsPol(2*lag+1:end);
                            eqQ = eqSym(1:end-lag, q);
                        else
                            ref = refBitsPol(1:end+2*lag);
                            eqQ = eqSym(-lag+1:end, q);
                        end
                        Nq = min(numel(eqQ), numel(ref)/2);
                        eqQ = eqQ(1:Nq);
                        ref = ref(1:2*Nq);
                        for kk = 0:15
                            rotated = eqQ .* exp(-1j * kk * pi/8);
                            decSym  = modem.decideSymbols(rotated);
                            decBits = modem.symbolsToBits(decSym);
                            polBER  = sum(ref ~= decBits) / numel(ref);
                            if polBER < bestPolBER
                                bestPolBER = polBER;
                            end
                        end
                    end
                end
                totalErrors = totalErrors + bestPolBER * numel(refBitsPol);
                totalBits   = totalBits + numel(refBitsPol);
            end
            bestBER = totalErrors / totalBits;
        end

        function plotBeforeAfter(testCase, rxSym, eqSym, titleStr, BER)
            figure('Name', titleStr, 'Position', [100 100 1200 500]);
            for p = 1:testCase.N_pol
                subplot(2, 2, (p-1)*2 + 1);
                plot(real(rxSym(:,p)), imag(rxSym(:,p)), '.', ...
                     'MarkerSize', 2);
                grid on; axis equal;
                title(sprintf('Before  -  Pol %d', p));
                xlabel('I'); ylabel('Q');

                subplot(2, 2, (p-1)*2 + 2);
                plot(real(eqSym(:,p)), imag(eqSym(:,p)), '.', ...
                     'MarkerSize', 2);
                grid on; axis equal;
                title(sprintf('After  -  Pol %d', p));
                xlabel('I'); ylabel('Q');
            end
            sgtitle(sprintf('%s  |  BER = %.2e', titleStr, BER));
        end

    end
end


% =====================================================================
function v = quadrantVariance(sym)
%QUADRANTVARIANCE  Mean cluster variance about per-quadrant centroids.
%   For multi-column input the per-polarisation variances are averaged.
%   Invariant to +/-pi/2 phase ambiguity; grows monotonically as timing
%   residual smears the clusters.

    if size(sym, 2) > 1
        accV = 0; nValid = 0;
        for p = 1:size(sym, 2)
            vp = quadrantVariance(sym(:, p));
            if isfinite(vp)
                accV   = accV + vp;
                nValid = nValid + 1;
            end
        end
        if nValid == 0
            v = NaN;
        else
            v = accV / nValid;
        end
        return;
    end

    sym = sym(:);
    P = mean(abs(sym).^2);
    if P > 0
        sym = sym / sqrt(P);
    end

    qIdx = (real(sym) > 0) + 2 * (imag(sym) > 0);   % 0..3
    accum = 0;
    cnt   = 0;
    for q = 0:3
        mask = (qIdx == q);
        n    = sum(mask);
        if n < 2
            continue;
        end
        cluster  = sym(mask);
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
