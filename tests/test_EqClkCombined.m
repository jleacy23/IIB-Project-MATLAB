classdef test_EqClkCombined < matlab.unittest.TestCase
    %TEST_EQCLKCOMBINED  Isolation harness for eq_clk.combined_cd_fd_godard_adaptive.
    %
    %   Focused on the production sweep's failing Godard config (NFFT = 128,
    %   NOverlap = 22, CD + PMD, dual pol, SFO, CFO = 0 with in-block coarse
    %   CFO, SNR 22) with the candidate fix applied: 3 CMA taps and 99th-pct
    %   normalisation.  The G-series toggles one suspect knob at a time:
    %     G0  candidate config as-is (NTaps = 3, normalise 99)
    %     G1  normalise 95  (clips ~5% of samples - distorts the metric)
    %     G2  normalise off
    %     G3  NTaps = 1     (the original failing single-tap baseline)
    %
    %   Each test sweeps the PI loop-filter gains and reports the best
    %   per-quadrant cluster variance (phase-ambiguity invariant) and the
    %   recovered BER, with a before/after constellation plot.  The in-block
    %   coarse-CFO integer-bin shift is removed afterwards (mirrors
    %   combined_eq_clk_sweep.removeKnownCFO) so the constellation resolves.
    %
    %   GAIN SCALING.  The combined block uses fft.fft_flp whose forward
    %   transform is 1/N-normalised, so the Godard metric S ~ 1/N^2.  The base
    %   gains bracket the standalone recovery_godard optimum (ki ~ 3.16e-5,
    %   kp ~ 1e-6); they are multiplied by NFFT^2 to keep tau invariant.

    properties (Constant)
        % --- System ---
        Rs      = 32          % [GBd]
        CWL     = 1550        % [nm]
        DGDSpec = 0.1         % PMD [ps/sqrt(km)]
        N_pmd   = 1
        SpS     = 2
        Ns      = 2^14        % symbols per polarisation

        % --- Pulse shaping ---
        Rolloff = 0.25
        Span    = 10

        % --- Adaptive EQ (held fixed; enough taps to mop up residual ISI
        %     so the metric reflects the timing loop, not CMA starvation) ---
        NTapsAdapt     = 7
        Mu             = 1e-3
        SingleSpike    = true
        N1             = 100
        NOutAdapt      = 200
        SignOnly       = true
        PLanes         = 32
        Mode           = 0
        BlockLen       = []
        SubframeBlocks = 0

        % --- Loop-filter gain sweep (BASE values, unnormalised) ----------
        %  Standalone recovery_godard optimum is ki ~ 3.16e-5, kp ~ 1e-6;
        %  these bracket it.  Combined scenarios apply a *NFFT^2 rescale.
        ki_base = logspace(-5.5, -3.5, 5)
        kp_base = logspace(-7.0, -5.0, 5)

        % --- Structure sweep ---
        BER_THRESHOLD = 5e-2
    end

    methods (TestMethodSetup)
        function seedRng(~)
            rng(42);
        end
    end

    % ================================================================
    %  Tests - ordered from most-isolated to most-realistic
    % ================================================================
    methods (Test)

        % --- G. Production sweep config (NFFT=128) with the candidate fix
        %     (3 CMA taps, 99.9-pct normalise), plus one knob toggled each --
        function test_G0_candidate(testCase)
            sc = sweepReplicaScenario(testCase);
            sc.label = 'G0. candidate  NFFT128 NOv22 NTaps3 norm99 CD+PMD';
            sweepCombined(testCase, sc);
        end

        function test_G1_normalise_95(testCase)
            sc = sweepReplicaScenario(testCase);
            sc.normPct = 95;            % old clipping (~5% of samples)
            sc.label = 'G1. candidate but normalise 95 pct';
            sweepCombined(testCase, sc);
        end

        function test_G2_normalise_off(testCase)
            sc = sweepReplicaScenario(testCase);
            sc.normPct = 0;             % no clipping at all
            sc.label = 'G2. candidate but no normalisation';
            sweepCombined(testCase, sc);
        end

        function test_G3_ntaps1_baseline(testCase)
            sc = sweepReplicaScenario(testCase);
            sc.NTaps = 1;               % original failing single-tap baseline
            sc.label = 'G3. candidate but NTaps=1 (failing baseline)';
            sweepCombined(testCase, sc);
        end

    end

    % ================================================================
    %  Scenario / sweep helpers
    % ================================================================
    methods (Access = private)

        function sc = baseScenario(testCase)
            % Cleanest baseline: single pol, no CD, no PMD, no CFO, SFO only.
            sc.N_pol     = 1;
            sc.Rs        = testCase.Rs;
            sc.D         = 0;
            sc.L_km      = 80;
            sc.usePMD    = false;
            sc.CFO_GHz   = 0;
            sc.cfoEnable = false;
            sc.SFO_ppm   = 40;
            sc.tau0      = 0;
            sc.SNR_dB    = 30;
            sc.NFFT      = 256;
            sc.NOverlap  = 64;
            sc.NTaps     = testCase.NTapsAdapt;
            sc.normPct   = 0;       % 0 -> no unit-box normalisation
            sc.label     = '';
        end

        function sc = sweepReplicaScenario(testCase)
            % combined_eq_clk_sweep's Godard tune point (CD + PMD, dual pol,
            % SFO, CFO = 0 with in-block coarse CFO, NFFT = 128, NOverlap = 22
            % from NCD = 22, SNR 22) with the candidate fix applied: 3 CMA
            % taps and 99.9th-pct unit-box normalisation.  G0 runs this as-is;
            % G1/G2 flip the normalise percentile, G3 drops back to 1 tap.
            sc = baseScenario(testCase);
            sc.Rs        = 30.5;
            sc.D         = 20;
            sc.L_km      = 80;
            sc.N_pol     = 2;
            sc.usePMD    = true;
            sc.CFO_GHz   = 0;
            sc.cfoEnable = true;
            sc.SFO_ppm   = 40;
            sc.SNR_dB    = 22;
            sc.NFFT      = 128;
            sc.NOverlap  = 22;
            sc.NTaps     = 3;
            sc.normPct   = 99;
        end

        function [rxSig, symbols] = genRx(testCase, sc)
            symbols = (2*randi([0 1], testCase.Ns, sc.N_pol) - 1) ...
                + 1j*(2*randi([0 1], testCase.Ns, sc.N_pol) - 1);

            rxSig = modem.rrcPulse(symbols, testCase.SpS, ...
                testCase.Rolloff, testCase.Span);

            if sc.D ~= 0
                rxSig = channel.add_chromatic_dispersion(rxSig, sc.L_km, ...
                    testCase.SpS, sc.Rs, sc.D, testCase.CWL);
            end
            if sc.usePMD && sc.N_pol == 2
                rxSig = channel.add_pmd(rxSig, sc.L_km, testCase.SpS, ...
                    sc.Rs, testCase.DGDSpec, testCase.N_pmd);
            end
            if sc.CFO_GHz ~= 0
                rxSig = channel.lo_freq_shift(rxSig, sc.CFO_GHz * 1000, ...
                    sc.Rs, testCase.SpS);
            end
            rxSig = channel.apply_timing_error(rxSig, sc.SFO_ppm, ...
                sc.tau0, testCase.SpS);
            rxSig = channel.add_awgn(rxSig, sc.SNR_dB);

            % Optional unit-box normalisation (CLIPS beyond the pct-th
            % percentile), mirroring combined_eq_clk_sweep's receiver input.
            if sc.normPct > 0
                rxSig = modem.normalise(rxSig, sc.normPct);
            end
        end

        function sweepCombined(testCase, sc)
            % Run the combined block over the NFFT^2-scaled gain grid.
            [rxSig, symbols] = genRx(testCase, sc);
            kiVec = testCase.ki_base * sc.NFFT^2;
            kpVec = testCase.kp_base * sc.NFFT^2;
            runFn = @(ki, kp) runCombined(testCase, sc, rxSig, ki, kp);
            % Combined block discards NOut symbols internally -> output
            % symbol 1 aligns to reference symbol NOutAdapt+1.
            runAndReport(testCase, runFn, kiVec, kpVec, rxSig, symbols, ...
                sc.label, testCase.NOutAdapt + 1);
        end

        function y = runCombined(testCase, sc, rxSig, ki, kp)
            [y, cfoBins] = eq_clk.combined_cd_fd_godard_adaptive( ...
                rxSig, testCase.SpS, sc.NFFT, sc.NOverlap, ...
                sc.D, sc.L_km, testCase.CWL, sc.Rs, testCase.Rolloff, ...
                ki, kp, testCase.Ns, adaptOpts(testCase, sc.NTaps), ...
                sc.cfoEnable, false);
            % Strip the residual CFO that the in-block coarse correction left
            % behind (the downstream carrier-recovery's job) so the
            % constellation resolves instead of forming a ring.  Runs whenever
            % the block applied a coarse shift (cfoEnable) OR a true CFO is
            % present: residual = true CFO - the integer-bin estimate applied.
            % Note coarse_cfo_fd applies a (small) integer-bin shift even when
            % the true CFO is 0, so this must run for cfoEnable too.  Mirrors
            % combined_eq_clk_sweep.removeKnownCFO (which runs unconditionally).
            if sc.cfoEnable || sc.CFO_GHz ~= 0
                binGHz   = testCase.SpS * sc.Rs / sc.NFFT;
                residGHz = sc.CFO_GHz - cfoBins * binGHz;
                n = (0 : size(y, 1) - 1).';
                y = y .* exp(-1j * 2*pi * (residGHz / sc.Rs) * n);
            end
        end

        function opts = adaptOpts(testCase, nTaps)
            opts = struct( ...
                'NTaps',          nTaps, ...
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

        function runAndReport(testCase, runFn, kiVec, kpVec, ...
                rxSig, symbols, label, firstRefIdx)
            [bestY, bestKi, bestKp, bestQVar] = ...
                runSweep(testCase, runFn, kiVec, kpVec, label);
            testCase.verifyTrue(~isempty(bestY) && isfinite(bestQVar), ...
                sprintf('%s: no finite-variance configuration found.', label));
            evaluate(testCase, rxSig, bestY, symbols, ...
                sprintf('%s  (best ki=%.3g kp=%.3g qVar=%.2e)', ...
                label, bestKi, bestKp, bestQVar), firstRefIdx);
        end

        function [bestY, bestKi, bestKp, bestQVar] = runSweep(testCase, ...
                runFn, ki_vec, kp_vec, label)
            bestY = []; bestKi = NaN; bestKp = NaN; bestQVar = Inf;
            fprintf('\n=== %s: (ki,kp) sweep [%d x %d] ===\n', ...
                label, numel(ki_vec), numel(kp_vec));
            for ii = 1:numel(ki_vec)
                for jj = 1:numel(kp_vec)
                    ki = ki_vec(ii); kp = kp_vec(jj);
                    try
                        y = runFn(ki, kp);
                    catch ME
                        fprintf('  ki=%.2e kp=%.2e -> ERROR: %s\n', ...
                            ki, kp, ME.message);
                        continue;
                    end
                    if isempty(y) || ~all(isfinite(y(:)))
                        fprintf('  ki=%.2e kp=%.2e -> diverged\n', ki, kp);
                        continue;
                    end
                    yEval = trimWarmup(testCase, downsampleSym(testCase, y));
                    qVar  = quadrantVariance(yEval);
                    fprintf('  ki=%.2e kp=%.2e -> quad-var = %.3e\n', ...
                        ki, kp, qVar);
                    if isfinite(qVar) && qVar < bestQVar
                        bestQVar = qVar; bestY = y; bestKi = ki; bestKp = kp;
                    end
                end
            end
            fprintf('  --> best (ki,kp) = (%.3g, %.3g)  min quad-var = %.3e\n', ...
                bestKi, bestKp, bestQVar);
        end

        function s = downsampleSym(testCase, y)
            % The combined block returns symbol-rate output; the standalone
            % reference returns 2-Sa/sym.  Downsample only the oversampled
            % case so the metric always sees one sample per symbol.
            if size(y, 1) > 1.5 * testCase.Ns
                s = y(1:testCase.SpS:end, :);
            else
                s = y;
            end
        end

        function yEval = trimWarmup(testCase, y)
            if size(y, 1) > testCase.NOutAdapt
                yEval = y(testCase.NOutAdapt + 1 : end, :);
            else
                yEval = y;
            end
        end

        function evaluate(testCase, rxSig, eqOut, symbols, label, firstRefIdx)
            eqSym = downsampleSym(testCase, eqOut);
            testCase.verifyTrue(all(isfinite(eqSym(:))), ...
                sprintf('%s: output contains NaN/Inf.', label));

            % Drop a common warm-up window so loop/CMA acquisition does not
            % pollute the BER, then align the reference: eqSym(1) maps to
            % symbol firstRefIdx (NOutAdapt+1 for the combined block, which
            % discards NOut internally; 1 for the standalone reference).
            W = testCase.NOutAdapt;
            if size(eqSym, 1) > W
                eqSym    = eqSym(W+1:end, :);
                firstRef = firstRefIdx + W;
            else
                firstRef = firstRefIdx;
            end

            Nout    = size(eqSym, 1);
            lastRef = min(firstRef + Nout - 1, size(symbols, 1));
            refSyms = symbols(firstRef:lastRef, :);

            nPol    = min(size(eqSym, 2), size(refSyms, 2));
            eqSym   = eqSym(:, 1:nPol);
            refSyms = refSyms(:, 1:nPol);
            Nuse    = min(size(eqSym, 1), size(refSyms, 1));
            eqSym   = eqSym(1:Nuse, :);
            refSyms = refSyms(1:Nuse, :);

            BER = computeBestBER(testCase, refSyms, eqSym);
            fprintf('%s\n   BER = %.3e   (used %d symbols, %d pol)\n', ...
                label, BER, Nuse, nPol);

            rxSym = rxSig(1:testCase.SpS:end, 1:nPol);
            plotBeforeAfter(rxSym, eqSym, label, BER);

            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD, ...
                sprintf('%s BER %.2e exceeds threshold %.2e.', ...
                label, BER, testCase.BER_THRESHOLD));
        end

        function bestBER = computeBestBER(~, refSyms, eqSym)
            % Per-pol phase rotation (16 angles), pol swap, and +/-5 sample
            % lag ambiguity; aggregate lowest BER over all polarisations.
            LAG_MAX = 5;
            nPol = size(eqSym, 2);
            totalErrors = 0; totalBits = 0;
            for p = 1:nPol
                refBitsPol = modem.symbolsToBits(refSyms(:,p));
                bestPolBER = Inf;
                for q = 1:nPol
                    for lag = -LAG_MAX:LAG_MAX
                        if lag >= 0
                            ref = refBitsPol(2*lag+1:end);
                            eqQ = eqSym(1:end-lag, q);
                        else
                            ref = refBitsPol(1:end+2*lag);
                            eqQ = eqSym(-lag+1:end, q);
                        end
                        Nq  = min(numel(eqQ), numel(ref)/2);
                        eqQ = eqQ(1:Nq);
                        ref = ref(1:2*Nq);
                        for kk = 0:15
                            rotated = eqQ .* exp(-1j * kk * pi/8);
                            decBits = modem.symbolsToBits(modem.decideSymbols(rotated));
                            polBER  = sum(ref ~= decBits) / numel(ref);
                            if polBER < bestPolBER, bestPolBER = polBER; end
                        end
                    end
                end
                totalErrors = totalErrors + bestPolBER * numel(refBitsPol);
                totalBits   = totalBits + numel(refBitsPol);
            end
            bestBER = totalErrors / totalBits;
        end

    end
end


% =====================================================================
%  Local functions
% =====================================================================
function plotBeforeAfter(rxSym, eqSym, titleStr, BER)
    nPol = size(eqSym, 2);
    figure('Name', titleStr, 'Position', [100 100 1100 450*nPol]);
    for p = 1:nPol
        subplot(nPol, 2, (p-1)*2 + 1);
        plot(real(rxSym(:,p)), imag(rxSym(:,p)), '.', 'MarkerSize', 2);
        grid on; axis equal;
        title(sprintf('Before - Pol %d', p)); xlabel('I'); ylabel('Q');

        subplot(nPol, 2, (p-1)*2 + 2);
        e = eqSym(:,p); e = e / sqrt(mean(abs(e).^2));
        plot(real(e), imag(e), '.', 'MarkerSize', 2);
        grid on; axis equal; xlim([-2 2]); ylim([-2 2]);
        title(sprintf('After - Pol %d', p)); xlabel('I'); ylabel('Q');
    end
    sgtitle(sprintf('%s  |  BER = %.2e', titleStr, BER), 'Interpreter', 'none');
end


function v = quadrantVariance(sym)
%QUADRANTVARIANCE  Mean cluster variance about per-quadrant centroids.
%   Multi-column input -> mean of per-polarisation variances.  Invariant to
%   +/-pi/2 phase ambiguity; grows as timing residual smears the clusters.
    if size(sym, 2) > 1
        accV = 0; nValid = 0;
        for p = 1:size(sym, 2)
            vp = quadrantVariance(sym(:, p));
            if isfinite(vp), accV = accV + vp; nValid = nValid + 1; end
        end
        if nValid == 0, v = NaN; else, v = accV / nValid; end
        return;
    end

    sym = sym(:);
    P = mean(abs(sym).^2);
    if P > 0, sym = sym / sqrt(P); end

    qIdx = (real(sym) > 0) + 2 * (imag(sym) > 0);
    accum = 0; cnt = 0;
    for q = 0:3
        mask = (qIdx == q);
        if sum(mask) < 2, continue; end
        cluster = sym(mask);
        accum = accum + sum(abs(cluster - mean(cluster)).^2);
        cnt   = cnt + numel(cluster);
    end
    if cnt == 0, v = NaN; else, v = accum / cnt; end
end
