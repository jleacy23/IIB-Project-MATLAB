classdef combined_eq_clk_sweep < matlab.unittest.TestCase
%COMBINED_EQ_CLK_SWEEP  Two-stage comparison of overlap-save combined
%   equalisation + clock-recovery implementations under worst-case CPON
%   tolerances.
%
%   Blocks under test:
%     1. cd_gardner_cma       eq_clk.combined_cd_fd_gardner_adaptive
%                             overlap-save CD/MF -> Gardner DPLL -> CMA
%
%     2. cd_godard_cma        eq_clk.combined_cd_fd_godard_adaptive
%                             overlap-save CD/MF with embedded Modified
%                             Godard loop (excess-band metric) -> CMA
%
%   Both apply the same overlap-save CD/matched-filter stage, so the
%   downstream CMA only sees PMD and residual ISI.
%
%   The class exposes two independent tests:
%
%     test_design_sweep
%        CFO is held at 0.  For every (block, NCD, NTaps) configuration,
%        the (ki, kp) loop-filter grid is searched to minimise BER at
%        SNR_dB_tune, then a BER-vs-SNR Monte-Carlo is run.  The output
%        records the per-configuration BER cube, the tuning grid, and a
%        suggested "best" design per block (lowest FEC SNR) for the user
%        to inspect.  Saved to combined_eq_clk_design_sweep.mat.
%
%     test_cfo_sweep
%        Stand-alone CFO sweep.  Takes the per-block (NCD, NTaps,
%        ki, kp) from the Cfo_* class constants -- which the user fills
%        in manually, typically after inspecting the design-sweep
%        results -- and sweeps BER vs (SNR, CFO).  Does not read the
%        design-sweep .mat.  Saved to combined_eq_clk_cfo_sweep.mat.
%
%   Channel chain (CPON tolerances):
%       Tx -> RRC -> CD -> PMD -> CFO -> SFO -> AWGN -> [block] ->
%       exact CFO correction (post-block) -> phase/pol resolution -> BER.
%   Each combined block now performs its own coarse FD CFO correction
%   (eq_clk.coarse_cfo_fd, integer-bin spectral-centroid shift) before
%   timing recovery; the post-block "exact" CFO removal only absorbs
%   the residual sub-bin CFO so the constellation can be resolved.
%
%   Run with:
%       runtests('combined_eq_clk_sweep')

    %% ================================================================
    %  Parameters
    %% ================================================================
    properties (Constant)

        % --- System (CPON downstream specification) -----------------
        Rs       = 30.5            % [GBd]
        SpS      = 2               % Gardner DPLL requires 2 Sa/sym
        N_pol    = 2
        D        = 20              % [ps/(nm*km)]
        CWL      = 1550            % [nm]
        DGDSpec  = 0.1             % [ps/sqrt(km)]
        N_pmd    = 1
        PMD_seed = 12345           % fixed seed: PMD draw is identical
                                   % across every trial, SNR, CFO, and
                                   % block under test

        % --- CPON tolerances (report worst case) --------------------
        L_km        = 80           % worst-case CPON reach
        SFO_ppm     = 40           % worst-case sample-frequency offset

        % CFO axes.  test_design_sweep runs at CFO = 0; test_cfo_sweep
        % uses CFO_GHz_vec as its CFO axis with the per-block optimal
        % design selected by test_design_sweep.  Loop-filter tuning is
        % always performed at CFO_GHz_tune.
        CFO_GHz_tune = 0
        CFO_GHz_vec  = 0:0.5:3.0

        % --- Pulse shaping ------------------------------------------
        Rolloff  = 0.25
        Span     = 10

        % --- Monte-Carlo --------------------------------------------
        Ns          = 37500         % symbols per polarisation per trial
        NTrials     = 5
        SNR_dB_vec  = 0 : 2 : 22
        SNR_dB_tune = 22           % SNR used for the (ki, kp) sweep

        % --- Static-equaliser sizing --------------------------------
        % NCD_FD_vec sweeps the CD filter length in test_design_sweep;
        % NOverlap is derived per-config as 2*ceil((NCD-1)/2) so the
        % overlap-save extension is even.
        NFFT        = 128
        NCD_FD_vec  = [22]

        % --- Adaptive-equaliser tap sweep ---------------------------
        NTaps_vec = [1]

        % --- Per-block FFT/IFFT twiddle quantisation sweep ----------
        % Each entry is passed as the po2Twiddle argument to the
        % combined block: false uses MATLAB's native fft/ifft via
        % fft.fft_flp; true snaps every twiddle component to the
        % nearest signed power of two (hardware shift-only).
        Po2Twiddle_vec = [false, true]

        % --- Adaptive equaliser (per-block convergence) -------------
        MuGardner   = 1e-3
        N1Gardner   = 500
        MuGodard    = 1e-3
        N1Godard    = 500

        % --- Adaptive equaliser (shared) ----------------------------
        NOut        = 1000
        SignOnly    = true
        SingleSpike = true
        PLanesAEQ   = 32

        % --- Clock-recovery loop-filter gain sweeps -----------------
        %  Gardner gains operate on the TIME-DOMAIN TED (post static
        %  IFFT); the forward(1/N) x inverse(1) round trip is unchanged by
        %  the FFT renormalisation, so these are unaffected.
        ki_gardner_vec = logspace(-7, -4, 4)
        kp_gardner_vec = logspace(-6, -3, 4)
        %  Godard gains act on the FREQUENCY-DOMAIN metric S, which is
        %  quadratic in the spectrum.  Since the forward FFT is now
        %  1/N-normalised (1/2 per stage), the spectrum scales 1/N and S
        %  scales 1/N^2, so the gains are scaled up by NFFT^2 (= 128^2) to
        %  keep ki*e / kp*e — and the tau trajectory — invariant.
        ki_godard_vec  = logspace(-6, -3, 4) * 128^2
        kp_godard_vec  = logspace(-5, -2, 4) * 128^2
        NLanesGard     = 32

        % --- Coarse FD CFO correction (eq_clk.coarse_cfo_fd) --------
        % One-shot pre-loop CFO correction: estimate from centroid of
        % first NFFT samples' FFT, apply continuous-valued time-domain
        % phasor to the whole input.  Set false to disable.
        CfoEnable = true

        % --- FEC threshold used to score designs --------------------
        FEC_BER = 2e-2

        % --- Manual designs for test_cfo_sweep ----------------------
        % Set these by hand from the test_design_sweep results.
        %   Cfo_blocks / Cfo_NCD / Cfo_NTaps are parallel 1-D arrays,
        %   one entry per design row.
        %   Cfo_ki / Cfo_kp are [NDESIGN x NPO2] matrices, indexed by
        %   (design_row, Po2Twiddle_vec_column), so each po2 variant
        %   can have its own loop-filter gains.  Default: same gains
        %   for both po2 columns.
        %  Godard rows are scaled by NFFT^2 (= 128^2) for the 1/N-normalised
        %  forward FFT (metric ~ 1/N^2); Gardner rows are unchanged.
        Cfo_blocks  = {'cd_gardner_cma', 'cd_godard_cma'}
        Cfo_NCD     = [22,   22]
        Cfo_NTaps   = [1,    1]
        Cfo_ki      = [1e-6       1e-7; ...
                       1e-4*128^2 1e-4*128^2]
        Cfo_kp      = [1e-4       1e-4; ...
                       1e-5*128^2 1e-5*128^2]
    end

    %% ================================================================
    %  Setup
    %% ================================================================
    methods (TestClassSetup)
        function setupPaths(~)
            here = fileparts(mfilename('fullpath'));
            addpath(genpath(fullfile(here, '..', '..', 'src')));
        end
    end

    methods (TestMethodSetup)
        function seedRng(~)
            rng(42);
        end
    end

    %% ================================================================
    %  Tests
    %% ================================================================
    methods (Test)

        function test_design_sweep(testCase)
            % Design sweep at CFO = 0: full (block, NCD, NTaps, ki, kp)
            % search to identify the best design per block.
            P    = combined_eq_clk_sweep.extractParams(testCase);
            cfgs = combined_eq_clk_sweep.buildConfigs(P);
            NCFG = numel(cfgs);
            NSNR = numel(P.SNR_dB_vec);
            cfo  = 0;

            % --- Pre-allocate table columns -----------------------
            block_name = strings(NCFG, 1);
            l_km       = nan(NCFG, 1);
            n_cd       = nan(NCFG, 1);
            n_overlap  = nan(NCFG, 1);
            n_fft      = nan(NCFG, 1);
            n_aeq      = nan(NCFG, 1);
            po2        = false(NCFG, 1);
            ki         = nan(NCFG, 1);
            kp         = nan(NCFG, 1);
            tune_ki    = cell(NCFG, 1);
            tune_kp    = cell(NCFG, 1);
            tune_ber   = cell(NCFG, 1);
            ber        = cell(NCFG, 1);
            fec_snr    = nan(NCFG, 1);

            % --- Stage 1: tune (ki, kp) per (block, NCD, NTaps, po2) ---
            fprintf('\n=== Design sweep, Stage 1: (ki, kp) tuning ===\n');
            for ci = 1:NCFG
                cfg = cfgs(ci);
                [bestKi, bestKp, kiVec, kpVec, berGrid] = ...
                    combined_eq_clk_sweep.tuneLoopFilter(P, cfg);
                cfgs(ci).ki  = bestKi;
                cfgs(ci).kp  = bestKp;
                tune_ki{ci}  = kiVec;
                tune_kp{ci}  = kpVec;
                tune_ber{ci} = berGrid;
                fprintf(['  %-16s NCD=%2d NTaps=%2d po2=%d -> ', ...
                    '(ki, kp) = (%.2e, %.2e), tune BER = %.2e\n'], ...
                    cfg.block, cfg.NCD, cfg.NTaps, cfg.Po2Twiddle, ...
                    bestKi, bestKp, min(berGrid(:)));
            end

            % --- Stage 2: BER vs SNR at CFO = 0 ---------------------
            fprintf('\n=== Design sweep, Stage 2: BER vs SNR ===\n');
            for ci = 1:NCFG
                ber{ci} = nan(P.NTrials, NSNR);
            end
            for tr = 1:P.NTrials
                fprintf('--- trial %d / %d ---\n', tr, P.NTrials);
                for si = 1:NSNR
                    snr = P.SNR_dB_vec(si);
                    [rxSig, symbols] = ...
                        combined_eq_clk_sweep.buildChannel( ...
                            P, snr, tr, cfo);
                    for ci = 1:NCFG
                        b = combined_eq_clk_sweep.runOnePoint( ...
                                P, cfgs(ci), rxSig, symbols, cfo);
                        ber{ci}(tr, si) = b;
                    end
                end
            end

            % --- Score by mean FEC SNR ----------------------------
            for ci = 1:NCFG
                fec_snr(ci) = combined_eq_clk_sweep.fecSnrFromBer( ...
                    P.SNR_dB_vec, mean(ber{ci}, 1, 'omitnan'), P.FEC_BER);
            end

            % --- Build output table -------------------------------
            for ci = 1:NCFG
                cfg = cfgs(ci);
                block_name(ci) = string(cfg.block);
                l_km(ci)       = P.L_km;
                ki(ci)         = cfg.ki;
                kp(ci)         = cfg.kp;
                n_aeq(ci)      = cfg.NTaps;
                n_cd(ci)       = cfg.NCD;
                n_overlap(ci)  = cfg.NOverlap;
                n_fft(ci)      = P.NFFT;
                po2(ci)        = cfg.Po2Twiddle;
            end
            tbl = table(block_name, l_km, n_cd, n_overlap, n_fft, ...
                n_aeq, po2, ki, kp, tune_ki, tune_kp, tune_ber, ber, ...
                fec_snr);

            % --- Select best (lowest FEC SNR) design per block -----
            best_designs = combined_eq_clk_sweep.selectBestDesigns(tbl);
            fprintf('\n--- Best design per block (lowest FEC SNR) ---\n');
            disp(best_designs(:, {'block_name', 'n_cd', 'n_aeq', ...
                'po2', 'ki', 'kp', 'fec_snr'}));

            % --- Save ---------------------------------------------
            S.tbl          = tbl;
            S.best_designs = best_designs;
            S.SNR_dB_vec   = P.SNR_dB_vec;
            S.SNR_dB_tune  = P.SNR_dB_tune;
            S.CFO_GHz      = cfo;
            S.params       = P;

            outFile = fullfile(fileparts(mfilename('fullpath')), ...
                'combined_eq_clk_design_sweep.mat');
            save(outFile, '-struct', 'S');
            fprintf('Saved design sweep to %s\n', outFile);

            testCase.verifyEqual(height(tbl), NCFG);
        end

        function test_cfo_sweep(testCase)
            % Stand-alone CFO sweep.  Takes the per-block (NCD, NTaps,
            % ki, kp) from the Cfo_* class constants, which the user
            % fills in manually -- typically from the test_design_sweep
            % output -- and sweeps BER over (SNR, CFO, Po2Twiddle).
            P = combined_eq_clk_sweep.extractParams(testCase);

            NDESIGN = numel(P.Cfo_blocks);
            NSNR = numel(P.SNR_dB_vec);
            NCFO = numel(P.CFO_GHz_vec);
            NPO2 = numel(P.Po2Twiddle_vec);
            NCFG = NDESIGN * NPO2;
            assert(numel(P.Cfo_NCD)   == NDESIGN && ...
                   numel(P.Cfo_NTaps) == NDESIGN && ...
                   isequal(size(P.Cfo_ki), [NDESIGN, NPO2]) && ...
                   isequal(size(P.Cfo_kp), [NDESIGN, NPO2]), ...
                ['Cfo_blocks/NCD/NTaps must have length NDESIGN; ', ...
                 'Cfo_ki/Cfo_kp must be [NDESIGN x NPO2].']);

            % Build the cfg array: each manual design is replicated
            % once per Po2Twiddle value, with its own (ki, kp) pair.
            cfgs = struct('block', {}, 'NTaps', {}, 'NCD', {}, ...
                'NOverlap', {}, 'Po2Twiddle', {}, 'ki', {}, 'kp', {});
            ci = 0;
            for di = 1:NDESIGN
                nCd = P.Cfo_NCD(di);
                for pi = 1:NPO2
                    ci = ci + 1;
                    cfgs(ci).block      = P.Cfo_blocks{di};
                    cfgs(ci).NTaps      = P.Cfo_NTaps(di);
                    cfgs(ci).NCD        = nCd;
                    cfgs(ci).NOverlap   = 2 * ceil((nCd - 1) / 2);
                    cfgs(ci).Po2Twiddle = logical(P.Po2Twiddle_vec(pi));
                    cfgs(ci).ki         = P.Cfo_ki(di, pi);
                    cfgs(ci).kp         = P.Cfo_kp(di, pi);
                end
            end

            ber = cell(NCFG, 1);
            for ci = 1:NCFG
                ber{ci} = nan(P.NTrials, NSNR, NCFO);
            end

            fprintf('\n=== CFO sweep: BER vs SNR x CFO x Po2Twiddle ===\n');
            for ci = 1:NCFG
                fprintf(['  %-16s NCD=%2d NTaps=%d po2=%d  ', ...
                    '(ki, kp) = (%.2e, %.2e)\n'], ...
                    cfgs(ci).block, cfgs(ci).NCD, cfgs(ci).NTaps, ...
                    cfgs(ci).Po2Twiddle, cfgs(ci).ki, cfgs(ci).kp);
            end

            for tr = 1:P.NTrials
                fprintf('--- trial %d / %d ---\n', tr, P.NTrials);
                for ic = 1:NCFO
                    cfo = P.CFO_GHz_vec(ic);
                    for si = 1:NSNR
                        snr = P.SNR_dB_vec(si);
                        [rxSig, symbols] = ...
                            combined_eq_clk_sweep.buildChannel( ...
                                P, snr, tr, cfo);
                        for ci = 1:NCFG
                            b = combined_eq_clk_sweep.runOnePoint( ...
                                    P, cfgs(ci), rxSig, symbols, cfo);
                            ber{ci}(tr, si, ic) = b;
                        end
                    end
                end
            end

            % --- Build output table -------------------------------
            block_name = strings(NCFG, 1);
            l_km       = repmat(P.L_km, NCFG, 1);
            n_cd       = nan(NCFG, 1);
            n_overlap  = nan(NCFG, 1);
            n_fft      = repmat(P.NFFT, NCFG, 1);
            n_aeq      = nan(NCFG, 1);
            po2        = false(NCFG, 1);
            ki         = nan(NCFG, 1);
            kp         = nan(NCFG, 1);
            for ci = 1:NCFG
                block_name(ci) = string(cfgs(ci).block);
                n_cd(ci)       = cfgs(ci).NCD;
                n_overlap(ci)  = cfgs(ci).NOverlap;
                n_aeq(ci)      = cfgs(ci).NTaps;
                po2(ci)        = cfgs(ci).Po2Twiddle;
                ki(ci)         = cfgs(ci).ki;
                kp(ci)         = cfgs(ci).kp;
            end
            tbl = table(block_name, l_km, n_cd, n_overlap, n_fft, ...
                n_aeq, po2, ki, kp, ber);

            S.tbl          = tbl;
            S.SNR_dB_vec   = P.SNR_dB_vec;
            S.CFO_GHz_vec  = P.CFO_GHz_vec;
            S.params       = P;

            outFile = fullfile(fileparts(mfilename('fullpath')), ...
                'combined_eq_clk_cfo_sweep.mat');
            save(outFile, '-struct', 'S');
            fprintf('Saved CFO sweep to %s\n', outFile);

            testCase.verifyEqual(height(tbl), NCFG);
        end

    end

    %% ================================================================
    %  Helpers
    %% ================================================================
    methods (Static)

        function P = extractParams(tc)
            P.Rs            = tc.Rs;
            P.SpS           = tc.SpS;
            P.N_pol         = tc.N_pol;
            P.D             = tc.D;
            P.CWL           = tc.CWL;
            P.DGDSpec       = tc.DGDSpec;
            P.N_pmd         = tc.N_pmd;
            P.PMD_seed      = tc.PMD_seed;
            P.L_km          = tc.L_km;
            P.CFO_GHz_tune  = tc.CFO_GHz_tune;
            P.CFO_GHz_vec   = tc.CFO_GHz_vec;
            P.SFO_ppm       = tc.SFO_ppm;
            P.Rolloff       = tc.Rolloff;
            P.Span          = tc.Span;
            P.Ns            = tc.Ns;
            P.NTrials       = tc.NTrials;
            P.SNR_dB_vec    = tc.SNR_dB_vec;
            P.SNR_dB_tune   = tc.SNR_dB_tune;
            P.NCD_FD_vec    = tc.NCD_FD_vec;
            P.NFFT          = tc.NFFT;
            P.NTaps_vec     = tc.NTaps_vec;
            P.Po2Twiddle_vec = tc.Po2Twiddle_vec;
            P.MuGardner     = tc.MuGardner;
            P.N1Gardner     = tc.N1Gardner;
            P.MuGodard      = tc.MuGodard;
            P.N1Godard      = tc.N1Godard;
            P.NOut          = tc.NOut;
            P.SignOnly      = tc.SignOnly;
            P.SingleSpike   = tc.SingleSpike;
            P.PLanesAEQ     = tc.PLanesAEQ;
            P.ki_gardner_vec = tc.ki_gardner_vec;
            P.kp_gardner_vec = tc.kp_gardner_vec;
            P.ki_godard_vec  = tc.ki_godard_vec;
            P.kp_godard_vec  = tc.kp_godard_vec;
            P.NLanesGard    = tc.NLanesGard;
            P.CfoEnable     = tc.CfoEnable;
            P.FEC_BER       = tc.FEC_BER;
            P.Cfo_blocks    = tc.Cfo_blocks;
            P.Cfo_NCD       = tc.Cfo_NCD;
            P.Cfo_NTaps     = tc.Cfo_NTaps;
            P.Cfo_ki        = tc.Cfo_ki;
            P.Cfo_kp        = tc.Cfo_kp;
        end

        function cfgs = buildConfigs(P)
            % Flat list of (block, NCD, NTaps, Po2Twiddle) configurations.
            % NOverlap is rounded up to the next even integer so
            % overlap-save can pad symmetrically.
            blocks = {'cd_gardner_cma', 'cd_godard_cma'};
            cfgs = struct('block', {}, 'NTaps', {}, 'NCD', {}, ...
                'NOverlap', {}, 'Po2Twiddle', {}, 'ki', {}, 'kp', {});
            for bi = 1:numel(blocks)
                for nCd = P.NCD_FD_vec
                    nOver = 2 * ceil((nCd - 1) / 2);
                    for n = P.NTaps_vec
                        for po2 = P.Po2Twiddle_vec
                            cfgs(end+1).block      = blocks{bi}; %#ok<AGROW>
                            cfgs(end).NTaps      = n;
                            cfgs(end).NCD        = nCd;
                            cfgs(end).NOverlap   = nOver;
                            cfgs(end).Po2Twiddle = logical(po2);
                            cfgs(end).ki         = NaN;
                            cfgs(end).kp         = NaN;
                        end
                    end
                end
            end
        end

        function [bestKi, bestKp, kiVec, kpVec, berGrid] = ...
                tuneLoopFilter(P, cfg)
            % Single-trial (ki, kp) grid search at SNR_dB_tune and
            % CFO_GHz_tune.
            switch cfg.block
                case 'cd_gardner_cma'
                    kiVec = P.ki_gardner_vec;
                    kpVec = P.kp_gardner_vec;
                case 'cd_godard_cma'
                    kiVec = P.ki_godard_vec;
                    kpVec = P.kp_godard_vec;
            end
            berGrid = nan(numel(kiVec), numel(kpVec));

            [rxSig, symbols] = combined_eq_clk_sweep.buildChannel( ...
                P, P.SNR_dB_tune, 0, P.CFO_GHz_tune);

            for ii = 1:numel(kiVec)
                for jj = 1:numel(kpVec)
                    cfgT = cfg;
                    cfgT.ki = kiVec(ii);
                    cfgT.kp = kpVec(jj);
                    berGrid(ii, jj) = ...
                        combined_eq_clk_sweep.runOnePoint( ...
                            P, cfgT, rxSig, symbols, P.CFO_GHz_tune);
                end
            end
            [bMin, idx] = min(berGrid(:));
            if ~isfinite(bMin)
                bestKi = kiVec(1);
                bestKp = kpVec(1);
                return;
            end
            [ii, jj] = ind2sub(size(berGrid), idx);
            bestKi = kiVec(ii);
            bestKp = kpVec(jj);
        end

        function best = selectBestDesigns(tbl)
            % Pick the row with the lowest finite fec_snr for each
            % unique block_name.  Falls back to lowest mean BER over
            % the swept SNRs if no row reaches the FEC threshold.
            blocks = unique(tbl.block_name, 'stable');
            keepIdx = nan(numel(blocks), 1);
            for bi = 1:numel(blocks)
                sub = find(tbl.block_name == blocks(bi));
                sFec = tbl.fec_snr(sub);
                if any(isfinite(sFec))
                    [~, k] = min(sFec);
                else
                    fb = arrayfun(@(r) ...
                        mean(tbl.ber{r}(:), 'omitnan'), sub);
                    [~, k] = min(fb);
                end
                keepIdx(bi) = sub(k);
            end
            best = tbl(keepIdx, :);
        end

        function snr = fecSnrFromBer(snrVec, berVec, fecBer)
            % Linear-on-log interpolation of the BER curve to find the
            % SNR at which BER == fecBer.  Returns NaN if the curve does
            % not bracket the FEC threshold.
            v = berVec(:).';
            v(v <= 0) = NaN;
            valid = isfinite(v);
            if nnz(valid) < 2
                snr = NaN;
                return;
            end
            x = snrVec(valid);
            y = log10(v(valid));
            target = log10(fecBer);
            if target < min(y) || target > max(y)
                snr = NaN;
                return;
            end
            % Monotone-ish: find first crossing point.
            for k = 1:numel(x) - 1
                if (y(k) - target) * (y(k+1) - target) <= 0
                    snr = x(k) + (target - y(k)) / (y(k+1) - y(k)) * ...
                        (x(k+1) - x(k));
                    return;
                end
            end
            snr = NaN;
        end

        function ber = runOnePoint(P, cfg, rxSig, symbols, cfo)
            % One (block, NTaps, ki, kp) BER measurement on a fixed
            % channel realisation at the supplied CFO (GHz).
            try
                [eqSym, cfoBinsApplied] = ...
                    combined_eq_clk_sweep.runBlock(P, cfg, rxSig);
            catch
                ber = NaN;
                return;
            end
            if isempty(eqSym) || ~all(isfinite(eqSym(:)))
                ber = NaN;
                return;
            end
            eqSym = combined_eq_clk_sweep.removeKnownCFO( ...
                eqSym, cfo, cfoBinsApplied, P);
            ber = combined_eq_clk_sweep.computeBER( ...
                eqSym, symbols, P.NOut);
        end

        function [rxSig, symbols] = buildChannel(P, SNR_dB, trialSeed, cfo)
            % Deterministic per-(trial, SNR, CFO) realisation.
            % trialSeed = 0 is the dedicated tuning realisation.
            rng(1000 * trialSeed + round(SNR_dB) + 13);

            symbols = (2*randi([0 1], P.Ns, P.N_pol) - 1) ...
                + 1j * (2*randi([0 1], P.Ns, P.N_pol) - 1);

            txSig = modem.rrcPulse(symbols, P.SpS, P.Rolloff, P.Span);

            rxSig = channel.add_chromatic_dispersion(txSig, P.L_km, ...
                P.SpS, P.Rs, P.D, P.CWL);

            % PMD: held fixed across every run by isolating its RNG
            % consumption with a dedicated seed.  The surrounding state
            % is preserved so the symbol/AWGN draws stay deterministic
            % in (trialSeed, SNR_dB).
            rngState = rng;
            rng(P.PMD_seed);
            rxSig = channel.add_pmd(rxSig, P.L_km, P.SpS, P.Rs, ...
                P.DGDSpec, P.N_pmd);
            rng(rngState);

            % cfo (GHz) * 1000 = MHz for lo_freq_shift
            rxSig = channel.lo_freq_shift(rxSig, cfo * 1000, ...
                P.Rs, P.SpS);
            rxSig = channel.apply_timing_error(rxSig, P.SFO_ppm, ...
                0, P.SpS);
            rxSig = channel.add_awgn(rxSig, SNR_dB);

            % Normalise into the unit box before the receiver, using the
            % 95th-percentile magnitude as the per-pol scale reference.
            rxSig = modem.normalise(rxSig, 95);
        end

        function [eqSym, cfoBinsApplied] = runBlock(P, cfg, rxSig)
            switch cfg.block
                case 'cd_gardner_cma'
                    mu = P.MuGardner;
                    n1 = P.N1Gardner;
                case 'cd_godard_cma'
                    mu = P.MuGodard;
                    n1 = P.N1Godard;
                otherwise
                    error('combined_eq_clk_sweep:badBlock', ...
                          'Unknown block: %s', cfg.block);
            end
            adaptOpts = struct( ...
                'NTaps',       cfg.NTaps, ...
                'Mu',          mu, ...
                'SingleSpike', P.SingleSpike, ...
                'N1',          n1, ...
                'NOut',        P.NOut, ...
                'SignOnly',    P.SignOnly, ...
                'PLanes',      P.PLanesAEQ);
            cfoBinsApplied = 0;
            switch cfg.block
                case 'cd_gardner_cma'
                    [eqSym, cfoBinsApplied] = ...
                        eq_clk.combined_cd_fd_gardner_adaptive( ...
                        rxSig, P.SpS, P.NFFT, cfg.NOverlap, P.D, P.L_km, ...
                        P.CWL, P.Rs, P.Rolloff, cfg.ki, cfg.kp, ...
                        P.Ns, P.NLanesGard, adaptOpts, P.CfoEnable, ...
                        cfg.Po2Twiddle);
                case 'cd_godard_cma'
                    [eqSym, cfoBinsApplied] = ...
                        eq_clk.combined_cd_fd_godard_adaptive( ...
                        rxSig, P.SpS, P.NFFT, cfg.NOverlap, P.D, P.L_km, ...
                        P.CWL, P.Rs, P.Rolloff, cfg.ki, cfg.kp, ...
                        P.Ns, adaptOpts, P.CfoEnable, cfg.Po2Twiddle);
                otherwise
                    eqSym = [];
            end
        end

        function y = removeKnownCFO(eqSym, CFO_GHz, cfoBinsApplied, P)
            % Combined block output is at the symbol rate.  Remove the
            % residual CFO that the block's one-shot coarse correction
            % could not absorb (estimation error): residual = true CFO
            % minus the continuous-valued estimate the block applied.
            % When the block correction is disabled, cfoBinsApplied = 0
            % and the full CFO is removed here.
            binGHz   = P.SpS * P.Rs / P.NFFT;
            residGHz = CFO_GHz - cfoBinsApplied * binGHz;
            n = (0 : size(eqSym, 1) - 1).';
            phi = 2*pi * (residGHz / P.Rs) * n;
            y = eqSym .* exp(-1j * phi);
        end

        function BER = computeBER(eqSym, symbols, NOut)
            % 32 pi/16 phase rotations + polarisation swap, picking the
            % alignment that minimises BER per polarisation.
            nlen   = size(eqSym, 1);
            refEnd = min(NOut + nlen, size(symbols, 1));
            ref    = symbols(NOut+1 : refEnd, :);
            m      = min(size(eqSym, 1), size(ref, 1));
            if m < 1
                BER = NaN;
                return;
            end
            eqSym = eqSym(1:m, :);
            ref   = ref(1:m, :);

            nPol = size(ref, 2);
            totErr = 0; totBits = 0;
            for p = 1:nPol
                refBits = modem.symbolsToBits(ref(:, p));
                best    = Inf;
                for q = 1:size(eqSym, 2)
                    for kk = 0:31
                        rotated = eqSym(:, q) .* exp(-1j * kk * pi/16);
                        dec     = modem.decideSymbols(rotated);
                        bits    = modem.symbolsToBits(dec);
                        k       = min(numel(refBits), numel(bits));
                        e       = sum(refBits(1:k) ~= bits(1:k)) / k;
                        if e < best, best = e; end
                    end
                end
                totErr  = totErr  + best * numel(refBits);
                totBits = totBits + numel(refBits);
            end
            BER = totErr / totBits;
        end

    end
end
