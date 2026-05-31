classdef combined_eq_clk_fxp_sweep < matlab.unittest.TestCase
%COMBINED_EQ_CLK_FXP_SWEEP  Fixed-point precision sweep for the two
%   combined CD-FD + clock-recovery + adaptive-CMA blocks.
%
%   Each algorithm:
%       cd_gardner_cma : eq_clk.combined_cd_fd_gardner_adaptive_fxp
%       cd_godard_cma  : eq_clk.combined_cd_fd_godard_adaptive_fxp
%   is exercised with the po2-twiddle option both on and off, over a
%   THREE-dimensional fixed-point precision grid that sets the fractional
%   length of each pipeline stage INDEPENDENTLY:
%
%       StaticFL_vec - static equaliser (T.Static): the SINGLE precision of
%                      the FFT/IFFT (input cast, butterflies, output), the
%                      twiddles, the CD response, and the frequency-domain
%                      CD + matched-filter multiply.  The forward FFT
%                      divides by 2 each stage (cumulative 1/N), losing up
%                      to log2(NFFT) LSBs, so this stage usually needs more
%                      fractional bits than the others; its axis is set
%                      separately.
%       AdaptFL_vec  - adaptive equaliser (T.AdaptEq).  In the struct path
%                      of adaptive_eq.equalize_fxp_types this is the
%                      *gradient* precision (T.grad); the data path is
%                      pinned high.
%       ClkFL_vec    - clock-recovery loop (T.Clk for Gardner / T.Godard
%                      for the modified-Godard PI loop).
%
%   The full sweep is the Cartesian product
%       blocks x Po2Twiddle_vec x StaticFL_vec x AdaptFL_vec x ClkFL_vec.
%   Integer bits are fixed at NIntBits; each section forms its struct as
%   struct('WL', NIntBits + FL, 'FL', FL).
%
%   CFO is held at 3 GHz for every configuration; the combined block
%   performs its one-shot coarse CFO correction in floating point
%   internally and the post-block exact-CFO removal absorbs the residual.
%
%   Codegen MEX dispatch
%       Each precision combo bakes the per-section fi numerictypes into a
%       separate MEX (Tcfg is a -args constant at codegen time).  The
%       po2-twiddle flag is a runtime argument and so does NOT spawn its
%       own MEX.  The TestClassSetup loops over
%       (block, StaticFL, AdaptFL, ClkFL) and produces one MEX per combo at
%           src/+eq_clk/<base>_a<AdaptFL>c<ClkFL>s<StaticFL>_mex.mexw64
%       where <base> is the fxp function's name.  The sweep itself
%       dispatches via feval at runtime, so per-config overhead is just
%       one MATLAB->MEX call.
%
%       NOTE: the number of MEX builds is
%           numel(blocks) * numel(StaticFL_vec) * numel(AdaptFL_vec)
%                         * numel(ClkFL_vec)
%       which grows quickly — trim the FL vectors to keep build time sane.
%
%   The user fills in the per-block manual designs (filter lengths and
%   loop-filter gains).  blocks{di} pairs with NCD(di), NTaps(di), and
%   ki/kp(di, po2_idx).
%
%   Saved to combined_eq_clk_fxp_sweep.mat for downstream processing.
%
%   Run with:
%       runtests('combined_eq_clk_fxp_sweep')

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
        PMD_seed = 12345
        L_km     = 80
        SFO_ppm  = 40

        % --- CFO held constant at 3 GHz -----------------------------
        CFO_GHz  = 3.0

        % --- Pulse shaping ------------------------------------------
        Rolloff  = 0.25
        Span     = 10

        % --- Monte-Carlo --------------------------------------------
        Ns          = 18750         % symbols per polarisation per trial
        NTrials     = 15
        SNR_dB_vec  = 0 : 2 : 30

        % --- Static-equaliser FFT size ------------------------------
        NFFT     = 128

        % --- Manual per-block designs --------------------------------
        %  blocks{di} pairs with NCD(di), NTaps(di); ki/kp are
        %  [NDESIGN x NPO2] matrices indexed by (design_row,
        %  Po2Twiddle_vec column).  Fill in from the design-sweep
        %  results.
        blocks         = {'cd_gardner_cma', 'cd_godard_cma'}
        Po2Twiddle_vec = [false, true]
        NCD            = [22,   22]
        NTaps          = [1,    1]
        %  Row 1 = Gardner DPLL gains.  Row 2 = Godard PI gains.  Both rows
        %  are only SEED / fallback values: when TuneGardner / TuneGodard is
        %  true the corresponding fixed-point tuning stage (tuneGardnerGains /
        %  tuneGodardGains) overwrites that row per po2 column with the
        %  (ki,kp) pair that minimises BER on the block's
        %  Ki*_vec x Kp*_vec grid.  Either tuner can be disabled
        %  independently to fall back on the hand-set seeds below.
        ki             = [1e-6      3e-7; ...
                          3e-1   1]
        kp             = [1e-3      1e-3; ...
                          1e-1   3e-3]

        % --- Gain source (sweep vs user-selected) -------------------
        %  Master switch over BOTH loop-filter gain tuners.
        %    false: run the fixed-point gain sweeps below (per-block
        %           TuneGardner / TuneGodard still apply individually).
        %    true : skip ALL tuning and use the hand-set ki/kp values in the
        %           properties above exactly as written (row 1 = Gardner,
        %           row 2 = Godard).  Takes precedence over TuneGardner /
        %           TuneGodard.
        UseManualGains = true

        % --- Gardner fixed-point gain tuning ------------------------
        %  Before the main precision sweep, the Gardner DPLL gains (row 1 of
        %  ki/kp) are retuned IN FIXED POINT exactly like the Godard stage
        %  below: for each po2 setting the Cartesian grid
        %  KiGardner_vec x KpGardner_vec is run through the
        %  precision-specialised Gardner MEX at the representative
        %  (TuneSNR_dB) operating point with EVERY stage held at high
        %  precision (all FLs = TuneFL), so the gains are tuned independently
        %  of the precision grid that is swept afterwards.  The pair giving
        %  the lowest BER is written back into
        %  P.ki(gardner,:) / P.kp(gardner,:).  The grids straddle the
        %  hand-set seeds (po2=off ki=1e-7/kp=1e-3, po2=on ki=1e-5/kp=1e-4),
        %  which already work fairly well.  Set TuneGardner=false to skip and
        %  use the seeds above.
        TuneGardner   = true
        KiGardner_vec = [1e-8 3e-8 1e-7 3e-7 1e-6 3e-6 1e-5 3e-5 1e-4]
        KpGardner_vec = [1e-5 3e-5 1e-4 3e-4 1e-3 3e-3 1e-2 3e-2 1e-1]

        % --- Godard fixed-point gain tuning -------------------------
        %  Before the main precision sweep, the Godard PI gains (row 2 of
        %  ki/kp) are retuned IN FIXED POINT: for each po2 setting the
        %  Cartesian grid KiGodard_vec x KpGodard_vec is run through the
        %  precision-specialised Godard MEX at the representative
        %  (TuneSNR_dB) operating point with EVERY stage held at high
        %  precision (all FLs = TuneFL), so the gains are tuned independently
        %  of the precision grid that is swept afterwards.  The pair giving
        %  the lowest BER is written back into
        %  P.ki(godard,:) / P.kp(godard,:).  Set TuneGodard=false to skip and
        %  use the hand-set seeds above.
        TuneGodard    = true
        KiGodard_vec  = [1e-4 3e-4 1e-3 3e-3 1e-2 3e-2 1e-1 3e-1 1e0]
        KpGodard_vec  = [1e-3 3e-3 1e-2 3e-2 1e-1 3e-1 1e0 3e0 1e1]
        TuneSNR_dB    = 20

        % --- Tuning precision ---------------------------------------
        %  Fractional length applied to ALL pipeline stages (static, adapt,
        %  clk) during BOTH gain-tuning stages.  Held high so the loop gains
        %  are chosen without quantisation noise from the swept precisions; a
        %  dedicated high-precision MEX (a<TuneFL>c<TuneFL>s<TuneFL>) is built
        %  per block for the tuners.
        TuneFL        = 30

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
        NLanesGard  = 32

        % --- Coarse FD CFO correction (floating point inside block) -
        CfoEnable   = true

        % --- Precision sweep (3-D, per-stage fractional lengths) ----
        %  THREE independent swept axes set the fractional length of each
        %  pipeline stage separately.  The full sweep is the Cartesian
        %  product blocks x Po2Twiddle_vec x StaticFL_vec x AdaptFL_vec x
        %  ClkFL_vec.  Integer bits are fixed at NIntBits; each section
        %  forms its struct as struct('WL', NIntBits + FL, 'FL', FL).
        %
        %    StaticFL_vec - T.Static: the single FFT/IFFT precision (input
        %                   cast, butterflies, output), twiddles, CD
        %                   response, FD CD+MF multiply.  Usually higher
        %                   than the other stages because the forward FFT's
        %                   per-stage /2 loses up to log2(NFFT) LSBs.
        %    AdaptFL_vec  - T.AdaptEq (the adaptive-equaliser gradient
        %                   precision; data path is pinned high).
        %    ClkFL_vec    - T.Clk (Gardner) / T.Godard (modified-Godard PI).
        NIntBits     = 16
        StaticFL_vec = [30]
        AdaptFL_vec  = [4,6,8,10,12]
        ClkFL_vec    = [30]

        % --- FEC threshold used to score designs --------------------
        FEC_BER = 2e-2

        % --- MEX build control --------------------------------------
        %  false: reuse any existing per-precision MEX
        %         (*_a<AdaptFL>c<ClkFL>s<StaticFL>_mex) — fast reruns.
        %  true : force a fresh codegen of EVERY combo at startup, ignoring
        %         the cache.  Set this after editing any *_fxp.m source (e.g.
        %         the recovery_fxp / Godard loop-filter changes) so stale
        %         cached binaries are regenerated.
        ForceRebuild = false
    end

    %% ================================================================
    %  Setup
    %% ================================================================
    methods (TestClassSetup)
        function setupAndBuildMex(testCase)
            here     = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(here));
            addpath(genpath(fullfile(repoRoot, 'src')));
            addpath(fullfile(repoRoot, 'build'));

            P = combined_eq_clk_fxp_sweep.extractParams(testCase);

            cfgCoder = coder.config('mex');
            cfgCoder.GenerateReport = false;

            NBlk    = numel(P.blocks);
            NStatic = numel(P.StaticFL_vec);
            NAdapt  = numel(P.AdaptFL_vec);
            NClk    = numel(P.ClkFL_vec);
            nTotal  = NBlk * NStatic * NAdapt * NClk;
            bi = 0;
            tStart = tic;
            if P.ForceRebuild
                fprintf('\n=== Building %d MEX combos (ForceRebuild ON: cache ignored) ===\n', nTotal);
            else
                fprintf('\n=== Building %d MEX combos (using cache where present) ===\n', nTotal);
            end
            for blk = 1:NBlk
                blkName = P.blocks{blk};
                for sfi = 1:NStatic
                    staticFL = P.StaticFL_vec(sfi);
                    for afi = 1:NAdapt
                        adaptFL = P.AdaptFL_vec(afi);
                        for cfi = 1:NClk
                            clkFL = P.ClkFL_vec(cfi);
                            bi = bi + 1;
                            mexBase = combined_eq_clk_fxp_sweep.mexBaseName( ...
                                blkName, adaptFL, clkFL, staticFL);
                            mexPath = fullfile(repoRoot, 'src', '+eq_clk', ...
                                [mexBase '.mexw64']);
                            if isfile(mexPath) && ~P.ForceRebuild
                                fprintf('  [%3d/%3d] %s -> cached\n', ...
                                    bi, nTotal, mexBase);
                                continue;
                            end
                            fprintf('  [%3d/%3d] %s ', bi, nTotal, mexBase);
                            t1 = tic;
                            combined_eq_clk_fxp_sweep.buildOneMex( ...
                                P, cfgCoder, blkName, adaptFL, clkFL, ...
                                staticFL, mexBase, repoRoot);
                            fprintf('(%.0fs)\n', toc(t1));
                        end
                    end
                end
            end

            % --- High-precision tuning MEX -----------------------------
            %  Both gain tuners run at TuneFL on every stage so the loop
            %  gains are chosen free of the swept quantisation.  That combo
            %  (a<TuneFL>c<TuneFL>s<TuneFL>) is not generally part of the grid
            %  above, so build it once per block here (cache-aware).  Skipped
            %  entirely when UseManualGains is set, since no tuning runs.
            tuneFL = P.TuneFL;
            if P.UseManualGains
                fprintf('=== Skipping tuning MEX (UseManualGains=true) ===\n');
            else
                fprintf('=== Building high-precision tuning MEX (all FLs = %d) ===\n', tuneFL);
                for blk = 1:NBlk
                    blkName = P.blocks{blk};
                    mexBase = combined_eq_clk_fxp_sweep.mexBaseName( ...
                        blkName, tuneFL, tuneFL, tuneFL);
                    mexPath = fullfile(repoRoot, 'src', '+eq_clk', ...
                        [mexBase '.mexw64']);
                    if isfile(mexPath) && ~P.ForceRebuild
                        fprintf('  [tune] %s -> cached\n', mexBase);
                        continue;
                    end
                    fprintf('  [tune] %s ', mexBase);
                    t1 = tic;
                    combined_eq_clk_fxp_sweep.buildOneMex( ...
                        P, cfgCoder, blkName, tuneFL, tuneFL, tuneFL, ...
                        mexBase, repoRoot);
                    fprintf('(%.0fs)\n', toc(t1));
                end
            end
            fprintf('Total build time: %.1fs\n', toc(tStart));
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

        function test_fxp_precision_sweep(testCase)
            P = combined_eq_clk_fxp_sweep.extractParams(testCase);

            % --- Initial fixed-point tuning of the loop-filter gains ---
            %  Each tuner overwrites its own block's row of P.ki/P.kp in
            %  place (Gardner = row 1, Godard = row 2).  With UseManualGains
            %  the tuners are bypassed entirely and the hand-set ki/kp values
            %  are used as-is.
            if P.UseManualGains
                fprintf(['\n=== Gain tuning SKIPPED (UseManualGains=true) ' ...
                    '— using hand-set ki/kp ===\n']);
            else
                P = combined_eq_clk_fxp_sweep.tuneGardnerGains(P);
                P = combined_eq_clk_fxp_sweep.tuneGodardGains(P);
            end

            NBlk    = numel(P.blocks);
            NPO2    = numel(P.Po2Twiddle_vec);
            NStatic = numel(P.StaticFL_vec);
            NAdapt  = numel(P.AdaptFL_vec);
            NClk    = numel(P.ClkFL_vec);
            NSNR    = numel(P.SNR_dB_vec);
            NCFG    = NBlk * NPO2 * NStatic * NAdapt * NClk;

            % Validate the manual-design arrays
            assert(numel(P.NCD)   == NBlk && ...
                   numel(P.NTaps) == NBlk && ...
                   isequal(size(P.ki), [NBlk, NPO2]) && ...
                   isequal(size(P.kp), [NBlk, NPO2]), ...
                'blocks / NCD / NTaps / ki / kp size mismatch.');

            % --- Build cfg array ---------------------------------------
            cfgs = struct('block',{}, 'NCD',{}, 'NOverlap',{}, ...
                'NTaps',{}, 'Po2Twiddle',{}, 'ki',{}, 'kp',{}, ...
                'StaticFL',{}, 'AdaptFL',{}, 'ClkFL',{}, ...
                'block_idx',{}, 'po2_idx',{}, ...
                'static_idx',{}, 'adapt_idx',{}, 'clk_idx',{});
            ci = 0;
            for bi = 1:NBlk
                nCd = P.NCD(bi);
                nOv = 2 * ceil((nCd - 1) / 2);
                nT  = P.NTaps(bi);
                for pi = 1:NPO2
                    ki_v = P.ki(bi, pi);
                    kp_v = P.kp(bi, pi);
                    for sfi = 1:NStatic
                        staticFL = P.StaticFL_vec(sfi);
                        for afi = 1:NAdapt
                            adaptFL = P.AdaptFL_vec(afi);
                            for cfi = 1:NClk
                                clkFL = P.ClkFL_vec(cfi);
                                ci = ci + 1;
                                cfgs(ci).block      = P.blocks{bi};
                                cfgs(ci).NCD        = nCd;
                                cfgs(ci).NOverlap   = nOv;
                                cfgs(ci).NTaps      = nT;
                                cfgs(ci).Po2Twiddle = logical(P.Po2Twiddle_vec(pi));
                                cfgs(ci).ki         = ki_v;
                                cfgs(ci).kp         = kp_v;
                                cfgs(ci).StaticFL   = staticFL;
                                cfgs(ci).AdaptFL    = adaptFL;
                                cfgs(ci).ClkFL      = clkFL;
                                cfgs(ci).block_idx  = bi;
                                cfgs(ci).po2_idx    = pi;
                                cfgs(ci).static_idx = sfi;
                                cfgs(ci).adapt_idx  = afi;
                                cfgs(ci).clk_idx    = cfi;
                            end
                        end
                    end
                end
            end

            fprintf('\n=== FXP precision sweep (3-D), CFO = %.2f GHz ===\n', P.CFO_GHz);
            fprintf(['Configurations: %d  (blocks=%d, po2=%d, ' ...
                'static=%d, adapt=%d, clk=%d)  IntBits=%d\n'], ...
                NCFG, NBlk, NPO2, NStatic, NAdapt, NClk, P.NIntBits);
            for ci = 1:NCFG
                fprintf(['  %-16s NCD=%2d NTaps=%d po2=%d  ki=%.2e kp=%.2e  ' ...
                    'AEQ FL=%d  CLK FL=%d  STAT FL=%d\n'], ...
                    cfgs(ci).block, cfgs(ci).NCD, cfgs(ci).NTaps, ...
                    cfgs(ci).Po2Twiddle, cfgs(ci).ki, cfgs(ci).kp, ...
                    cfgs(ci).AdaptFL, cfgs(ci).ClkFL, cfgs(ci).StaticFL);
            end

            % --- BER vs SNR per config ---------------------------------
            ber = cell(NCFG, 1);
            for ci = 1:NCFG
                ber{ci} = nan(P.NTrials, NSNR);
            end

            for tr = 1:P.NTrials
                fprintf('--- trial %d / %d ---\n', tr, P.NTrials);
                for si = 1:NSNR
                    snr = P.SNR_dB_vec(si);
                    [rxSig, symbols] = ...
                        combined_eq_clk_fxp_sweep.buildChannel( ...
                            P, snr, tr, P.CFO_GHz);
                    % Plot the normalised pre-equalisation signal once,
                    % for the very first channel realisation.
                    if tr == 1 && si == 1
                        combined_eq_clk_fxp_sweep.plotNormalisedSignal( ...
                            rxSig, snr, P.CFO_GHz);
                    end
                    for ci = 1:NCFG
                        b = combined_eq_clk_fxp_sweep.runOnePoint( ...
                                P, cfgs(ci), rxSig, symbols, P.CFO_GHz);
                        ber{ci}(tr, si) = b;
                    end
                end
            end

            % --- FEC SNR per config -----------------------------------
            fec_snr = nan(NCFG, 1);
            for ci = 1:NCFG
                fec_snr(ci) = combined_eq_clk_fxp_sweep.fecSnrFromBer( ...
                    P.SNR_dB_vec, mean(ber{ci}, 1, 'omitnan'), P.FEC_BER);
            end

            % --- Build output table -----------------------------------
            block_name = strings(NCFG, 1);
            n_cd       = nan(NCFG, 1);
            n_overlap  = nan(NCFG, 1);
            n_aeq      = nan(NCFG, 1);
            po2        = false(NCFG, 1);
            ki_col     = nan(NCFG, 1);
            kp_col     = nan(NCFG, 1);
            adapt_wl   = nan(NCFG, 1);
            adapt_fl   = nan(NCFG, 1);
            clk_wl     = nan(NCFG, 1);
            clk_fl     = nan(NCFG, 1);
            static_wl  = nan(NCFG, 1);
            static_fl  = nan(NCFG, 1);
            block_idx  = nan(NCFG, 1);
            po2_idx    = nan(NCFG, 1);
            static_idx = nan(NCFG, 1);
            adapt_idx  = nan(NCFG, 1);
            clk_idx    = nan(NCFG, 1);
            for ci = 1:NCFG
                cfg = cfgs(ci);
                block_name(ci) = string(cfg.block);
                n_cd(ci)       = cfg.NCD;
                n_overlap(ci)  = cfg.NOverlap;
                n_aeq(ci)      = cfg.NTaps;
                po2(ci)        = cfg.Po2Twiddle;
                ki_col(ci)     = cfg.ki;
                kp_col(ci)     = cfg.kp;
                adapt_fl(ci)   = cfg.AdaptFL;
                adapt_wl(ci)   = P.NIntBits + cfg.AdaptFL;
                clk_fl(ci)     = cfg.ClkFL;
                clk_wl(ci)     = P.NIntBits + cfg.ClkFL;
                static_fl(ci)  = cfg.StaticFL;
                static_wl(ci)  = P.NIntBits + cfg.StaticFL;
                block_idx(ci)  = cfg.block_idx;
                po2_idx(ci)    = cfg.po2_idx;
                static_idx(ci) = cfg.static_idx;
                adapt_idx(ci)  = cfg.adapt_idx;
                clk_idx(ci)    = cfg.clk_idx;
            end
            tbl = table(block_name, n_cd, n_overlap, n_aeq, po2, ...
                ki_col, kp_col, adapt_wl, adapt_fl, clk_wl, clk_fl, ...
                static_wl, static_fl, block_idx, po2_idx, ...
                static_idx, adapt_idx, clk_idx, ber, fec_snr);
            tbl.Properties.VariableNames{'ki_col'} = 'ki';
            tbl.Properties.VariableNames{'kp_col'} = 'kp';

            % --- Save ------------------------------------------------
            S.tbl          = tbl;
            S.SNR_dB_vec   = P.SNR_dB_vec;
            S.CFO_GHz      = P.CFO_GHz;
            S.NIntBits     = P.NIntBits;
            S.StaticFL_vec = P.StaticFL_vec;
            S.AdaptFL_vec  = P.AdaptFL_vec;
            S.ClkFL_vec    = P.ClkFL_vec;
            S.params       = P;

            outFile = fullfile(fileparts(mfilename('fullpath')), ...
                'combined_eq_clk_fxp_sweep.mat');
            save(outFile, '-struct', 'S');
            fprintf('Saved fxp precision sweep to %s\n', outFile);

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
            P.CFO_GHz       = tc.CFO_GHz;
            P.SFO_ppm       = tc.SFO_ppm;
            P.Rolloff       = tc.Rolloff;
            P.Span          = tc.Span;
            P.Ns            = tc.Ns;
            P.NTrials       = tc.NTrials;
            P.SNR_dB_vec    = tc.SNR_dB_vec;
            P.NFFT          = tc.NFFT;
            P.blocks        = tc.blocks;
            P.Po2Twiddle_vec = tc.Po2Twiddle_vec;
            P.NCD           = tc.NCD;
            P.NTaps         = tc.NTaps;
            P.ki            = tc.ki;
            P.kp            = tc.kp;
            P.UseManualGains = tc.UseManualGains;
            P.TuneGardner   = tc.TuneGardner;
            P.KiGardner_vec = tc.KiGardner_vec;
            P.KpGardner_vec = tc.KpGardner_vec;
            P.TuneGodard    = tc.TuneGodard;
            P.KiGodard_vec  = tc.KiGodard_vec;
            P.KpGodard_vec  = tc.KpGodard_vec;
            P.TuneSNR_dB    = tc.TuneSNR_dB;
            P.TuneFL        = tc.TuneFL;
            P.MuGardner     = tc.MuGardner;
            P.N1Gardner     = tc.N1Gardner;
            P.MuGodard      = tc.MuGodard;
            P.N1Godard      = tc.N1Godard;
            P.NOut          = tc.NOut;
            P.SignOnly      = tc.SignOnly;
            P.SingleSpike   = tc.SingleSpike;
            P.PLanesAEQ     = tc.PLanesAEQ;
            P.NLanesGard    = tc.NLanesGard;
            P.CfoEnable     = tc.CfoEnable;
            P.NIntBits      = tc.NIntBits;
            P.StaticFL_vec  = tc.StaticFL_vec;
            P.AdaptFL_vec   = tc.AdaptFL_vec;
            P.ClkFL_vec     = tc.ClkFL_vec;
            P.FEC_BER       = tc.FEC_BER;
            P.ForceRebuild  = tc.ForceRebuild;
        end

        function P = tuneGardnerGains(P)
            % Initial fixed-point tuning of the Gardner DPLL loop gains.
            %  Mirrors tuneGodardGains: for each po2 setting, sweep the
            %  KiGardner_vec x KpGardner_vec grid through the
            %  precision-specialised Gardner MEX at a single representative
            %  operating point (TuneSNR_dB with every stage at TuneFL) and
            %  pick the (ki,kp) pair with the lowest BER.  The winning pair is
            %  written back into P.ki / P.kp for the Gardner block row; the
            %  seeds are kept if no grid point beats them (or the grid is all
            %  NaN).  Godard is not tuned here.
            bi = find(strcmp(P.blocks, 'cd_gardner_cma'), 1);
            if isempty(bi) || ~P.TuneGardner
                if ~isempty(bi)
                    fprintf(['\n=== Gardner gain tuning SKIPPED ' ...
                        '(TuneGardner=false) — using seed ki/kp ===\n']);
                end
                return;
            end

            staticFL = P.TuneFL;
            adaptFL  = P.TuneFL;
            clkFL    = P.TuneFL;

            nCd = P.NCD(bi);
            nOv = 2 * ceil((nCd - 1) / 2);

            NKi  = numel(P.KiGardner_vec);
            NKp  = numel(P.KpGardner_vec);
            NPO2 = numel(P.Po2Twiddle_vec);

            fprintf(['\n=== Gardner DPLL gain tuning (fixed point) ===\n' ...
                'SNR=%g dB  STAT FL=%d  AEQ FL=%d  CLK FL=%d  ' ...
                'grid=%dx%d  po2 cols=%d\n'], ...
                P.TuneSNR_dB, staticFL, adaptFL, clkFL, NKi, NKp, NPO2);

            % Single channel realisation for the whole tuning grid so the
            % comparison is apples-to-apples.
            [rxSig, symbols] = combined_eq_clk_fxp_sweep.buildChannel( ...
                P, P.TuneSNR_dB, 1, P.CFO_GHz);

            for pj = 1:NPO2
                po2 = logical(P.Po2Twiddle_vec(pj));

                bestBer = Inf;
                bestKi  = P.ki(bi, pj);
                bestKp  = P.kp(bi, pj);

                for ii = 1:NKi
                    for jj = 1:NKp
                        cfg = struct( ...
                            'block',      'cd_gardner_cma', ...
                            'NCD',        nCd, ...
                            'NOverlap',   nOv, ...
                            'NTaps',      P.NTaps(bi), ...
                            'Po2Twiddle', po2, ...
                            'ki',         P.KiGardner_vec(ii), ...
                            'kp',         P.KpGardner_vec(jj), ...
                            'StaticFL',   staticFL, ...
                            'AdaptFL',    adaptFL, ...
                            'ClkFL',      clkFL);
                        ber = combined_eq_clk_fxp_sweep.runOnePoint( ...
                            P, cfg, rxSig, symbols, P.CFO_GHz);
                        fprintf('  po2=%d ki=%.2e kp=%.2e -> BER=%.3e\n', ...
                            po2, cfg.ki, cfg.kp, ber);
                        if isfinite(ber) && ber < bestBer
                            bestBer = ber;
                            bestKi  = cfg.ki;
                            bestKp  = cfg.kp;
                        end
                    end
                end

                P.ki(bi, pj) = bestKi;
                P.kp(bi, pj) = bestKp;
                if isfinite(bestBer)
                    fprintf(['  -> po2=%d BEST ki=%.2e kp=%.2e ' ...
                        '(BER=%.3e)\n'], po2, bestKi, bestKp, bestBer);
                else
                    fprintf(['  -> po2=%d no valid grid point; ' ...
                        'keeping seed ki=%.2e kp=%.2e\n'], ...
                        po2, bestKi, bestKp);
                end
            end
            fprintf('=== Gardner gain tuning complete ===\n');
        end

        function P = tuneGodardGains(P)
            % Initial fixed-point tuning of the Godard PI loop gains.
            %  For each po2 setting, sweep the KiGodard_vec x KpGodard_vec
            %  grid through the precision-specialised Godard MEX at a single
            %  representative operating point (TuneSNR_dB with every stage at
            %  TuneFL) and pick the (ki,kp) pair with the lowest BER.  The
            %  winning pair is written back into P.ki / P.kp for the Godard
            %  block row; the seeds are kept if no grid point beats them (or
            %  the grid is all NaN).  Gardner is not tuned.
            bi = find(strcmp(P.blocks, 'cd_godard_cma'), 1);
            if isempty(bi) || ~P.TuneGodard
                if ~isempty(bi)
                    fprintf(['\n=== Godard gain tuning SKIPPED ' ...
                        '(TuneGodard=false) — using seed ki/kp ===\n']);
                end
                return;
            end

            staticFL = P.TuneFL;
            adaptFL  = P.TuneFL;
            clkFL    = P.TuneFL;

            nCd = P.NCD(bi);
            nOv = 2 * ceil((nCd - 1) / 2);

            NKi  = numel(P.KiGodard_vec);
            NKp  = numel(P.KpGodard_vec);
            NPO2 = numel(P.Po2Twiddle_vec);

            fprintf(['\n=== Godard PI gain tuning (fixed point) ===\n' ...
                'SNR=%g dB  STAT FL=%d  AEQ FL=%d  CLK FL=%d  ' ...
                'grid=%dx%d  po2 cols=%d\n'], ...
                P.TuneSNR_dB, staticFL, adaptFL, clkFL, NKi, NKp, NPO2);

            % Single channel realisation for the whole tuning grid so the
            % comparison is apples-to-apples.
            [rxSig, symbols] = combined_eq_clk_fxp_sweep.buildChannel( ...
                P, P.TuneSNR_dB, 1, P.CFO_GHz);

            for pj = 1:NPO2
                po2 = logical(P.Po2Twiddle_vec(pj));

                bestBer = Inf;
                bestKi  = P.ki(bi, pj);
                bestKp  = P.kp(bi, pj);

                for ii = 1:NKi
                    for jj = 1:NKp
                        cfg = struct( ...
                            'block',      'cd_godard_cma', ...
                            'NCD',        nCd, ...
                            'NOverlap',   nOv, ...
                            'NTaps',      P.NTaps(bi), ...
                            'Po2Twiddle', po2, ...
                            'ki',         P.KiGodard_vec(ii), ...
                            'kp',         P.KpGodard_vec(jj), ...
                            'StaticFL',   staticFL, ...
                            'AdaptFL',    adaptFL, ...
                            'ClkFL',      clkFL);
                        ber = combined_eq_clk_fxp_sweep.runOnePoint( ...
                            P, cfg, rxSig, symbols, P.CFO_GHz);
                        fprintf('  po2=%d ki=%.2e kp=%.2e -> BER=%.3e\n', ...
                            po2, cfg.ki, cfg.kp, ber);
                        if isfinite(ber) && ber < bestBer
                            bestBer = ber;
                            bestKi  = cfg.ki;
                            bestKp  = cfg.kp;
                        end
                    end
                end

                P.ki(bi, pj) = bestKi;
                P.kp(bi, pj) = bestKp;
                if isfinite(bestBer)
                    fprintf(['  -> po2=%d BEST ki=%.2e kp=%.2e ' ...
                        '(BER=%.3e)\n'], po2, bestKi, bestKp, bestBer);
                else
                    fprintf(['  -> po2=%d no valid grid point; ' ...
                        'keeping seed ki=%.2e kp=%.2e\n'], ...
                        po2, bestKi, bestKp);
                end
            end
            fprintf('=== Godard gain tuning complete ===\n');
        end

        function [rxSig, symbols] = buildChannel(P, SNR_dB, trialSeed, cfo)
            % Deterministic per-(trial, SNR, CFO) realisation.
            rng(1000 * trialSeed + round(SNR_dB) + 13);

            symbols = (2*randi([0 1], P.Ns, P.N_pol) - 1) ...
                + 1j * (2*randi([0 1], P.Ns, P.N_pol) - 1);

            txSig = modem.rrcPulse(symbols, P.SpS, P.Rolloff, P.Span);

            rxSig = channel.add_chromatic_dispersion(txSig, P.L_km, ...
                P.SpS, P.Rs, P.D, P.CWL);

            % PMD: isolated rng draw so it stays fixed across trials.
            rngState = rng;
            rng(P.PMD_seed);
            rxSig = channel.add_pmd(rxSig, P.L_km, P.SpS, P.Rs, ...
                P.DGDSpec, P.N_pmd);
            rng(rngState);

            rxSig = channel.lo_freq_shift(rxSig, cfo * 1000, ...
                P.Rs, P.SpS);
            rxSig = channel.apply_timing_error(rxSig, P.SFO_ppm, ...
                0, P.SpS);
            rxSig = channel.add_awgn(rxSig, SNR_dB);

            % Normalise into the unit box before the receiver, using the
            % 95th-percentile magnitude as the per-pol scale reference.
            rxSig = modem.normalise(rxSig, 99.9);
        end

        function plotNormalisedSignal(rxSig, SNR_dB, cfo)
            % Scatter of the normalised, pre-equalisation signal (the
            % receiver input after modem.normalise) for each polarisation.
            % Samples lie in the unit box [-1,1] on each axis.
            nPol = size(rxSig, 2);
            fig  = figure('Name', ...
                'Normalised pre-eq signal (first run)', ...
                'Position', [80 80 360 * nPol 360]);
            for p = 1:nPol
                ax = subplot(1, nPol, p, 'Parent', fig);
                plot(ax, real(rxSig(:, p)), imag(rxSig(:, p)), '.', ...
                    'MarkerSize', 2);
                axis(ax, 'equal'); grid(ax, 'on'); box(ax, 'on');
                xlim(ax, [-1.05, 1.05]); ylim(ax, [-1.05, 1.05]);
                xlabel(ax, 'In-phase'); ylabel(ax, 'Quadrature');
                title(ax, sprintf('Pol %d', p));
            end
            sgtitle(fig, sprintf(['Normalised pre-eq signal ', ...
                '(SNR = %g dB, CFO = %.2f GHz)'], SNR_dB, cfo));

            outFile = fullfile(fileparts(mfilename('fullpath')), ...
                'combined_eq_clk_fxp_normalised_input.png');
            exportgraphics(fig, outFile, 'Resolution', 200);
            fprintf('Saved normalised pre-eq signal plot to %s\n', outFile);
        end

        function ber = runOnePoint(P, cfg, rxSig, symbols, cfo)
            try
                [eqSym, cfoBinsApplied] = ...
                    combined_eq_clk_fxp_sweep.runFxpBlock(P, cfg, rxSig);
            catch ME
                fprintf(['    %s NCD=%d po2=%d  AEQ FL=%d CLK FL=%d Stat FL=%d ', ...
                    '-> ERROR: %s\n'], ...
                    cfg.block, cfg.NCD, cfg.Po2Twiddle, ...
                    cfg.AdaptFL, cfg.ClkFL, cfg.StaticFL, ME.message);
                ber = NaN;
                return;
            end
            if isempty(eqSym) || ~all(isfinite(eqSym(:)))
                ber = NaN;
                return;
            end
            eqSym = combined_eq_clk_fxp_sweep.removeKnownCFO( ...
                eqSym, cfo, cfoBinsApplied, P);
            ber = combined_eq_clk_fxp_sweep.computeBER( ...
                eqSym, symbols, P.NOut);
        end

        function [eqSym, cfoBinsApplied] = runFxpBlock(P, cfg, rxSig)
            % Per-block mu / n1 selection.
            switch cfg.block
                case 'cd_gardner_cma'
                    muVal = P.MuGardner;
                    n1Val = P.N1Gardner;
                case 'cd_godard_cma'
                    muVal = P.MuGodard;
                    n1Val = P.N1Godard;
                otherwise
                    error('combined_eq_clk_fxp_sweep:badBlock', ...
                          'Unknown block: %s', cfg.block);
            end

            % Per-stage precisions (independent fractional lengths).
            % Integer bits are fixed at P.NIntBits.
            StaticPrec = struct('WL', P.NIntBits + cfg.StaticFL, 'FL', cfg.StaticFL);
            AdaptPrec  = struct('WL', P.NIntBits + cfg.AdaptFL,  'FL', cfg.AdaptFL);
            ClkPrec    = struct('WL', P.NIntBits + cfg.ClkFL,    'FL', cfg.ClkFL);

            % Build T at runtime: cheap, and needed to know T.Static.x
            % (input cast) and T.AdaptEq.y (Pilots cast) for the MEX.
            switch cfg.block
                case 'cd_gardner_cma'
                    Tcfg = struct( ...
                        'Static',  StaticPrec, ...
                        'Clk',     ClkPrec, ...
                        'AdaptEq', AdaptPrec);
                    T = eq_clk.combined_cd_fd_gardner_adaptive_fxp_types(Tcfg);
                case 'cd_godard_cma'
                    Tcfg = struct( ...
                        'Static',  StaticPrec, ...
                        'Godard',  ClkPrec, ...
                        'AdaptEq', AdaptPrec);
                    T = eq_clk.combined_cd_fd_godard_adaptive_fxp_types(Tcfg);
                otherwise
                    error('combined_eq_clk_fxp_sweep:badBlock', ...
                          'Unknown block: %s', cfg.block);
            end

            % Pilots must be a fi at T.AdaptEq.y precision (the MEX's
            % typed prototype is variable-row [Inf, 2]).
            PilotsEmpty = cast(complex(zeros(0, 2)), 'like', T.AdaptEq.y);

            adaptOpts = struct( ...
                'NTaps',          double(cfg.NTaps), ...
                'Mu',             double(muVal), ...
                'SingleSpike',    logical(P.SingleSpike), ...
                'N1',             double(n1Val), ...
                'NOut',           double(P.NOut), ...
                'SignOnly',       logical(P.SignOnly), ...
                'UpdateStep',     double(1), ...
                'PLanes',         double(P.PLanesAEQ), ...
                'Mode',           double(0), ...
                'Pilots',         PilotsEmpty, ...
                'BlockLen',       double(P.PLanesAEQ), ...
                'SubframeBlocks', double(0));

            rxSig_fi = cast(rxSig, 'like', T.Static.x);

            % Dispatch to the precision-specific MEX built in TestClassSetup.
            mexBase = combined_eq_clk_fxp_sweep.mexBaseName( ...
                cfg.block, cfg.AdaptFL, cfg.ClkFL, cfg.StaticFL);
            mexFn = str2func(['eq_clk.' mexBase]);

            switch cfg.block
                case 'cd_gardner_cma'
                    [yFi, cfoBinsApplied] = mexFn( ...
                        rxSig_fi, double(P.SpS), double(P.NFFT), ...
                        double(cfg.NOverlap), double(P.D), double(P.L_km), ...
                        double(P.CWL), double(P.Rs), double(P.Rolloff), ...
                        double(cfg.ki), double(cfg.kp), double(P.Ns), ...
                        double(P.NLanesGard), adaptOpts, ...
                        logical(P.CfoEnable), logical(cfg.Po2Twiddle), T);

                case 'cd_godard_cma'
                    [yFi, cfoBinsApplied] = mexFn( ...
                        rxSig_fi, double(P.SpS), double(P.NFFT), ...
                        double(cfg.NOverlap), double(P.D), double(P.L_km), ...
                        double(P.CWL), double(P.Rs), double(P.Rolloff), ...
                        double(cfg.ki), double(cfg.kp), double(P.Ns), ...
                        adaptOpts, logical(P.CfoEnable), ...
                        logical(cfg.Po2Twiddle), T);
            end

            eqSym = double(yFi);
        end

        function y = removeKnownCFO(eqSym, CFO_GHz, cfoBinsApplied, P)
            % Combined block output is at the symbol rate.  Remove the
            % residual CFO that the block's coarse correction could not
            % absorb.
            binGHz   = P.SpS * P.Rs / P.NFFT;
            residGHz = CFO_GHz - cfoBinsApplied * binGHz;
            n = (0 : size(eqSym, 1) - 1).';
            phi = 2*pi * (residGHz / P.Rs) * n;
            y = eqSym .* exp(-1j * phi);
        end

        function BER = computeBER(eqSym, symbols, NOut)
            % 32 pi/16 phase rotations + polarisation swap.
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

        function snr = fecSnrFromBer(snrVec, berVec, fecBer)
            % Linear-on-log interpolation to find the SNR at which
            % BER == fecBer.  Returns NaN if the curve does not bracket
            % the FEC threshold.
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
            for k = 1:numel(x) - 1
                if (y(k) - target) * (y(k+1) - target) <= 0
                    snr = x(k) + (target - y(k)) / (y(k+1) - y(k)) * ...
                        (x(k+1) - x(k));
                    return;
                end
            end
            snr = NaN;
        end

        function name = mexBaseName(blkName, adaptFL, clkFL, staticFL)
            % Unique MEX file name per (block, AdaptFL, ClkFL, StaticFL)
            % combo.  The po2-twiddle flag is a runtime argument and is NOT
            % part of the name.
            switch blkName
                case 'cd_gardner_cma'
                    base = 'combined_cd_fd_gardner_adaptive_fxp';
                case 'cd_godard_cma'
                    base = 'combined_cd_fd_godard_adaptive_fxp';
                otherwise
                    error('combined_eq_clk_fxp_sweep:badBlock', ...
                          'Unknown block: %s', blkName);
            end
            name = sprintf('%s_a%02dc%02ds%02d_mex', ...
                base, adaptFL, clkFL, staticFL);
        end

        function buildOneMex(P, cfgCoder, blkName, adaptFL, clkFL, ...
                             staticFL, mexBase, repoRoot)
            % Codegen one precision-specialised MEX into src/+eq_clk/.
            StaticPrec = struct('WL', P.NIntBits + staticFL, 'FL', staticFL);
            AdaptPrec  = struct('WL', P.NIntBits + adaptFL,  'FL', adaptFL);
            ClkPrec    = struct('WL', P.NIntBits + clkFL,    'FL', clkFL);
            outPath = fullfile(repoRoot, 'src', '+eq_clk', mexBase);

            % Pick nominal NCD/NTaps/ki/kp from the first block design;
            % all of these are runtime args, only their *types* (double
            % scalars) are baked in by codegen.
            nCdProto = P.NCD(1);
            nOvProto = 2 * ceil((nCdProto - 1) / 2);

            switch blkName
                case 'cd_gardner_cma'
                    Tcfg = struct('Static',  StaticPrec, ...
                                  'Clk',     ClkPrec, ...
                                  'AdaptEq', AdaptPrec);
                    T = eq_clk.combined_cd_fd_gardner_adaptive_fxp_types(Tcfg);
                    muVal = P.MuGardner;
                    n1Val = P.N1Gardner;
                case 'cd_godard_cma'
                    Tcfg = struct('Static',  StaticPrec, ...
                                  'Godard',  ClkPrec, ...
                                  'AdaptEq', AdaptPrec);
                    T = eq_clk.combined_cd_fd_godard_adaptive_fxp_types(Tcfg);
                    muVal = P.MuGodard;
                    n1Val = P.N1Godard;
            end

            x_proto = fi(complex(0,0), numerictype(T.Static.x), fimath(T.Static.x));
            In_type = coder.typeof(x_proto, [Inf, P.N_pol], [true, false]);

            p_proto     = fi(complex(0,0), numerictype(T.AdaptEq.y), fimath(T.AdaptEq.y));
            Pilots_type = coder.typeof(p_proto, [Inf, 2], [true, false]);

            AdaptOptsProto = struct( ...
                'NTaps',          double(P.NTaps(1)), ...
                'Mu',             double(muVal), ...
                'SingleSpike',    logical(P.SingleSpike), ...
                'N1',             double(n1Val), ...
                'NOut',           double(P.NOut), ...
                'SignOnly',       logical(P.SignOnly), ...
                'UpdateStep',     double(1), ...
                'PLanes',         double(P.PLanesAEQ), ...
                'Mode',           double(0), ...
                'Pilots',         p_proto, ...
                'BlockLen',       double(P.PLanesAEQ), ...
                'SubframeBlocks', double(0));
            AdaptOptsType = coder.typeof(AdaptOptsProto);
            AdaptOptsType.Fields.Pilots = Pilots_type;

            switch blkName
                case 'cd_gardner_cma'
                    args = { ...
                        In_type, ...
                        double(P.SpS), double(P.NFFT), double(nOvProto), ...
                        double(P.D), double(P.L_km), double(P.CWL), ...
                        double(P.Rs), double(P.Rolloff), ...
                        double(P.ki(1,1)), double(P.kp(1,1)), ...
                        double(P.Ns), double(P.NLanesGard), ...
                        AdaptOptsType, ...
                        logical(P.CfoEnable), false, ...
                        T};
                    codegen('-config', cfgCoder, ...
                        'eq_clk.combined_cd_fd_gardner_adaptive_fxp', ...
                        '-args', args, '-o', outPath);
                case 'cd_godard_cma'
                    args = { ...
                        In_type, ...
                        double(P.SpS), double(P.NFFT), double(nOvProto), ...
                        double(P.D), double(P.L_km), double(P.CWL), ...
                        double(P.Rs), double(P.Rolloff), ...
                        double(P.ki(1,1)), double(P.kp(1,1)), ...
                        double(P.Ns), ...
                        AdaptOptsType, ...
                        logical(P.CfoEnable), false, ...
                        T};
                    codegen('-config', cfgCoder, ...
                        'eq_clk.combined_cd_fd_godard_adaptive_fxp', ...
                        '-args', args, '-o', outPath);
            end
        end

    end
end
