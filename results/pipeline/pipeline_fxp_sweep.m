classdef pipeline_fxp_sweep < matlab.unittest.TestCase
%PIPELINE_FXP_SWEEP  Full receiver pipeline BER under fxp precision sweeps.
%
%   Pipeline:
%       CPON Tx -> channel -> normalise -> ADC (ENOB) -> Godard-combined eq ->
%       pilot lag-align to subframe boundary -> FR -> pilots-only CR ->
%       decision -> BER.
%
%   Blocks under test (all fixed-point):
%       eq_clk.combined_cd_fd_godard_adaptive_fxp_old  (ARCHIVED circular-wrap
%                                                    static CD/MF + Godard
%                                                    timing + adaptive CMA)
%       freq_recovery.differential_kay_fxp          (training-aided FR)
%       carrier_recovery.pilots_only_fxp            (pilot-only CR)
%
%   Precision sweep
%     The fixed-point precision is split into two groups that are swept
%     INDEPENDENTLY, plus the power-of-two twiddle flag as its own axis:
%
%        P.Po2_vec        - power-of-two twiddle on (1) / off (0), eq FFTs
%        P.EqClkDesigns   - rows [StaticFL, ClkFL, AdaptFL]  (equalisation +
%                           clock recovery: static CD/MF, Godard timing loop,
%                           adaptive equaliser)
%        P.CarrierDesigns - rows [FRFL, CRFL]                (carrier recovery:
%                           frequency recovery + phase recovery)
%
%     where each *FL is the fractional length of that stage.  Integer bits
%     are held at IntBits, so each stage word length is IntBits + FL.
%
%     The two groups are swept under an OUTER sweep over the ADC effective
%     number of bits, P.ENOB_vec.  The full config list is therefore the
%     Cartesian product
%        ENOB_vec x Po2_vec x rows(EqClkDesigns) x rows(CarrierDesigns),
%     enumerated ENOB-major.  Sweeping the eq/clk and carrier groups
%     independently lets a coarse equaliser pair with a fine carrier
%     recovery and vice-versa.
%
%   CPON pipeline integration
%     The received signal is amplitude-normalised into the unit box
%     (modem.normalise) before the receiver.  The post/pre-normalisation
%     energy ratio is forwarded to the adaptive equaliser as the CMA radius
%     scale
%     (AdaptOpts.RScale).  The equaliser sees pilots / training at the data
%     energy of +/-1+/-1j (the output of modem.modulate).  After the
%     equaliser the stream is locked onto a CPON subframe boundary by a
%     pilot-phase search (alignToSubframe) over (lag x {direct, X/Y-swapped}):
%     the blind butterfly CMA converges with the polarisations in order or
%     swapped, so the search resolves that swap (de-swapping the columns) at
%     the same time as the boundary, then truncates to an integer number of
%     subframes, so the training-aided FR and the carrier-recovery pilots line
%     up.  Pilots / training remain at +/-1+/-1j throughout (no pilot
%     rescaling).
%
%   Codegen MEX dispatch
%     Each distinct precision per block is built into its own MEX up front in
%     TestClassSetup (po2 and ENOB are runtime / channel-side and do NOT
%     spawn MEX variants).  Existing MEX files at the per-precision suffix
%     are skipped unless ForceRebuild.  The sweep dispatches via str2func at
%     runtime, so per-config overhead is just one MATLAB -> MEX call.
%
%   Output
%     pipeline_fxp_sweep.mat with a table whose rows correspond to
%     (enob, design row) and per-trial BERs in a cell column.

    %% ================================================================
    %  Parameters
    %% ================================================================
    properties (Constant)

        % --- System (CPON downstream) ---------------------------------
        Rs       = 30.5
        SpS      = 2
        N_pol    = 2
        D        = 20            % [ps/(nm*km)]
        CWL      = 1550          % [nm]
        DGDSpec  = 0.1           % [ps/sqrt(km)]
        N_pmd    = 1
        PMD_seed = 42
        L_km     = 80
        SFO_ppm  = 40
        CFO_GHz  = 3          % worst-case CPON CFO
        LW_Hz    = 1e6           % 1 MHz combined laser linewidth
        % Phase noise is OFF to match combined_eq_clk_fxp_sweep (which never
        % adds it) — the Godard timing loop is verified to work in that
        % regime.  Set true to re-introduce the 1 MHz linewidth.
        PhaseNoiseOn = true

        % --- Pulse shaping --------------------------------------------
        Rolloff = 0.25
        Span    = 10

        % --- Monte-Carlo ----------------------------------------------
        N_sub_target = 5         % CPON subframes per trial (Tx)
        NTrials      = 10
        SNR_dB_vec   = 0 : 2 : 30
        % --- Equaliser (overlap-save + Godard + adaptive) ------------
        NFFT        = 128
        NCD         = 22

        % Godard PI gains, one pair per po2-twiddle value.  These are SEED /
        % fallback values (row 2 of ki/kp in combined_eq_clk_fxp_sweep.m); the
        % gains that combined sweep uses are tuned for SHORT records with no
        % phase noise and do NOT transfer to the pipeline's long CPON record +
        % 1 MHz phase noise (the Godard timing loop then closes the eye).
        % When TuneGodard is true they are overwritten per po2 by the
        % tuneGodardGains stage below before the main sweep.
        ki_godard_po2_off = 3e-1
        kp_godard_po2_off = 1e-1
        ki_godard_po2_on  = 1
        kp_godard_po2_on  = 3e-3

        % --- Godard PI gain tuning (before the precision sweep) --------
        %  Gardner timing recovers cleanly on this channel but the Godard PI
        %  loop does not with the seed gains, so retune them for the
        %  pipeline's operating point first.  For each po2 setting in
        %  Po2_vec, the grid KiGodard_vec x KpGodard_vec is run through the
        %  FULL pipeline (eq + lag-align + FR + CR) at TuneSNR_dB with every
        %  stage held at TuneFL, and the (ki,kp) giving the lowest BER is
        %  written into ki_godard_po2_off/on.  Seeds are kept if no grid point
        %  is valid or TuneGodard is false.
        %  Disabled by default so the pipeline uses the SAME fixed manual
        %  gains as combined_eq_clk_fxp_sweep (ki=3e-1, kp=1e-1 for po2 off),
        %  i.e. an identical Godard operating point.  Set true to re-tune.
        TuneGodard   = false
        KiGodard_vec = [1e-2 3e-2 1e-1 3e-1 1]
        KpGodard_vec = [1e-2 3e-2 1e-1 3e1 1 3]
        TuneSNR_dB   = 20
        TuneTrial    = 1
        TuneFL       = 24       % high precision so gains are tuned ~lossless

        % Adaptive eq (matches the existing combined_eq_clk_sweep)
        NTapsAEQ    = 1
        MuAEQ       = 1e-3
        N1AEQ       = 500
        NOutAEQ     = 1000
        SignOnly    = true
        SingleSpike = true
        PLanesAEQ   = 32

        % --- Frequency recovery (differential_kay) --------------------
        MaxFreq       = 0.1
        TrainingLen   = 11

        % --- Carrier recovery (pilots_only) ---------------------------
        BlockLen_CR = 32

        % --- CORDIC iterations ----------------------------------------
        %  CordicIts is NOT a fixed constant: the CORDIC angular resolution
        %  inside frequency recovery (differential_kay_fxp) and carrier
        %  recovery (pilots_only_fxp) is tied to that block's swept
        %  fractional length — FR uses CordicIts = FRFL, CR uses
        %  CordicIts = CRFL.  Each per-precision MEX bakes the matching
        %  iteration count (coder.Constant), so no separate parameter.

        % --- Subframe boundary alignment ------------------------------
        %  The combined equaliser's net symbol delay is not guaranteed to be
        %  exactly NOutAEQ (the Godard timing loop can add a constant offset),
        %  so the receiver locks onto the CPON subframe boundary with a
        %  pilot-phase lag search over +/- AlignLagMax symbols around the
        %  nominal drop, before frequency / carrier recovery.
        AlignLagMax = 64

        % --- Coarse FD CFO correction (floating point inside eq) ------
        CfoEnable   = true

        % --- Receiver front-end (normalisation) -----------------------
        %  The received signal is amplitude-normalised into the unit box
        %  (modem.normalise at NormPct) before entering the receiver.  The
        %  post/pre-normalisation energy ratio is forwarded to the adaptive
        %  equaliser as the CMA radius scale (AdaptOpts.RScale).
        NormPct = 99.9           % modem.normalise percentile

        % --- Fixed point ----------------------------------------------
        IntBits = 16

        % --- MEX build control ----------------------------------------
        %  false: reuse any existing per-precision MEX — fast reruns.
        %  true : force a fresh codegen of EVERY fxp MEX at startup,
        %         ignoring the cache.  Set this after editing any *_fxp.m
        %         source so the stale cached binaries are regenerated.
        ForceRebuild = false

        % --- Precision design: two independently-swept groups ----------
        %  Po2_vec is the power-of-two FFT-twiddle flag (0/1) for the eq
        %  FFTs (shared by the Godard timing loop).  EqClkDesigns and
        %  CarrierDesigns are swept as an independent Cartesian product, so
        %  each eq/clk precision tier can pair with either carrier tier.
        %  The low / high tiers match the equalisation and recovery chapter
        %  tables of report/full/full.tex.
        Po2_vec = [0, 1]
        % Equalisation + clock recovery: [StaticFL, ClkFL, AdaptFL].
        EqClkDesigns = [ ...
            8,  6, 4; ...      % low  precision
            10, 8, 6           % high precision
            ]
        % Carrier recovery: [FRFL, CRFL] (frequency + phase recovery).
        CarrierDesigns = [ ...
            10, 4; ...         % low  precision
            12, 6              % high precision
            ]
        % --- Outer ADC ENOB sweep -------------------------------------
        %  channel.adc effective number of bits.  The full config list is
        %  the Cartesian product
        %  ENOB_vec x Po2_vec x rows(EqClkDesigns) x rows(CarrierDesigns).
        ENOB_vec = [4,6,8]

        % --- CPON subframe constant (must match modem.modulate) --------
        SUBFRAME_SYMS = 3712
    end

    %% ================================================================
    %  Setup: build one MEX per distinct per-block precision
    %% ================================================================
    methods (TestClassSetup)
        function setupAndBuildMex(testCase)
            here     = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(here));
            addpath(genpath(fullfile(repoRoot, 'src')));
            addpath(fullfile(repoRoot, 'build'));

            P = pipeline_fxp_sweep.extractParams(testCase);

            % --- Gather distinct precisions per builder ---------------
            %  EQ: unique (StaticFL, ClkFL, AdaptFL) triples from EqClkDesigns.
            %  FR: unique FRFL (CarrierDesigns col 1).  CR: unique CRFL (col 2).
            %  po2 is a runtime arg and ENOB is channel-side, so neither
            %  spawns MEX variants.
            eqTriples = unique(P.EqClkDesigns, 'rows');
            FL_fr     = unique(P.CarrierDesigns(:, 1));
            FL_cr     = unique(P.CarrierDesigns(:, 2));

            cfgCoder = coder.config('mex');
            cfgCoder.GenerateReport = false;

            nTotal = size(eqTriples, 1) + numel(FL_fr) + numel(FL_cr);
            bi = 0; tStart = tic;
            if P.ForceRebuild
                fprintf('\n=== Building %d MEX combos (ForceRebuild ON: cache ignored) ===\n', nTotal);
            else
                fprintf('\n=== Building %d MEX combos (using cache where present) ===\n', nTotal);
            end

            % --- Equaliser builds (Godard combined) -------------------
            for ti = 1:size(eqTriples, 1)
                staticFL = eqTriples(ti, 1);
                clkFL    = eqTriples(ti, 2);
                adaptFL  = eqTriples(ti, 3);
                bi = bi + 1;
                mexBase = pipeline_fxp_sweep.mexEqName(staticFL, clkFL, adaptFL);
                mexPath = fullfile(repoRoot, 'src', '+eq_clk', ...
                    [mexBase '.mexw64']);
                if isfile(mexPath) && ~P.ForceRebuild
                    fprintf('  [%2d/%2d] %s -> cached\n', bi, nTotal, mexBase);
                    continue;
                end
                fprintf('  [%2d/%2d] %s ', bi, nTotal, mexBase);
                t1 = tic;
                pipeline_fxp_sweep.buildEqMex( ...
                    P, cfgCoder, staticFL, clkFL, adaptFL, mexBase, repoRoot);
                fprintf('(%.0fs)\n', toc(t1));
            end

            % --- differential_kay FR builds ---------------------------
            for fl = FL_fr(:).'
                bi = bi + 1;
                mexBase = pipeline_fxp_sweep.mexFRName(fl);
                mexPath = fullfile(repoRoot, 'src', '+freq_recovery', ...
                    [mexBase '.mexw64']);
                if isfile(mexPath) && ~P.ForceRebuild
                    fprintf('  [%2d/%2d] %s -> cached\n', bi, nTotal, mexBase);
                    continue;
                end
                fprintf('  [%2d/%2d] %s ', bi, nTotal, mexBase);
                t1 = tic;
                pipeline_fxp_sweep.buildFRMex(P, cfgCoder, fl, mexBase, repoRoot);
                fprintf('(%.0fs)\n', toc(t1));
            end

            % --- pilots_only CR builds --------------------------------
            for fl = FL_cr(:).'
                bi = bi + 1;
                mexBase = pipeline_fxp_sweep.mexCRName(fl);
                mexPath = fullfile(repoRoot, 'src', '+carrier_recovery', ...
                    [mexBase '.mexw64']);
                if isfile(mexPath) && ~P.ForceRebuild
                    fprintf('  [%2d/%2d] %s -> cached\n', bi, nTotal, mexBase);
                    continue;
                end
                fprintf('  [%2d/%2d] %s ', bi, nTotal, mexBase);
                t1 = tic;
                pipeline_fxp_sweep.buildCRMex( ...
                    P, cfgCoder, fl, mexBase, repoRoot);
                fprintf('(%.0fs)\n', toc(t1));
            end

            % --- High-precision tuning MEX (eq + FR + CR at TuneFL) ----
            %  tuneGodardGains runs the FULL pipeline at TuneFL, so build
            %  that precision's eq/FR/CR once here (cache-aware) unless they
            %  already exist from the design grid.  Skipped when tuning is off.
            if P.TuneGodard
                tfl = P.TuneFL;
                fprintf('=== Building Godard tuning MEX (all stages FL=%d) ===\n', tfl);
                tuneBuilds = { ...
                    fullfile(repoRoot,'src','+eq_clk'), ...
                        pipeline_fxp_sweep.mexEqName(tfl, tfl, tfl), 'eq'; ...
                    fullfile(repoRoot,'src','+freq_recovery'), ...
                        pipeline_fxp_sweep.mexFRName(tfl), 'fr'; ...
                    fullfile(repoRoot,'src','+carrier_recovery'), ...
                        pipeline_fxp_sweep.mexCRName(tfl), 'cr'};
                for ti = 1:size(tuneBuilds,1)
                    mexBase = tuneBuilds{ti,2};
                    mexPath = fullfile(tuneBuilds{ti,1}, [mexBase '.mexw64']);
                    if isfile(mexPath) && ~P.ForceRebuild
                        fprintf('  [tune] %s -> cached\n', mexBase);
                        continue;
                    end
                    fprintf('  [tune] %s ', mexBase);
                    t1 = tic;
                    switch tuneBuilds{ti,3}
                        case 'eq'
                            pipeline_fxp_sweep.buildEqMex(P, cfgCoder, tfl, tfl, tfl, mexBase, repoRoot);
                        case 'fr'
                            pipeline_fxp_sweep.buildFRMex(P, cfgCoder, tfl, mexBase, repoRoot);
                        case 'cr'
                            pipeline_fxp_sweep.buildCRMex(P, cfgCoder, tfl, mexBase, repoRoot);
                    end
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
    %  Test
    %% ================================================================
    methods (Test)
        function test_pipeline_fxp_sweep(testCase)
            P = pipeline_fxp_sweep.extractParams(testCase);

            % Retune the Godard PI loop gains for the pipeline's operating
            % point before the precision sweep (writes P.ki/kp_godard_po2_*).
            P = pipeline_fxp_sweep.tuneGodardGains(P);

            cfgs = pipeline_fxp_sweep.buildCfgs(P);
            NCFG = numel(cfgs);
            NSNR = numel(P.SNR_dB_vec);

            fprintf(['\n=== Pipeline fxp sweep: %d configs x %d SNRs ', ...
                     'x %d trials ===\n'], NCFG, NSNR, P.NTrials);

            ber = nan(NCFG, NSNR, P.NTrials);
            for ci = 1:NCFG
                cfg = cfgs(ci);
                fprintf(['[%2d/%2d] ENOB=%d po2=%d  StatFL=%2d ClkFL=%2d ', ...
                    'AdaptFL=%2d FRFL=%2d CRFL=%2d\n'], ...
                    ci, NCFG, cfg.enob, cfg.po2, cfg.staticFL, cfg.clkFL, ...
                    cfg.adaptFL, cfg.frFL, cfg.crFL);
                for si = 1:NSNR
                    snr = P.SNR_dB_vec(si);
                    fprintf('    SNR=%4.1f dB  ', snr);
                    for tr = 1:P.NTrials
                        b = pipeline_fxp_sweep.runOneTrial(P, cfg, snr, tr);
                        ber(ci, si, tr) = b;
                        if isfinite(b)
                            fprintf('%.2e ', b);
                        else
                            fprintf(' NaN  ');
                        end
                    end
                    fprintf('\n');
                end
            end

            % --- Output table -----------------------------------------
            enob        = nan(NCFG, 1);
            design_idx  = nan(NCFG, 1);
            eqclk_idx   = nan(NCFG, 1);
            carrier_idx = nan(NCFG, 1);
            po2        = false(NCFG, 1);
            static_fl  = nan(NCFG, 1);
            clk_fl     = nan(NCFG, 1);
            adapt_fl   = nan(NCFG, 1);
            fr_fl      = nan(NCFG, 1);
            cr_fl      = nan(NCFG, 1);
            static_wl  = nan(NCFG, 1);
            clk_wl     = nan(NCFG, 1);
            adapt_wl   = nan(NCFG, 1);
            fr_wl      = nan(NCFG, 1);
            cr_wl      = nan(NCFG, 1);
            ber_cell   = cell(NCFG, 1);
            ber_mean   = cell(NCFG, 1);
            for ci = 1:NCFG
                cfg = cfgs(ci);
                enob(ci)        = cfg.enob;
                design_idx(ci)  = cfg.design_idx;
                eqclk_idx(ci)   = cfg.eqclk_idx;
                carrier_idx(ci) = cfg.carrier_idx;
                po2(ci)        = cfg.po2;
                static_fl(ci)  = cfg.staticFL;
                clk_fl(ci)     = cfg.clkFL;
                adapt_fl(ci)   = cfg.adaptFL;
                fr_fl(ci)      = cfg.frFL;
                cr_fl(ci)      = cfg.crFL;
                static_wl(ci)  = P.IntBits + cfg.staticFL;
                clk_wl(ci)     = P.IntBits + cfg.clkFL;
                adapt_wl(ci)   = P.IntBits + cfg.adaptFL;
                fr_wl(ci)      = P.IntBits + cfg.frFL;
                cr_wl(ci)      = P.IntBits + cfg.crFL;
                % ber(ci, :, :) is [NSNR x NTrials]
                ber_cell{ci} = squeeze(ber(ci, :, :));
                ber_mean{ci} = mean(squeeze(ber(ci, :, :)), 2, 'omitnan').';
            end
            tbl = table(enob, design_idx, eqclk_idx, carrier_idx, po2, ...
                static_fl, clk_fl, adapt_fl, ...
                fr_fl, cr_fl, static_wl, clk_wl, adapt_wl, fr_wl, cr_wl, ...
                ber_cell, ber_mean);
            tbl.Properties.VariableNames{'ber_cell'} = 'ber';
            tbl.Properties.VariableNames{'ber_mean'} = 'ber_mean_per_snr';

            S.tbl        = tbl;
            S.SNR_dB_vec = P.SNR_dB_vec;
            S.CFO_GHz    = P.CFO_GHz;
            S.NTrials    = P.NTrials;
            S.ENOB_vec       = P.ENOB_vec;
            S.Po2_vec        = P.Po2_vec;
            S.EqClkDesigns   = P.EqClkDesigns;
            S.CarrierDesigns = P.CarrierDesigns;
            S.params     = P;
            outFile = fullfile(fileparts(mfilename('fullpath')), ...
                'pipeline_fxp_sweep.mat');
            save(outFile, '-struct', 'S');
            fprintf('Saved to %s\n', outFile);

            testCase.verifyEqual(height(tbl), NCFG);
        end
    end

    %% ================================================================
    %  Static helpers
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
            P.SFO_ppm       = tc.SFO_ppm;
            P.CFO_GHz       = tc.CFO_GHz;
            P.LW_Hz         = tc.LW_Hz;
            P.PhaseNoiseOn  = tc.PhaseNoiseOn;
            P.Rolloff       = tc.Rolloff;
            P.Span          = tc.Span;
            P.N_sub_target  = tc.N_sub_target;
            P.NTrials       = tc.NTrials;
            P.SNR_dB_vec    = tc.SNR_dB_vec;
            P.NFFT          = tc.NFFT;
            P.NCD           = tc.NCD;
            P.ki_godard_po2_off = tc.ki_godard_po2_off;
            P.kp_godard_po2_off = tc.kp_godard_po2_off;
            P.ki_godard_po2_on  = tc.ki_godard_po2_on;
            P.kp_godard_po2_on  = tc.kp_godard_po2_on;
            P.TuneGodard    = tc.TuneGodard;
            P.KiGodard_vec  = tc.KiGodard_vec;
            P.KpGodard_vec  = tc.KpGodard_vec;
            P.TuneSNR_dB    = tc.TuneSNR_dB;
            P.TuneTrial     = tc.TuneTrial;
            P.TuneFL        = tc.TuneFL;
            P.NTapsAEQ      = tc.NTapsAEQ;
            P.MuAEQ         = tc.MuAEQ;
            P.N1AEQ         = tc.N1AEQ;
            P.NOutAEQ       = tc.NOutAEQ;
            P.SignOnly      = tc.SignOnly;
            P.SingleSpike   = tc.SingleSpike;
            P.PLanesAEQ     = tc.PLanesAEQ;
            P.MaxFreq       = tc.MaxFreq;
            P.TrainingLen   = tc.TrainingLen;
            P.BlockLen_CR   = tc.BlockLen_CR;
            P.AlignLagMax   = tc.AlignLagMax;
            P.CfoEnable     = tc.CfoEnable;
            P.NormPct       = tc.NormPct;
            P.IntBits       = tc.IntBits;
            P.ForceRebuild  = tc.ForceRebuild;
            P.Po2_vec        = tc.Po2_vec;
            P.EqClkDesigns   = tc.EqClkDesigns;
            P.CarrierDesigns = tc.CarrierDesigns;
            P.ENOB_vec      = tc.ENOB_vec;
            P.SUBFRAME_SYMS = tc.SUBFRAME_SYMS;
        end

        function cfgs = buildCfgs(P)
            % Cartesian product
            %   ENOB_vec x Po2_vec x EqClkDesigns x CarrierDesigns,
            % ENOB-major (the outer ADC sweep).  design_idx is a flat index
            % over the (po2, eqclk, carrier) combinations within one ENOB;
            % eqclk_idx / carrier_idx identify the group rows.
            cfgs = struct('enob', {}, 'design_idx', {}, ...
                'eqclk_idx', {}, 'carrier_idx', {}, 'po2', {}, ...
                'staticFL', {}, 'clkFL', {}, 'adaptFL', {}, ...
                'frFL', {}, 'crFL', {});
            for ei = 1:numel(P.ENOB_vec)
                di = 0;
                for pi = 1:numel(P.Po2_vec)
                    for qi = 1:size(P.EqClkDesigns, 1)
                        for ri = 1:size(P.CarrierDesigns, 1)
                            di = di + 1;
                            eqRow = P.EqClkDesigns(qi, :);
                            crRow = P.CarrierDesigns(ri, :);
                            cfgs(end + 1).enob = P.ENOB_vec(ei);  %#ok<AGROW>
                            cfgs(end).design_idx  = di;
                            cfgs(end).eqclk_idx   = qi;
                            cfgs(end).carrier_idx = ri;
                            cfgs(end).po2      = logical(P.Po2_vec(pi));
                            cfgs(end).staticFL = eqRow(1);
                            cfgs(end).clkFL    = eqRow(2);
                            cfgs(end).adaptFL  = eqRow(3);
                            cfgs(end).frFL     = crRow(1);
                            cfgs(end).crFL     = crRow(2);
                        end
                    end
                end
            end
        end

        function P = tuneGodardGains(P)
            % Retune the Godard PI loop gains for the pipeline before the
            % precision sweep.  For each po2 setting in Po2_vec, run
            % the KiGodard_vec x KpGodard_vec grid through the FULL pipeline
            % (runOneTrial with a cfg.kiSel/kpSel override) at TuneSNR_dB with
            % every stage at TuneFL, and write the lowest-BER (ki,kp) back into
            % P.ki/kp_godard_po2_*.  The deterministic channel (TuneTrial,
            % TuneSNR_dB) is regenerated identically per grid point, so the
            % comparison is apples-to-apples.
            if ~P.TuneGodard
                fprintf(['\n=== Godard gain tuning SKIPPED ' ...
                    '(TuneGodard=false) — using seed ki/kp ===\n']);
                return;
            end

            tfl     = P.TuneFL;
            NKi     = numel(P.KiGodard_vec);
            NKp     = numel(P.KpGodard_vec);
            po2vals = unique(P.Po2_vec(:)).';      % po2 settings in use

            fprintf(['\n=== Godard PI gain tuning (pipeline) ===\n' ...
                'SNR=%g dB  FL=%d  grid=%dx%d  po2=%s\n'], ...
                P.TuneSNR_dB, tfl, NKi, NKp, mat2str(po2vals));

            for po2 = po2vals
                po2L = logical(po2);
                if po2L
                    bestKi = P.ki_godard_po2_on;  bestKp = P.kp_godard_po2_on;
                else
                    bestKi = P.ki_godard_po2_off; bestKp = P.kp_godard_po2_off;
                end
                bestBer = Inf;

                for ii = 1:NKi
                    for jj = 1:NKp
                        cfg = struct( ...
                            'po2',      po2L, ...
                            'staticFL', tfl, 'clkFL', tfl, 'adaptFL', tfl, ...
                            'frFL',     tfl, 'crFL',  tfl, ...
                            'kiSel',    P.KiGodard_vec(ii), ...
                            'kpSel',    P.KpGodard_vec(jj));
                        ber = pipeline_fxp_sweep.runOneTrial( ...
                            P, cfg, P.TuneSNR_dB, P.TuneTrial);
                        fprintf('  po2=%d ki=%.2e kp=%.2e -> BER=%.3e\n', ...
                            po2L, cfg.kiSel, cfg.kpSel, ber);
                        if isfinite(ber) && ber < bestBer
                            bestBer = ber; bestKi = cfg.kiSel; bestKp = cfg.kpSel;
                        end
                    end
                end

                if po2L
                    P.ki_godard_po2_on  = bestKi;  P.kp_godard_po2_on  = bestKp;
                else
                    P.ki_godard_po2_off = bestKi;  P.kp_godard_po2_off = bestKp;
                end
                if isfinite(bestBer)
                    fprintf('  -> po2=%d BEST ki=%.2e kp=%.2e (BER=%.3e)\n', ...
                        po2L, bestKi, bestKp, bestBer);
                else
                    fprintf(['  -> po2=%d no valid grid point; keeping seed ' ...
                        'ki=%.2e kp=%.2e\n'], po2L, bestKi, bestKp);
                end
            end
            fprintf('=== Godard gain tuning complete ===\n');
        end

        % ------------------------------------------------------------------
        %  Per-trial run: channel + eq + drop subframe + FR + CR + BER
        % ------------------------------------------------------------------
        function ber = runOneTrial(P, cfg, SNR_dB, trialSeed)
            try
                % ADC effective bits for this config.  The gain-tuning cfgs
                % carry no enob field -> no ADC quantisation (tune lossless).
                if isfield(cfg, 'enob')
                    enob = cfg.enob;
                else
                    enob = [];
                end
                [rxSig, txSymbols, training, pilotsRef, ~, RScale] = ...
                    pipeline_fxp_sweep.genChannel(P, SNR_dB, trialSeed, enob);

                % --- 1. Equaliser MEX (Godard combined) ---------------
                %  Per-section precisions: Static / Godard / AdaptEq.
                StaticPrec = struct('WL', P.IntBits + cfg.staticFL, 'FL', cfg.staticFL);
                GodardPrec = struct('WL', P.IntBits + cfg.clkFL,    'FL', cfg.clkFL);
                AdaptPrec  = struct('WL', P.IntBits + cfg.adaptFL,  'FL', cfg.adaptFL);
                Tcfg = struct('Static',  StaticPrec, ...
                              'Godard',  GodardPrec, ...
                              'AdaptEq', AdaptPrec);
                T_eq = eq_clk.combined_cd_fd_godard_adaptive_fxp_types(Tcfg);

                rxSig_fi = cast(rxSig, 'like', T_eq.Static.x);
                PilotsEmpty = cast(complex(zeros(0, 2)), 'like', T_eq.AdaptEq.y);

                adaptOpts = struct( ...
                    'NTaps',          double(P.NTapsAEQ), ...
                    'Mu',             double(P.MuAEQ), ...
                    'SingleSpike',    logical(P.SingleSpike), ...
                    'N1',             double(P.N1AEQ), ...
                    'NOut',           double(P.NOutAEQ), ...
                    'SignOnly',       logical(P.SignOnly), ...
                    'UpdateStep',     double(1), ...
                    'PLanes',         double(P.PLanesAEQ), ...
                    'Mode',           double(0), ...
                    'Pilots',         PilotsEmpty, ...
                    'BlockLen',       double(P.PLanesAEQ), ...
                    'SubframeBlocks', double(0), ...
                    'RScale',         double(RScale));

                nOv = 2 * ceil((P.NCD - 1) / 2);
                eqMex = str2func(['eq_clk.' ...
                    pipeline_fxp_sweep.mexEqName(cfg.staticFL, cfg.clkFL, cfg.adaptFL)]);

                % Godard PI gains.  cfg.kiSel/kpSel override (used by the
                % tuneGodardGains grid); otherwise use the per-po2 gains in P.
                if isfield(cfg, 'kiSel')
                    kiSel = cfg.kiSel;
                    kpSel = cfg.kpSel;
                elseif cfg.po2
                    kiSel = P.ki_godard_po2_on;
                    kpSel = P.kp_godard_po2_on;
                else
                    kiSel = P.ki_godard_po2_off;
                    kpSel = P.kp_godard_po2_off;
                end

                % NsymTx = real-subframe count (txSymbols excludes the guard
                % symbols appended in genChannel).  The equaliser keeps only
                % the first NsymTx output symbols, so the guard tail — where
                % the old Godard's cyclic seam lands — is discarded here,
                % before alignment / FR / CR see the stream.
                NsymTx = size(txSymbols, 1);
                [yFi, ~] = eqMex(rxSig_fi, double(P.SpS), ...
                    double(P.NFFT), double(nOv), ...
                    double(P.D), double(P.L_km), double(P.CWL), ...
                    double(P.Rs), double(P.Rolloff), ...
                    double(kiSel), double(kpSel), ...
                    double(NsymTx), adaptOpts, ...
                    logical(P.CfoEnable), logical(cfg.po2), T_eq);
                eqOut = double(yFi);

                % --- 2. Align to a subframe boundary (pilot lag search) --
                %  The nominal drop (one subframe minus the equaliser's NOut
                %  transient) assumes zero net symbol delay, but the Godard
                %  timing loop can add a constant offset.  Search a lag window
                %  around the nominal boundary and keep the start whose
                %  block-start symbols are most phase-coherent with the known
                %  pilots, then truncate to an integer number of subframes.
                %  No rescaling (pilots are at +/-1+/-1j).
                nDropNom = P.SUBFRAME_SYMS - P.NOutAEQ;
                if nDropNom < 0
                    nDropNom = 0;
                end
                [eqOut, nSubUsable] = pipeline_fxp_sweep.alignToSubframe( ...
                    eqOut, nDropNom, pilotsRef, P);
                if nSubUsable < 1
                    ber = NaN;
                    return;
                end

                % --- 3. Frequency recovery (differential_kay) --------
                %  Pilots / training are at the data energy (+/-1+/-1j), so no
                %  pilot rescaling is applied before FR/CR.
                FxpFR = struct('WL', P.IntBits + cfg.frFL, 'FL', cfg.frFL);
                T_fr = freq_recovery.fxp_types(FxpFR);
                eqOut_fi    = cast(eqOut,    'like', T_fr.x);
                training_fi = cast(training, 'like', T_fr.x);

                frMex = str2func(['freq_recovery.' ...
                    pipeline_fxp_sweep.mexFRName(cfg.frFL)]);
                % CORDIC iterations = FR fractional length (matches the
                % constant baked into this precision's MEX).
                [frOut_fi, ~] = frMex(eqOut_fi, training_fi, ...
                    double(P.Rs), double(cfg.frFL), T_fr, ...
                    true, double(0), double(P.MaxFreq));

                % --- 4. Carrier recovery (pilots_only) ---------------
                FxpCR = struct('WL', P.IntBits + cfg.crFL, 'FL', cfg.crFL);
                T_cr = carrier_recovery.fxp_types(FxpCR);
                frOut_for_cr = cast(double(frOut_fi), 'like', T_cr.x);

                pilotsAll = repmat(pilotsRef, nSubUsable, 1);
                pilots_fi = cast(pilotsAll, 'like', T_cr.x);

                crMex = str2func(['carrier_recovery.' ...
                    pipeline_fxp_sweep.mexCRName(cfg.crFL)]);

                % CORDIC iterations = CR fractional length (matches the
                % constant baked into this precision's MEX).
                [crOut_fi, ~] = crMex(frOut_for_cr, double(P.N_pol), ...
                    double(P.BlockLen_CR), pilots_fi, ...
                    double(cfg.crFL), T_cr);
                crOut = double(crOut_fi);

                % --- 5. BER vs reference (subframes 2..end of Tx) ----
                txOffset = P.SUBFRAME_SYMS;       % matches the discard
                refEnd   = min(txOffset + size(crOut, 1), size(txSymbols, 1));
                refSyms  = txSymbols(txOffset + 1 : refEnd, :);
                m        = min(size(crOut, 1), size(refSyms, 1));
                if m < 1
                    ber = NaN;
                    return;
                end
                ber = pipeline_fxp_sweep.computeBER(crOut(1:m, :), ...
                                                    refSyms(1:m, :));
            catch ME
                fprintf('\n    ERROR (trial %d): %s\n', trialSeed, ME.message);
                ber = NaN;
            end
        end

        % ------------------------------------------------------------------
        %  Channel: CPON modulator (+ guard symbols for the old circular-wrap
        %  Godard) + RRC + CD + PMD + CFO + SFO + AWGN (+ optional phase
        %  noise, off by default to match combined_eq_clk_fxp_sweep), then
        %  normalise + ADC quantisation to ENOB bits (the equaliser sees
        %  pilots/training at +/-1+/-1j).  The returned `symbols` / `nSubTx`
        %  are the REAL subframes only; the guard symbols live only in rxSig
        %  and are trimmed by the equaliser.
        % ------------------------------------------------------------------
        function [rxSig, symbols, training, pilotsRef, nSubTx, RScale] = ...
                genChannel(P, SNR_dB, trialSeed, ENOB)
            if nargin < 4
                ENOB = [];
            end
            rng(1000 * trialSeed + round(SNR_dB) + 13);

            CPON_BITS_PER_SF = 3586 * P.N_pol * 2;
            nBits = P.N_sub_target * CPON_BITS_PER_SF;
            bits  = randi([0 1], nBits, 1);

            [symbols, pilotsRef, training, nSubTx] = modem.modulate(bits);

            % --- Guard symbols for the OLD circular-wrap Godard ----------
            %  The archived Godard cyclically wraps the record (tail glued to
            %  head); the seam is timing-consistent only when the total SFO
            %  drift across the record is an INTEGER number of samples.  Pad
            %  the whole-subframe stream with random QPSK guard symbols so the
            %  record reaches the next integer-drift length.  The seam's
            %  timing discontinuity then falls on these throwaway TAIL symbols,
            %  which the equaliser's NSymb trim (NsymTx = real-subframe count
            %  in runOneTrial) drops before FR / CR.  The returned `symbols`
            %  (BER reference) and `nSubTx` cover the real subframes only.
            nGuard = pipeline_fxp_sweep.guardSyms(P, size(symbols, 1));
            if nGuard > 0
                guard  = (2*randi([0 1], nGuard, P.N_pol) - 1) ...
                       + 1j*(2*randi([0 1], nGuard, P.N_pol) - 1);
                txSyms = [symbols; guard];
            else
                txSyms = symbols;
            end

            txSig = modem.rrcPulse(txSyms, P.SpS, P.Rolloff, P.Span);
            rxSig = channel.add_chromatic_dispersion(txSig, P.L_km, ...
                P.SpS, P.Rs, P.D, P.CWL);

            % PMD held constant by isolating its RNG draw.
            rngState = rng;
            rng(P.PMD_seed);
            rxSig = channel.add_pmd(rxSig, P.L_km, P.SpS, P.Rs, ...
                P.DGDSpec, P.N_pmd);
            rng(rngState);

            % CFO (GHz -> MHz for channel.lo_freq_shift) + SFO + AWGN +
            % laser phase noise.
            rxSig = channel.lo_freq_shift(rxSig, P.CFO_GHz * 1000, ...
                P.Rs, P.SpS);
            rxSig = channel.apply_timing_error(rxSig, P.SFO_ppm, 0, P.SpS);
            rxSig = channel.add_awgn(rxSig, SNR_dB);
            if P.PhaseNoiseOn
                rxSig = channel.add_phase_noise(rxSig, P.Rs * P.SpS, P.LW_Hz);
            end

            % --- Receiver front-end: normalise --------------------------
            %  Energy before/after normalisation gives the CMA radius scale
            %  (RScale) forwarded to the adaptive equaliser.
            eBefore = sum(abs(rxSig(:)).^2);
            rxSig   = modem.normalise(rxSig, P.NormPct);
            eAfter  = sum(abs(rxSig(:)).^2);
            if eBefore > 0
                RScale = eAfter / eBefore;
            else
                RScale = 1;
            end

            % --- ADC: quantise to ENOB effective bits -------------------
            %  Applied AFTER normalisation so the signal fills the ADC
            %  full-scale [-1, 1] and the ENOB bits span that range — this is
            %  what makes the outer ENOB sweep meaningful.  channel.adc
            %  quantises the I and Q of each polarisation independently.
            %  ENOB empty (tuning path) -> no quantisation.
            if ~isempty(ENOB) && isfinite(ENOB)
                rxSig = channel.adc(rxSig, ENOB);
            end
        end

        function nGuard = guardSyms(P, nSubSyms)
            % Number of guard symbols to append so the TOTAL record produces
            % an integer number of samples of SFO drift — the old circular-
            % wrap Godard's convergence condition.  The accumulated drift is
            %   drift_samples = nSyms * SpS * SFO_ppm * 1e-6,
            % so the record length must be a multiple of
            %   symsPerSample = 1 / (SpS * SFO_ppm * 1e-6).
            % Pad up to the next such multiple.  (No drift -> no guard.)
            if P.SFO_ppm == 0
                nGuard = 0;
                return;
            end
            symsPerSample = 1 / (P.SpS * P.SFO_ppm * 1e-6);
            totalSyms     = ceil(nSubSyms / symsPerSample) * symsPerSample;
            nGuard        = round(totalSyms - nSubSyms);
        end

        function [eqOut, nSubUsable, swapped] = alignToSubframe(eqOutFull, nDropNom, pilotsRef, P)
            %ALIGNTOSUBFRAME  Lock the equaliser output onto a CPON subframe
            %   boundary AND resolve the butterfly-CMA polarisation-swap
            %   ambiguity, then return an integer number of subframes.
            %
            %   The block-start symbols of a correctly aligned stream are the
            %   known pilots; multiplying by the conjugate pilot leaves the
            %   (slowly varying) carrier phase, so the block-to-block product
            %   z_b*conj(z_{b-1}) is phase-coherent.  On a wrong lag the
            %   block-start symbols are random data and that product averages
            %   to zero.  Pilots are at +/-1+/-1j (same amplitude as data), so
            %   magnitude cannot localise them — only the phase coherence can.
            %
            %   The butterfly CMA is blind and converges with the two
            %   polarisations either in order or SWAPPED; nothing upstream
            %   resolves this, and a swap decoheres the per-pol pilot phase
            %   (each block gets a different conj(pilotX)*pilotY phasor), so
            %   it floors the BER.  The search therefore runs over
            %   (lag in [-AlignLagMax,+AlignLagMax]) x ({direct, X/Y-swapped}
            %   pilot columns) and keeps the highest-coherence combination.
            %   When the swapped order wins, the output columns are swapped
            %   back so downstream training-aided FR and pilot CR (fixed X/Y
            %   order) line up.
            N = size(eqOutFull, 1);
            orders = {[1 2], [2 1]};      % direct, X/Y-swapped pilot columns
            bestC = -inf; bestLag = 0; bestSwap = false;
            for so = 1:numel(orders)
                pr = pilotsRef(:, orders{so});
                for lag = -P.AlignLagMax : P.AlignLagMax
                    start = nDropNom + 1 + lag;
                    if start < 1, continue; end
                    nSub = floor((N - start + 1) / P.SUBFRAME_SYMS);
                    if nSub < 1, continue; end
                    C = pipeline_fxp_sweep.pilotCoherence( ...
                        eqOutFull, start, nSub, pr, P);
                    if isfinite(C) && C > bestC
                        bestC = C; bestLag = lag; bestSwap = (so == 2);
                    end
                end
            end
            swapped = bestSwap;
            start = nDropNom + 1 + bestLag;
            if start < 1
                eqOut = eqOutFull([], :); nSubUsable = 0; return;
            end
            nSubUsable = floor((N - start + 1) / P.SUBFRAME_SYMS);
            if nSubUsable < 1
                eqOut = eqOutFull([], :); return;
            end
            eqOut = eqOutFull(start : start + nSubUsable * P.SUBFRAME_SYMS - 1, :);
            if bestSwap
                eqOut = eqOut(:, [2 1]);   % undo the CMA polarisation swap
            end
        end

        function C = pilotCoherence(x, start, nSub, pilotsRef, P)
            %PILOTCOHERENCE  Differential pilot-phase coherence at a candidate
            %   subframe-boundary start (1-based sample index).  Returns a
            %   value near 1 when the block-start symbols land on the known
            %   pilots and near 0 on random data.
            NBsub = size(pilotsRef, 1);          % blocks per subframe (=116)
            NB    = nSub * NBsub;
            z = complex(zeros(NB, 1));
            for b = 1:NB
                bb = mod(b - 1, NBsub) + 1;          % pilot index in subframe
                p  = start + (b - 1) * P.BlockLen_CR; % block-start sample
                acc = complex(0, 0);
                for pol = 1:P.N_pol
                    acc = acc + conj(pilotsRef(bb, pol)) * x(p, pol);
                end
                z(b) = acc;
            end
            d = z(2:end) .* conj(z(1:end-1));
            denom = sum(abs(d));
            if denom > 0
                C = abs(sum(d)) / denom;
            else
                C = NaN;
            end
        end

        function BER = computeBER(decoded, refSyms)
            % 32 pi/16 phase rotations + polarisation swap, picking the
            % alignment that minimises BER per polarisation.  Same scheme
            % as combined_eq_clk_sweep.computeBER.
            nPol   = size(refSyms, 2);
            totErr = 0; totBits = 0;
            for p = 1:nPol
                refBits = modem.symbolsToBits(refSyms(:, p));
                best    = Inf;
                for q = 1:size(decoded, 2)
                    for kk = 0:31
                        rotated = decoded(:, q) .* exp(-1j * kk * pi/16);
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

        % ------------------------------------------------------------------
        %  MEX naming and per-precision builds
        % ------------------------------------------------------------------
        function name = mexEqName(staticFL, clkFL, adaptFL)
            % Unique MEX per (StaticFL, ClkFL, AdaptFL).  po2 is a runtime
            % argument and is NOT part of the name.  This pipeline uses the
            % ARCHIVED circular-wrap Godard block
            % (eq_clk.combined_cd_fd_godard_adaptive_fxp_old); the "_old" tag
            % keeps its binaries distinct from the streaming version's cache.
            name = sprintf( ...
                'combined_cd_fd_godard_adaptive_fxp_old_a%02dc%02ds%02d_mex', ...
                adaptFL, clkFL, staticFL);
        end

        function name = mexFRName(fl)
            name = sprintf('differential_kay_fxp_fc%02d_mex', fl);
        end

        function name = mexCRName(fl)
            name = sprintf('pilots_only_fxp_fc%02d_mex', fl);
        end

        function buildEqMex(P, cfgCoder, staticFL, clkFL, adaptFL, mexBase, repoRoot)
            % Reuses build_eq_clk_combined_cd_fd_godard_adaptive_fxp_old_mex
            % (the ARCHIVED circular-wrap Godard block) then renames the output
            % to the precision-suffixed name.  Per-section precisions are
            % passed as a composite struct config.
            B = struct();
            B.FxpConfig_CombGodard = struct( ...
                'Static',  struct('WL', P.IntBits + staticFL, 'FL', staticFL), ...
                'Godard',  struct('WL', P.IntBits + clkFL,    'FL', clkFL), ...
                'AdaptEq', struct('WL', P.IntBits + adaptFL,  'FL', adaptFL));
            B.SpS       = P.SpS;
            B.NFFT      = P.NFFT;
            B.NOverlap  = 2 * ceil((P.NCD - 1) / 2);
            B.D         = P.D;
            B.L         = P.L_km;
            B.CWL       = P.CWL;
            B.Rs        = P.Rs;
            B.Rolloff   = P.Rolloff;
            B.N_pol     = P.N_pol;
            B.po2Twiddle = false;
            B.cfoEnable  = logical(P.CfoEnable);
            B.Ns        = P.N_sub_target * P.SUBFRAME_SYMS;
            % ki/kp here only type the codegen prototype as a double scalar;
            % the runtime value is picked per design row's po2 flag.
            B.CR_ki     = P.ki_godard_po2_off;
            B.CR_kp     = P.kp_godard_po2_off;
            B.AEQ_NTaps          = P.NTapsAEQ;
            B.AEQ_Mu             = P.MuAEQ;
            B.AEQ_SingleSpike    = P.SingleSpike;
            B.AEQ_N1             = P.N1AEQ;
            B.AEQ_NOut           = P.NOutAEQ;
            B.AEQ_SignOnly       = P.SignOnly;
            B.AEQ_UpdateStep     = 1;
            B.AEQ_PLanes         = P.PLanesAEQ;
            B.AEQ_Mode           = 0;
            B.AEQ_BlockLen       = P.PLanesAEQ;
            B.AEQ_SubframeBlocks = 0;
            B.AEQ_RScale         = 1;

            build_eq_clk_combined_cd_fd_godard_adaptive_fxp_old_mex(B, cfgCoder);

            srcDir = fullfile(repoRoot, 'src', '+eq_clk');
            ext    = ['.' mexext];
            srcFile = fullfile(srcDir, ['combined_cd_fd_godard_adaptive_fxp_old_mex' ext]);
            dstFile = fullfile(srcDir, [mexBase ext]);
            if isfile(srcFile)
                if isfile(dstFile), delete(dstFile); end
                movefile(srcFile, dstFile);
            end
        end

        function buildFRMex(P, cfgCoder, fl, mexBase, repoRoot)
            B = struct();
            B.Rs            = P.Rs;
            B.N_pol         = P.N_pol;
            B.TrainingLen   = P.TrainingLen;
            B.FxpConfig_FR  = struct('WL', P.IntBits + fl, 'FL', fl);
            B.CordicIts     = fl;   % CORDIC iterations = FR fractional length
            B.MaxFreq       = P.MaxFreq;

            build_freq_recovery_differential_kay_fxp_mex(B, cfgCoder);

            srcDir  = fullfile(repoRoot, 'src', '+freq_recovery');
            ext     = ['.' mexext];
            srcFile = fullfile(srcDir, ['differential_kay_fxp_mex' ext]);
            dstFile = fullfile(srcDir, [mexBase ext]);
            if isfile(srcFile)
                if isfile(dstFile), delete(dstFile); end
                movefile(srcFile, dstFile);
            end
        end

        function buildCRMex(P, cfgCoder, fl, mexBase, repoRoot)
            B = struct();
            B.N_pol        = P.N_pol;
            B.BlockLen     = P.BlockLen_CR;
            B.FxpConfig_PO = struct('WL', P.IntBits + fl, 'FL', fl);
            B.CordicIts    = fl;   % CORDIC iterations = CR fractional length

            build_carrier_recovery_pilots_only_fxp_mex(B, cfgCoder);

            srcDir  = fullfile(repoRoot, 'src', '+carrier_recovery');
            ext     = ['.' mexext];
            srcFile = fullfile(srcDir, ['pilots_only_fxp_mex' ext]);
            dstFile = fullfile(srcDir, [mexBase ext]);
            if isfile(srcFile)
                if isfile(dstFile), delete(dstFile); end
                movefile(srcFile, dstFile);
            end
        end

    end
end
