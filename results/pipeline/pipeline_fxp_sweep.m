classdef pipeline_fxp_sweep < matlab.unittest.TestCase
%PIPELINE_FXP_SWEEP  Full receiver pipeline BER under fxp precision sweeps.
%
%   Pipeline:
%       CPON Tx -> channel -> Gardner-combined eq -> drop 1 subframe ->
%       pilot rescale -> FR -> pilots-only CR -> decision -> BER.
%
%   Equaliser block under test:
%       eq_clk.combined_cd_fd_gardner_adaptive_fxp  (po2 false and true).
%   FR blocks under test:
%       freq_recovery.differential_kay_fxp  (training-aided)
%       freq_recovery.fft_search_fxp        (blind)
%   CR block:
%       carrier_recovery.pilots_only_fxp.
%
%   For each (po2, FR algo) combination, the user supplies a 2D array of
%   (FL_eq, FL_frcr) pairs (one row per precision combination).  FL_eq
%   applies jointly to the static CD/MF, the in-block clock recovery, and
%   the adaptive equaliser type tables.  FL_frcr applies jointly to the FR
%   and CR type tables.  Integer bits are held at IntBits, so the per-stage
%   word length is IntBits + FL.
%
%   CPON pipeline integration
%     The equaliser sees pilots / training at the data energy of ±1±1j (the
%     output of modem.modulate).  After the equaliser, one full CPON
%     subframe (3712 symbols) is discarded so the remaining stream starts
%     at a subframe boundary and contains an integer number of subframes;
%     the pilot positions in the discarded-aligned stream are then scaled
%     ×3 (to ±3±3j, per the CPON spec) so that the pilots_only CR sees
%     the design-intended pilot amplitude.
%
%   Codegen MEX dispatch
%     Each (function, FL) combination is built into its own MEX up front
%     in TestClassSetup.  Existing MEX files at the per-precision suffix
%     are skipped (delete to force a rebuild).  The sweep itself dispatches
%     via str2func at runtime, so per-config overhead is just one MATLAB
%     -> MEX call.
%
%   Output
%     pipeline_fxp_sweep.mat with a table whose rows correspond to
%     (algo_idx, prec_row) and per-trial BERs in a cell column.

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
        PMD_seed = 12345
        L_km     = 80
        SFO_ppm  = 40
        CFO_GHz  = 3          % worst-case CPON CFO
        LW_Hz    = 1e6           % 1 MHz combined laser linewidth

        % --- Pulse shaping --------------------------------------------
        Rolloff = 0.25
        Span    = 10

        % --- Monte-Carlo ----------------------------------------------
        N_sub_target = 8         % CPON subframes per trial (Tx)
        NTrials      = 5
        SNR_dB_vec   = 0 : 1 : 20
        % --- Equaliser (overlap-save + Gardner + adaptive) -----------
        NFFT        = 128
        NCD         = 22
        NLanesGard  = 32

        % Gardner PI gains, one pair per po2-twiddle value.  Selected at
        % run-time by cfg.po2; tune these in the design sweep first.
        ki_gardner_po2_off = 1e-7
        kp_gardner_po2_off = 1e-6
        ki_gardner_po2_on  = 1e-6
        kp_gardner_po2_on  = 1e-6

        % Adaptive eq (matches the existing combined_eq_clk_sweep)
        NTapsAEQ    = 1
        MuAEQ       = 1e-3
        N1AEQ       = 500
        NOutAEQ     = 1000
        SignOnly    = true
        SingleSpike = true
        PLanesAEQ   = 32

        % --- Frequency recovery ---------------------------------------
        FR_Nfft       = 512
        FR_Po2Twiddle = false
        FR_BlindD     = 512
        MaxFreq       = 0.1
        TrainingLen   = 11

        % --- Carrier recovery (pilots_only) ---------------------------
        CordicIts   = 16
        BlockLen_CR = 32

        % --- Fixed point ----------------------------------------------
        IntBits = 16

        % --- MEX build control ----------------------------------------
        %  false: reuse any existing per-precision MEX (*_eq<FL>_mex /
        %         *_fc<FL>_mex) — fast reruns.
        %  true : force a fresh codegen of EVERY fxp MEX at startup,
        %         ignoring the cache.  Set this after editing any *_fxp.m
        %         source (e.g. the per-pol carrier_recovery.pilots_only_fxp
        %         change) so the stale cached binaries are regenerated.
        ForceRebuild = true

        % --- Algorithm combinations --------------------------------------
        %  Parallel arrays.  Entry i = (AlgoPo2(i), AlgoFR{i}) is exercised
        %  on the full Cartesian grid FL_eq_vec x FL_frcr_vec.
        AlgoPo2 = [false, false, true, true]
        AlgoFR  = {'differential_kay', 'fft_search_blind', ...
                   'differential_kay', 'fft_search_blind'}

        % --- Precision grid (1D per block, full Cartesian product) -----
        %  FL_eq_vec   sweeps the equaliser fxp fractional length
        %              (jointly applied to T.Static, T.Clk, T.AdaptEq).
        %  FL_frcr_vec sweeps the FR/CR fxp fractional length
        %              (jointly applied to FR T and CR T).
        FL_eq_vec   = [4,6,8]
        FL_frcr_vec = [2,4,6]

        % --- CPON subframe constants (must match modem.modulate) -------
        SUBFRAME_SYMS = 3712
        BLOCK_LEN     = 32
        N_BLOCKS      = 116
    end

    %% ================================================================
    %  Setup: build one MEX per (function, FL) combination
    %% ================================================================
    methods (TestClassSetup)
        function setupAndBuildMex(testCase)
            here     = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(here));
            addpath(genpath(fullfile(repoRoot, 'src')));
            addpath(fullfile(repoRoot, 'build'));

            P = pipeline_fxp_sweep.extractParams(testCase);

            % --- Gather unique FL values per builder ------------------
            allEqFL = unique(P.FL_eq_vec(:));
            allCrFL = unique(P.FL_frcr_vec(:));
            FL_dk   = [];
            FL_fft  = [];
            for ai = 1:numel(P.AlgoFR)
                if strcmp(P.AlgoFR{ai}, 'differential_kay')
                    FL_dk = [FL_dk; P.FL_frcr_vec(:)];   %#ok<AGROW>
                elseif strcmp(P.AlgoFR{ai}, 'fft_search_blind')
                    FL_fft = [FL_fft; P.FL_frcr_vec(:)]; %#ok<AGROW>
                end
            end
            FL_dk  = unique(FL_dk);
            FL_fft = unique(FL_fft);

            cfgCoder = coder.config('mex');
            cfgCoder.GenerateReport = false;

            nTotal = numel(allEqFL) + numel(FL_dk) + numel(FL_fft) + numel(allCrFL);
            bi = 0; tStart = tic;
            if P.ForceRebuild
                fprintf('\n=== Building %d MEX combos (ForceRebuild ON: cache ignored) ===\n', nTotal);
            else
                fprintf('\n=== Building %d MEX combos (using cache where present) ===\n', nTotal);
            end

            % --- Equaliser builds --------------------------------------
            for fl = allEqFL(:).'
                bi = bi + 1;
                mexBase = pipeline_fxp_sweep.mexEqName(fl);
                mexPath = fullfile(repoRoot, 'src', '+eq_clk', ...
                    [mexBase '.mexw64']);
                if isfile(mexPath) && ~P.ForceRebuild
                    fprintf('  [%2d/%2d] %s -> cached\n', bi, nTotal, mexBase);
                    continue;
                end
                fprintf('  [%2d/%2d] %s ', bi, nTotal, mexBase);
                t1 = tic;
                pipeline_fxp_sweep.buildEqMex( ...
                    P, cfgCoder, fl, mexBase, repoRoot);
                fprintf('(%.0fs)\n', toc(t1));
            end

            % --- differential_kay FR builds ---------------------------
            for fl = FL_dk(:).'
                bi = bi + 1;
                mexBase = pipeline_fxp_sweep.mexFRName('differential_kay', fl);
                mexPath = fullfile(repoRoot, 'src', '+freq_recovery', ...
                    [mexBase '.mexw64']);
                if isfile(mexPath) && ~P.ForceRebuild
                    fprintf('  [%2d/%2d] %s -> cached\n', bi, nTotal, mexBase);
                    continue;
                end
                fprintf('  [%2d/%2d] %s ', bi, nTotal, mexBase);
                t1 = tic;
                pipeline_fxp_sweep.buildFRMex(P, cfgCoder, ...
                    'differential_kay', fl, mexBase, repoRoot);
                fprintf('(%.0fs)\n', toc(t1));
            end

            % --- fft_search_blind FR builds ---------------------------
            for fl = FL_fft(:).'
                bi = bi + 1;
                mexBase = pipeline_fxp_sweep.mexFRName('fft_search_blind', fl);
                mexPath = fullfile(repoRoot, 'src', '+freq_recovery', ...
                    [mexBase '.mexw64']);
                if isfile(mexPath) && ~P.ForceRebuild
                    fprintf('  [%2d/%2d] %s -> cached\n', bi, nTotal, mexBase);
                    continue;
                end
                fprintf('  [%2d/%2d] %s ', bi, nTotal, mexBase);
                t1 = tic;
                pipeline_fxp_sweep.buildFRMex(P, cfgCoder, ...
                    'fft_search_blind', fl, mexBase, repoRoot);
                fprintf('(%.0fs)\n', toc(t1));
            end

            % --- pilots_only CR builds --------------------------------
            for fl = allCrFL(:).'
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
            cfgs = pipeline_fxp_sweep.buildCfgs(P);
            NCFG = numel(cfgs);
            NSNR = numel(P.SNR_dB_vec);

            fprintf(['\n=== Pipeline fxp sweep: %d configs x %d SNRs ', ...
                     'x %d trials ===\n'], NCFG, NSNR, P.NTrials);

            ber = nan(NCFG, NSNR, P.NTrials);
            for ci = 1:NCFG
                cfg = cfgs(ci);
                fprintf('[%2d/%2d] po2=%d fr=%-18s FL_eq=%2d FL_frcr=%2d\n', ...
                    ci, NCFG, cfg.po2, cfg.fr_algo, cfg.fl_eq, cfg.fl_frcr);
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
            algo_idx   = nan(NCFG, 1);
            po2        = false(NCFG, 1);
            fr_algo    = strings(NCFG, 1);
            fl_eq      = nan(NCFG, 1);
            eq_wl      = nan(NCFG, 1);
            fl_frcr    = nan(NCFG, 1);
            frcr_wl    = nan(NCFG, 1);
            ber_cell   = cell(NCFG, 1);
            ber_mean   = cell(NCFG, 1);
            for ci = 1:NCFG
                cfg = cfgs(ci);
                algo_idx(ci) = cfg.algo_idx;
                po2(ci)      = cfg.po2;
                fr_algo(ci)  = string(cfg.fr_algo);
                fl_eq(ci)    = cfg.fl_eq;
                eq_wl(ci)    = P.IntBits + cfg.fl_eq;
                fl_frcr(ci)  = cfg.fl_frcr;
                frcr_wl(ci)  = P.IntBits + cfg.fl_frcr;
                % ber(ci, :, :) is [NSNR x NTrials]
                ber_cell{ci} = squeeze(ber(ci, :, :));
                ber_mean{ci} = mean(squeeze(ber(ci, :, :)), 2, 'omitnan').';
            end
            tbl = table(algo_idx, po2, fr_algo, fl_eq, eq_wl, ...
                fl_frcr, frcr_wl, ber_cell, ber_mean);
            tbl.Properties.VariableNames{'ber_cell'} = 'ber';
            tbl.Properties.VariableNames{'ber_mean'} = 'ber_mean_per_snr';

            S.tbl        = tbl;
            S.SNR_dB_vec = P.SNR_dB_vec;
            S.CFO_GHz    = P.CFO_GHz;
            S.NTrials    = P.NTrials;
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
            P.Rolloff       = tc.Rolloff;
            P.Span          = tc.Span;
            P.N_sub_target  = tc.N_sub_target;
            P.NTrials       = tc.NTrials;
            P.SNR_dB_vec    = tc.SNR_dB_vec;
            P.NFFT          = tc.NFFT;
            P.NCD           = tc.NCD;
            P.NLanesGard    = tc.NLanesGard;
            P.ki_gardner_po2_off = tc.ki_gardner_po2_off;
            P.kp_gardner_po2_off = tc.kp_gardner_po2_off;
            P.ki_gardner_po2_on  = tc.ki_gardner_po2_on;
            P.kp_gardner_po2_on  = tc.kp_gardner_po2_on;
            P.NTapsAEQ      = tc.NTapsAEQ;
            P.MuAEQ         = tc.MuAEQ;
            P.N1AEQ         = tc.N1AEQ;
            P.NOutAEQ       = tc.NOutAEQ;
            P.SignOnly      = tc.SignOnly;
            P.SingleSpike   = tc.SingleSpike;
            P.PLanesAEQ     = tc.PLanesAEQ;
            P.FR_Nfft       = tc.FR_Nfft;
            P.FR_Po2Twiddle = tc.FR_Po2Twiddle;
            P.FR_BlindD     = tc.FR_BlindD;
            P.MaxFreq       = tc.MaxFreq;
            P.TrainingLen   = tc.TrainingLen;
            P.CordicIts     = tc.CordicIts;
            P.BlockLen_CR   = tc.BlockLen_CR;
            P.IntBits       = tc.IntBits;
            P.ForceRebuild  = tc.ForceRebuild;
            P.AlgoPo2       = tc.AlgoPo2;
            P.AlgoFR        = tc.AlgoFR;
            P.FL_eq_vec     = tc.FL_eq_vec;
            P.FL_frcr_vec   = tc.FL_frcr_vec;
            P.SUBFRAME_SYMS = tc.SUBFRAME_SYMS;
            P.BLOCK_LEN     = tc.BLOCK_LEN;
            P.N_BLOCKS      = tc.N_BLOCKS;
        end

        function cfgs = buildCfgs(P)
            % Full Cartesian product of (algo combo) x FL_eq_vec x FL_frcr_vec.
            cfgs = struct('algo_idx', {}, 'po2', {}, 'fr_algo', {}, ...
                'fl_eq', {}, 'fl_frcr', {});
            for ai = 1:numel(P.AlgoFR)
                for fe = P.FL_eq_vec(:).'
                    for ff = P.FL_frcr_vec(:).'
                        cfgs(end + 1).algo_idx = ai;             %#ok<AGROW>
                        cfgs(end).po2         = logical(P.AlgoPo2(ai));
                        cfgs(end).fr_algo     = P.AlgoFR{ai};
                        cfgs(end).fl_eq       = fe;
                        cfgs(end).fl_frcr     = ff;
                    end
                end
            end
        end

        % ------------------------------------------------------------------
        %  Per-trial run: channel + eq + drop subframe + scale + FR + CR + BER
        % ------------------------------------------------------------------
        function ber = runOneTrial(P, cfg, SNR_dB, trialSeed)
            try
                [rxSig, txSymbols, training, pilotsRef, ~] = ...
                    pipeline_fxp_sweep.genChannel(P, SNR_dB, trialSeed);

                % --- 1. Equaliser MEX ---------------------------------
                FxpEq = struct('WL', P.IntBits + cfg.fl_eq, 'FL', cfg.fl_eq);
                Tcfg = struct('Static',  FxpEq, ...
                              'Clk',     FxpEq, ...
                              'AdaptEq', FxpEq);
                T_eq = eq_clk.combined_cd_fd_gardner_adaptive_fxp_types(Tcfg);

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
                    'SubframeBlocks', double(0));

                nOv = 2 * ceil((P.NCD - 1) / 2);
                eqMex = str2func(['eq_clk.' ...
                    pipeline_fxp_sweep.mexEqName(cfg.fl_eq)]);

                % Pick the Gardner PI gains tuned for this po2 setting.
                if cfg.po2
                    kiSel = P.ki_gardner_po2_on;
                    kpSel = P.kp_gardner_po2_on;
                else
                    kiSel = P.ki_gardner_po2_off;
                    kpSel = P.kp_gardner_po2_off;
                end

                NsymTx = size(txSymbols, 1);
                [yFi, ~] = eqMex(rxSig_fi, double(P.SpS), ...
                    double(P.NFFT), double(nOv), ...
                    double(P.D), double(P.L_km), double(P.CWL), ...
                    double(P.Rs), double(P.Rolloff), ...
                    double(kiSel), double(kpSel), ...
                    double(NsymTx), double(P.NLanesGard), ...
                    adaptOpts, true, logical(cfg.po2), T_eq);
                eqOut = double(yFi);

                % --- 2. Drop one subframe; align + integer subframes --
                nDropEq = P.SUBFRAME_SYMS - P.NOutAEQ;
                if nDropEq < 0
                    nDropEq = 0;
                end
                if size(eqOut, 1) <= nDropEq
                    ber = NaN;
                    return;
                end
                eqOut = eqOut(nDropEq + 1 : end, :);
                nSubUsable = floor(size(eqOut, 1) / P.SUBFRAME_SYMS);
                if nSubUsable < 1
                    ber = NaN;
                    return;
                end
                nUsable = nSubUsable * P.SUBFRAME_SYMS;
                eqOut   = eqOut(1 : nUsable, :);

                % --- 3. Scale pilot + training positions ×3 -----------
                %  Pilots: every block-start symbol of each subframe.
                %  Training: the first TrainingLen symbols of every
                %  subframe (TS1 coincides with the block-1 pilot, the
                %  remaining 10 are pure training).
                rescaleIdx = pipeline_fxp_sweep.pilotTrainingIndices( ...
                    nSubUsable, P);
                eqOut(rescaleIdx, :) = 3 * eqOut(rescaleIdx, :);

                % --- 4. Frequency recovery ---------------------------
                FxpFRCR = struct('WL', P.IntBits + cfg.fl_frcr, ...
                                 'FL', cfg.fl_frcr);
                T_fr = freq_recovery.fxp_types(FxpFRCR);
                eqOut_fi    = cast(eqOut,        'like', T_fr.x);
                training_fi = cast(3 * training, 'like', T_fr.x);

                frMex = str2func(['freq_recovery.' ...
                    pipeline_fxp_sweep.mexFRName(cfg.fr_algo, cfg.fl_frcr)]);

                switch cfg.fr_algo
                    case 'differential_kay'
                        [frOut_fi, ~] = frMex(eqOut_fi, training_fi, ...
                            double(P.Rs), double(P.CordicIts), T_fr, ...
                            true, double(0), double(P.MaxFreq));
                    case 'fft_search_blind'
                        [frOut_fi, ~] = frMex(eqOut_fi, training_fi, ...
                            double(P.Rs), double(P.FR_Nfft), ...
                            logical(P.FR_Po2Twiddle), double(P.CordicIts), ...
                            double(P.MaxFreq), T_fr, false, ...
                            double(P.FR_BlindD));
                end

                % --- 5. Carrier recovery (pilots_only) ---------------
                T_cr = carrier_recovery.fxp_types(FxpFRCR);
                frOut_for_cr = cast(double(frOut_fi), 'like', T_cr.x);

                pilotsAll = 3 * repmat(pilotsRef, nSubUsable, 1);
                pilots_fi = cast(pilotsAll, 'like', T_cr.x);

                crMex = str2func(['carrier_recovery.' ...
                    pipeline_fxp_sweep.mexCRName(cfg.fl_frcr)]);

                [crOut_fi, ~] = crMex(frOut_for_cr, double(P.N_pol), ...
                    double(P.BlockLen_CR), pilots_fi, ...
                    double(P.CordicIts), T_cr);
                crOut = double(crOut_fi);

                % --- 6. BER vs reference (subframes 2..end of Tx) ----
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
        %  Channel: CPON modulator + RRC + CD + PMD + CFO + SFO + AWGN +
        %  phase noise (the equaliser sees pilots/training at ±1±1j).
        % ------------------------------------------------------------------
        function [rxSig, symbols, training, pilotsRef, nSubTx] = ...
                genChannel(P, SNR_dB, trialSeed)
            rng(1000 * trialSeed + round(SNR_dB) + 13);

            CPON_BITS_PER_SF = 3586 * P.N_pol * 2;
            nBits = P.N_sub_target * CPON_BITS_PER_SF;
            bits  = randi([0 1], nBits, 1);

            [symbols, pilotsRef, training, nSubTx] = modem.modulate(bits);

            txSig = modem.rrcPulse(symbols, P.SpS, P.Rolloff, P.Span);
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
            rxSig = channel.add_phase_noise(rxSig, P.Rs * P.SpS, P.LW_Hz);
        end

        function idx = pilotTrainingIndices(nSubUsable, P)
            %  Indices of pilot + training symbols in the aligned stream.
            %  Per subframe: positions 1..TrainingLen are training (TS1 is
            %  also pilot index 1), then positions (b-1)*BLOCK_LEN + 1 for
            %  b = 2..N_BLOCKS are the remaining pilots.
            nPerSub = P.TrainingLen + (P.N_BLOCKS - 1);
            idx = zeros(nSubUsable * nPerSub, 1);
            k = 0;
            for sf = 0 : nSubUsable - 1
                base = sf * P.SUBFRAME_SYMS;
                % Training (positions 1..TrainingLen of the subframe)
                for t = 1 : P.TrainingLen
                    k = k + 1;
                    idx(k) = base + t;
                end
                % Pilots for blocks 2..N_BLOCKS (block 1's pilot is TS1).
                for b = 2 : P.N_BLOCKS
                    k = k + 1;
                    idx(k) = base + (b - 1) * P.BLOCK_LEN + 1;
                end
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
        function name = mexEqName(fl)
            name = sprintf('combined_cd_fd_gardner_adaptive_fxp_eq%02d_mex', fl);
        end

        function name = mexFRName(algo, fl)
            switch algo
                case 'differential_kay'
                    name = sprintf('differential_kay_fxp_fc%02d_mex', fl);
                case 'fft_search_blind'
                    name = sprintf('fft_search_fxp_fc%02d_mex', fl);
                otherwise
                    error('pipeline_fxp_sweep:badFR', ...
                          'Unknown FR algo: %s', algo);
            end
        end

        function name = mexCRName(fl)
            name = sprintf('pilots_only_fxp_fc%02d_mex', fl);
        end

        function buildEqMex(P, cfgCoder, fl, mexBase, repoRoot)
            % Reuses build_eq_clk_combined_cd_fd_gardner_adaptive_fxp_mex
            % then renames the output to the precision-suffixed name.
            B = struct();
            B.FxpConfig_CombGardner = struct( ...
                'WL', P.IntBits + fl, 'FL', fl);
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
            B.cfoEnable  = false;
            B.Ns        = P.N_sub_target * P.SUBFRAME_SYMS;
            % ki/kp here are only used to type the codegen prototype as a
            % double scalar; the runtime value is picked per cfg.po2.
            B.CR_ki     = P.ki_gardner_po2_off;
            B.CR_kp     = P.kp_gardner_po2_off;
            B.CR_NLanes = P.NLanesGard;
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

            build_eq_clk_combined_cd_fd_gardner_adaptive_fxp_mex(B, cfgCoder);

            srcDir = fullfile(repoRoot, 'src', '+eq_clk');
            ext    = ['.' mexext];
            srcFile = fullfile(srcDir, ['combined_cd_fd_gardner_adaptive_fxp_mex' ext]);
            dstFile = fullfile(srcDir, [mexBase ext]);
            if isfile(srcFile)
                if isfile(dstFile), delete(dstFile); end
                movefile(srcFile, dstFile);
            end
        end

        function buildFRMex(P, cfgCoder, algo, fl, mexBase, repoRoot)
            B = struct();
            B.Rs            = P.Rs;
            B.N_pol         = P.N_pol;
            B.TrainingLen   = P.TrainingLen;
            B.FxpConfig_FR  = struct('WL', P.IntBits + fl, 'FL', fl);
            B.CordicIts     = P.CordicIts;
            B.MaxFreq       = P.MaxFreq;
            B.FR_Nfft       = P.FR_Nfft;
            B.FR_Po2Twiddle = P.FR_Po2Twiddle;
            B.FR_BlindD     = P.FR_BlindD;

            srcDir = fullfile(repoRoot, 'src', '+freq_recovery');
            ext    = ['.' mexext];
            switch algo
                case 'differential_kay'
                    build_freq_recovery_differential_kay_fxp_mex(B, cfgCoder);
                    srcFile = fullfile(srcDir, ['differential_kay_fxp_mex' ext]);
                case 'fft_search_blind'
                    build_freq_recovery_fft_search_fxp_mex(B, cfgCoder);
                    srcFile = fullfile(srcDir, ['fft_search_fxp_mex' ext]);
            end
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
            B.CordicIts    = P.CordicIts;

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
