classdef combined_eq_clk_fxp_sweep < matlab.unittest.TestCase
%COMBINED_EQ_CLK_FXP_SWEEP  Fixed-point precision sweep for the two
%   combined CD-FD + clock-recovery + adaptive-CMA blocks.
%
%   Each algorithm:
%       cd_gardner_cma : eq_clk.combined_cd_fd_gardner_adaptive_fxp
%       cd_godard_cma  : eq_clk.combined_cd_fd_godard_adaptive_fxp
%   is exercised with the po2-twiddle option both on and off, over a
%   grid of fixed-point precisions:
%
%       EqPrec_vec  - applied jointly to T.Static and T.AdaptEq
%       ClkPrec_vec - applied to T.Clk (Gardner) or T.Godard
%
%   CFO is held at 3 GHz for every configuration; the combined block
%   performs its one-shot coarse CFO correction in floating point
%   internally and the post-block exact-CFO removal absorbs the residual.
%
%   Codegen MEX dispatch
%       Each precision combo bakes the per-section fi numerictypes into
%       a separate MEX (Tcfg is a -args constant at codegen time).  The
%       TestClassSetup loops over (block, EqFL, ClkFL) and produces one
%       MEX per combo at
%           src/+eq_clk/<base>_e<EqFL>c<ClkFL>_mex.mexw64
%       where <base> is the fxp function's name.  Existing MEX files
%       are NOT rebuilt; delete them by hand to force a refresh.  The
%       sweep itself dispatches via feval at runtime, so per-config
%       overhead is just one MATLAB->MEX call.
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
        Ns          = 37500         % symbols per polarisation per trial
        NTrials     = 5
        SNR_dB_vec  = 0 : 2 : 20

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
        ki             = [1e-6 1e-7; ...
                          1e-4 1e-5]
        kp             = [1e-4 1e-4; ...
                          1e-4 1e-5]

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

        % --- Precision sweep ----------------------------------------
        %  Integer bits are fixed at NIntBits; each entry in *_FL_vec is
        %  the fractional length to sweep.  At each point the test forms
        %  the per-section struct as struct('WL', NIntBits + FL, 'FL', FL)
        %  and passes it through the composite fxp types builder.
        %
        %  Note on the adaptive equaliser: the struct path of
        %  adaptive_eq.equalize_fxp_types interprets WL/FL as the
        %  *gradient* precision (T.grad) and pins the data path at high
        %  precision.  EqFL_vec therefore sweeps the joint
        %  static-equaliser precision and adaptive-equaliser gradient
        %  precision.
        %  Default grids kept to 3 x 3 = 9 combos per block (18 builds
        %  total) so the up-front codegen phase stays roughly within
        %  ~10 minutes on a typical workstation.  Extend as needed.
        NIntBits   = 16
        EqFL_vec   = [2, 4, 6, 8, 10]
        ClkFL_vec  = [2, 4, 6, 8, 10]

        % --- FEC threshold used to score designs --------------------
        FEC_BER = 2e-2
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

            NBlk = numel(P.blocks);
            NEQ  = numel(P.EqFL_vec);
            NCLK = numel(P.ClkFL_vec);
            nTotal = NBlk * NEQ * NCLK;
            bi = 0;
            tStart = tic;
            fprintf('\n=== Building %d MEX combos ===\n', nTotal);
            for blk = 1:NBlk
                blkName = P.blocks{blk};
                for ei = 1:NEQ
                    for cli = 1:NCLK
                        bi = bi + 1;
                        eqFL  = P.EqFL_vec(ei);
                        clkFL = P.ClkFL_vec(cli);
                        mexBase = combined_eq_clk_fxp_sweep.mexBaseName( ...
                            blkName, eqFL, clkFL);
                        mexPath = fullfile(repoRoot, 'src', '+eq_clk', ...
                            [mexBase '.mexw64']);
                        if isfile(mexPath)
                            fprintf('  [%2d/%2d] %s -> cached\n', ...
                                bi, nTotal, mexBase);
                            continue;
                        end
                        fprintf('  [%2d/%2d] %s ', bi, nTotal, mexBase);
                        t1 = tic;
                        combined_eq_clk_fxp_sweep.buildOneMex( ...
                            P, cfgCoder, blkName, eqFL, clkFL, ...
                            mexBase, repoRoot);
                        fprintf('(%.0fs)\n', toc(t1));
                    end
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

            NBlk  = numel(P.blocks);
            NPO2  = numel(P.Po2Twiddle_vec);
            NEQ   = numel(P.EqFL_vec);
            NCLK  = numel(P.ClkFL_vec);
            NSNR  = numel(P.SNR_dB_vec);
            NCFG  = NBlk * NPO2 * NEQ * NCLK;

            % Validate the manual-design arrays
            assert(numel(P.NCD)   == NBlk && ...
                   numel(P.NTaps) == NBlk && ...
                   isequal(size(P.ki), [NBlk, NPO2]) && ...
                   isequal(size(P.kp), [NBlk, NPO2]), ...
                'blocks / NCD / NTaps / ki / kp size mismatch.');

            % --- Build cfg array ---------------------------------------
            cfgs = struct('block',{}, 'NCD',{}, 'NOverlap',{}, ...
                'NTaps',{}, 'Po2Twiddle',{}, 'ki',{}, 'kp',{}, ...
                'EqFL',{}, 'ClkFL',{}, ...
                'block_idx',{}, 'po2_idx',{}, ...
                'eq_idx',{}, 'clk_idx',{});
            ci = 0;
            for bi = 1:NBlk
                nCd = P.NCD(bi);
                nOv = 2 * ceil((nCd - 1) / 2);
                nT  = P.NTaps(bi);
                for pi = 1:NPO2
                    ki_v = P.ki(bi, pi);
                    kp_v = P.kp(bi, pi);
                    for ei = 1:NEQ
                        for cli = 1:NCLK
                            ci = ci + 1;
                            cfgs(ci).block      = P.blocks{bi};
                            cfgs(ci).NCD        = nCd;
                            cfgs(ci).NOverlap   = nOv;
                            cfgs(ci).NTaps      = nT;
                            cfgs(ci).Po2Twiddle = logical(P.Po2Twiddle_vec(pi));
                            cfgs(ci).ki         = ki_v;
                            cfgs(ci).kp         = kp_v;
                            cfgs(ci).EqFL       = P.EqFL_vec(ei);
                            cfgs(ci).ClkFL      = P.ClkFL_vec(cli);
                            cfgs(ci).block_idx  = bi;
                            cfgs(ci).po2_idx    = pi;
                            cfgs(ci).eq_idx     = ei;
                            cfgs(ci).clk_idx    = cli;
                        end
                    end
                end
            end

            fprintf('\n=== FXP precision sweep, CFO = %.2f GHz ===\n', P.CFO_GHz);
            fprintf('Configurations: %d  (blocks=%d, po2=%d, EqFL=%d, ClkFL=%d)  IntBits=%d\n', ...
                NCFG, NBlk, NPO2, NEQ, NCLK, P.NIntBits);
            for ci = 1:NCFG
                fprintf(['  %-16s NCD=%2d NTaps=%d po2=%d  ki=%.2e kp=%.2e  ' ...
                    'EQ WL=%d FL=%d  Clk WL=%d FL=%d\n'], ...
                    cfgs(ci).block, cfgs(ci).NCD, cfgs(ci).NTaps, ...
                    cfgs(ci).Po2Twiddle, cfgs(ci).ki, cfgs(ci).kp, ...
                    P.NIntBits + cfgs(ci).EqFL, cfgs(ci).EqFL, ...
                    P.NIntBits + cfgs(ci).ClkFL, cfgs(ci).ClkFL);
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
            eq_wl      = nan(NCFG, 1);
            eq_fl      = nan(NCFG, 1);
            clk_wl     = nan(NCFG, 1);
            clk_fl     = nan(NCFG, 1);
            block_idx  = nan(NCFG, 1);
            po2_idx    = nan(NCFG, 1);
            eq_idx     = nan(NCFG, 1);
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
                eq_fl(ci)      = cfg.EqFL;
                eq_wl(ci)      = P.NIntBits + cfg.EqFL;
                clk_fl(ci)     = cfg.ClkFL;
                clk_wl(ci)     = P.NIntBits + cfg.ClkFL;
                block_idx(ci)  = cfg.block_idx;
                po2_idx(ci)    = cfg.po2_idx;
                eq_idx(ci)     = cfg.eq_idx;
                clk_idx(ci)    = cfg.clk_idx;
            end
            tbl = table(block_name, n_cd, n_overlap, n_aeq, po2, ...
                ki_col, kp_col, eq_wl, eq_fl, clk_wl, clk_fl, ...
                block_idx, po2_idx, eq_idx, clk_idx, ber, fec_snr);
            tbl.Properties.VariableNames{'ki_col'} = 'ki';
            tbl.Properties.VariableNames{'kp_col'} = 'kp';

            % --- Save ------------------------------------------------
            S.tbl         = tbl;
            S.SNR_dB_vec  = P.SNR_dB_vec;
            S.CFO_GHz     = P.CFO_GHz;
            S.NIntBits    = P.NIntBits;
            S.EqFL_vec    = P.EqFL_vec;
            S.ClkFL_vec   = P.ClkFL_vec;
            S.params      = P;

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
            P.EqFL_vec      = tc.EqFL_vec;
            P.ClkFL_vec     = tc.ClkFL_vec;
            P.FEC_BER       = tc.FEC_BER;
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
        end

        function ber = runOnePoint(P, cfg, rxSig, symbols, cfo)
            try
                [eqSym, cfoBinsApplied] = ...
                    combined_eq_clk_fxp_sweep.runFxpBlock(P, cfg, rxSig);
            catch ME
                fprintf(['    %s NCD=%d po2=%d  EQ FL=%d Clk FL=%d ', ...
                    '-> ERROR: %s\n'], ...
                    cfg.block, cfg.NCD, cfg.Po2Twiddle, ...
                    cfg.EqFL, cfg.ClkFL, ME.message);
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

            % Build per-section precision structs from the swept FL
            % entries (integer bits fixed at P.NIntBits).
            EqPrec  = struct('WL', P.NIntBits + cfg.EqFL,  'FL', cfg.EqFL);
            ClkPrec = struct('WL', P.NIntBits + cfg.ClkFL, 'FL', cfg.ClkFL);

            % Build T at runtime: cheap, and needed to know T.Static.x
            % (input cast) and T.AdaptEq.y (Pilots cast) for the MEX.
            switch cfg.block
                case 'cd_gardner_cma'
                    Tcfg = struct( ...
                        'Static',  EqPrec, ...
                        'Clk',     ClkPrec, ...
                        'AdaptEq', EqPrec);
                    T = eq_clk.combined_cd_fd_gardner_adaptive_fxp_types(Tcfg);
                case 'cd_godard_cma'
                    Tcfg = struct( ...
                        'Static',  EqPrec, ...
                        'Godard',  ClkPrec, ...
                        'AdaptEq', EqPrec);
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
                cfg.block, cfg.EqFL, cfg.ClkFL);
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

        function name = mexBaseName(blkName, eqFL, clkFL)
            % Unique MEX file name per (block, EqFL, ClkFL) combo.
            switch blkName
                case 'cd_gardner_cma'
                    base = 'combined_cd_fd_gardner_adaptive_fxp';
                case 'cd_godard_cma'
                    base = 'combined_cd_fd_godard_adaptive_fxp';
                otherwise
                    error('combined_eq_clk_fxp_sweep:badBlock', ...
                          'Unknown block: %s', blkName);
            end
            name = sprintf('%s_e%02dc%02d_mex', base, eqFL, clkFL);
        end

        function buildOneMex(P, cfgCoder, blkName, eqFL, clkFL, ...
                             mexBase, repoRoot)
            % Codegen one precision-specialised MEX into src/+eq_clk/.
            EqPrec  = struct('WL', P.NIntBits + eqFL,  'FL', eqFL);
            ClkPrec = struct('WL', P.NIntBits + clkFL, 'FL', clkFL);
            outPath = fullfile(repoRoot, 'src', '+eq_clk', mexBase);

            % Pick nominal NCD/NTaps/ki/kp from the first block design;
            % all of these are runtime args, only their *types* (double
            % scalars) are baked in by codegen.
            nCdProto = P.NCD(1);
            nOvProto = 2 * ceil((nCdProto - 1) / 2);

            switch blkName
                case 'cd_gardner_cma'
                    Tcfg = struct('Static',  EqPrec, ...
                                  'Clk',     ClkPrec, ...
                                  'AdaptEq', EqPrec);
                    T = eq_clk.combined_cd_fd_gardner_adaptive_fxp_types(Tcfg);
                    muVal = P.MuGardner;
                    n1Val = P.N1Gardner;
                case 'cd_godard_cma'
                    Tcfg = struct('Static',  EqPrec, ...
                                  'Godard',  ClkPrec, ...
                                  'AdaptEq', EqPrec);
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
