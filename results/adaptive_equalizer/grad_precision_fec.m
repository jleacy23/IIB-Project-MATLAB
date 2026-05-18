classdef grad_precision_fec < matlab.unittest.TestCase
%GRAD_PRECISION_FEC  FEC-SNR vs gradient-estimation precision sweep (MEX).
%
%   Sweeps the fixed-point precision of the CMA gradient estimate
%   (T.grad fractional bits) and the weight-update decimation
%   (UpdateStep) of src/+adaptive_eq/equalize_fxp.m, for both
%   SignOnly = false and SignOnly = true, and records the SNR that
%   achieves the FEC BER limit (2e-2) for every grid point and every
%   Monte-Carlo run.
%
%   Speed:  both equalize_fxp and viterbiViterbi_fxp are compiled to
%   MEX.  UpdateStep and SignOnly are runtime arguments, so only one
%   equalize MEX is built per GradFL value (the types table is baked
%   in at compile time).  The high-precision VV MEX is built once.
%
%   Signal chain
%     random unit-energy DP-QPSK  ->  SpS upsample
%       ->  PMD (constant realisation for every run)
%       ->  AWGN (SNR sweep)
%       ->  equalize_fxp_mex  (gradient precision swept)
%       ->  viterbiViterbi_fxp_mex  (wide fixed-point = high precision)
%       ->  QPSK pi/2 + polarisation ambiguity resolved by min BER
%       ->  BER  ->  FEC-SNR crossing (linear interpolation)
%
%   Grid (integer bits fixed at 16 for the gradient type)
%     GradFL_vec     = [2 4 6 8 10 12 14 16]
%     UpdateStep_vec = [1 2 4 16 32 64]
%     SignOnly       = {false, true}
%     MuScaling      = {fixed, scaled}  (scaled: Mu *= UpdateStep, so a
%                       decimated update of N uses N*Mu to keep the mean
%                       per-sample step size constant)
%     NTaps_vec      = [3 5 7]
%     L_km_vec       = [80 20]   (network configs: 80 km/1:16,
%                                 20 km/1:512 — see K_vec)
%
%   Output  ->  grad_precision_fec_sweep.mat
%     fecSNR  - [NGradFL x NUpdateStep x 2 x 2 x NNTaps x NL x NTrials]
%               FEC SNR [dB]  (3rd dim: 1 = SignOnly false,
%               2 = SignOnly true; 4th dim: 1 = fixed Mu,
%               2 = Mu scaled by UpdateStep; 5th dim: NTaps_vec;
%               6th dim: L_km_vec).
%               NaN where the BER curve never crosses the FEC limit.
%     berAll  - [NGradFL x NUpdateStep x 2 x 2 x NNTaps x NL x NTrials x NSNR]
%               raw BER.
%     params  - struct of the sweep parameter vectors.
%
%   Run with:
%       runtests('grad_precision_fec')
%       runtests('grad_precision_fec','ProcedureName','test_grid_sweep')
%
%   NOTE  This rebuilds src/+adaptive_eq/equalize_fxp_mex once per
%   GradFL and src/+carrier_recovery/viterbiViterbi_fxp_mex once,
%   overwriting any existing binaries (same pattern as bit_width_full).

    %% ================================================================
    %  Parameters
    %% ================================================================
    properties (Constant)

        % System
        N_pol   = 2
        Nsym    = 16000          % QPSK symbols per polarisation
        SpS     = 2             % samples per symbol
        Rs      = 30.5            % [GBd]

        % Monte-Carlo
        NTrials = 100             % runs used to recompute the FEC SNR

        % SNR sweep [dB]
        SNR_dB_vec = 0 : 2 : 20

        % Gradient-precision sweep — integer bits fixed at 16
        IntBits        = 16
        GradFL_vec     = [2, 6, 10]
        UpdateStep_vec = [1]

        % Shared equalise fimath + non-swept numerictypes.
        %   fi binary ops require equal fimath, and equalize_fxp mixes
        %   T.w/T.x and T.x/T.err/T.y, so every field shares one fimath;
        %   only T.grad's numerictype changes across the sweep.
        DpWL   = 48             % shared product/sum word length
        DpFL   = 32             % shared product/sum fraction length
        HiWL   = 32             % x / y / acc / err / R_CMA numerictype
        HiFL   = 24
        WWL    = 40             % weight numerictype (fine: mu*grad != 0)
        WFL    = 32

        % Adaptive EQ settings
        %   NTaps is swept (NTaps_vec).  The scalar NTaps is only the
        %   representative example arg for the equalize MEX build; it is
        %   a runtime argument, so any swept value works without rebuild.
        NTaps       = 3
        NTaps_vec   = [1, 3, 5]
        Mu          = 1e-3
        SingleSpike = true
        N1          = 4000
        NOut        = 8000
        PLanes      = 32            % parallel lanes for equalize_fxp

        % Channel — two network configs (fibre length / splitting ratio)
        L_km_vec = [80, 20]     % fibre lengths [km] (80 km / 20 km)
        K_vec    = [16, 512]    % 1:K split per L (1:16 / 1:512)
        DGDSpec  = 0.1          % PMD coeff [ps/sqrt(km)]
        N_pmd    = 1

        % High-precision Viterbi-Viterbi carrier recovery
        VV_NTaps       = 5
        VV_BlockLen    = 100
        VV_PilotThresh = 1e6                       % pilot correction off
        VV_FxpConfig   = struct('WL', 64, 'FL', 32)% wide => high precision
        CordicIts      = 16

        % FEC limit
        FEC_BER = 2e-2

        % Deterministic seeds
        SymSeed   = 1000        % per-trial symbol seed offset
        PmdSeed   = 777         % CONSTANT PMD realisation
        NoiseSeed = 5000        % per-trial / per-SNR noise seed offset

    end

    %% ================================================================
    %  Setup
    %% ================================================================
    methods (TestClassSetup)
        function setupPath(~)
            root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
            addpath(fullfile(root, 'src'));
            addpath(fullfile(root, 'build'));
        end

        function buildVV(testCase)
            % High-precision VV MEX, built once into src/+carrier_recovery.
            grad_precision_fec.buildVVMexInSrc(testCase);
        end
    end

    %% ================================================================
    %  Test
    %% ================================================================
    methods (Test)

        function test_grid_sweep(testCase)
            P   = testCase;
            NFL = numel(P.GradFL_vec);
            NUS = numel(P.UpdateStep_vec);
            NSO = 2;
            NMS = 2;                 % 1 = fixed Mu, 2 = Mu * UpdateStep
            NNT = numel(P.NTaps_vec);
            NL  = numel(P.L_km_vec);
            NTR = P.NTrials;
            NSN = numel(P.SNR_dB_vec);

            fecSNR = nan(NFL, NUS, NSO, NMS, NNT, NL, NTR);
            berAll = nan(NFL, NUS, NSO, NMS, NNT, NL, NTR, NSN);
            signOnlyVals = [false, true];

            T_vv = carrier_recovery.fxp_types(P.VV_FxpConfig);

            for fl = 1:NFL
                gradFL = P.GradFL_vec(fl);
                T_eq   = grad_precision_fec.buildEqTypes(P, gradFL);

                % One equalize MEX per GradFL (UpdateStep/SignOnly are
                % runtime args, so they need no rebuild).
                fprintf('=== building equalize MEX  GradFL = %2d  (%d/%d) ===\n', ...
                    gradFL, fl, NFL);
                grad_precision_fec.buildEqMexInSrc(P, T_eq);

                for so = 1:NSO
                    SignOnly = signOnlyVals(so);
                    for us = 1:NUS
                        UpdateStep = P.UpdateStep_vec(us);
                        for ms = 1:NMS
                            % ms == 1 : fixed step size.
                            % ms == 2 : step size scaled in proportion to
                            %           the update decimation, so an
                            %           UpdateStep of N uses N*Mu (with
                            %           N == 1 this equals the fixed case).
                            if ms == 1
                                muVal = P.Mu;
                            else
                                muVal = P.Mu * UpdateStep;
                            end
                            for nt = 1:NNT
                                ntapsVal = P.NTaps_vec(nt);
                                fprintf(['  GradFL=%2d SignOnly=%d ', ...
                                         'UpdateStep=%2d Mu=%.3g ', ...
                                         'NTaps=%d (so %d/%d, us %d/%d, ', ...
                                         'ms %d/%d, nt %d/%d)\n'], ...
                                    gradFL, SignOnly, UpdateStep, muVal, ...
                                    ntapsVal, so, NSO, us, NUS, ...
                                    ms, NMS, nt, NNT);
                                for nl = 1:NL
                                    lVal = P.L_km_vec(nl);
                                    for tr = 1:NTR
                                        ber = grad_precision_fec.snrSweep( ...
                                            P, T_eq, T_vv, UpdateStep, ...
                                            SignOnly, muVal, ntapsVal, ...
                                            lVal, tr);
                                        berAll(fl, us, so, ms, nt, nl, tr, :) = ber;
                                        fecSNR(fl, us, so, ms, nt, nl, tr) = ...
                                            grad_precision_fec.fecCrossing( ...
                                                P.SNR_dB_vec, ber, P.FEC_BER);
                                    end
                                end
                            end
                        end
                    end
                end
            end

            params = struct( ...
                'GradFL_vec',     P.GradFL_vec, ...
                'UpdateStep_vec', P.UpdateStep_vec, ...
                'SNR_dB_vec',     P.SNR_dB_vec, ...
                'SignOnly_dim',   {{'false', 'true'}}, ...
                'MuScaling_dim',  {{'fixed', 'scaled_by_UpdateStep'}}, ...
                'NTaps_vec',      P.NTaps_vec, ...
                'L_km_vec',       P.L_km_vec, ...
                'K_vec',          P.K_vec, ...
                'Mu',             P.Mu, ...
                'FEC_BER',        P.FEC_BER);             %#ok<NASGU>

            outDir  = fileparts(mfilename('fullpath'));
            outFile = fullfile(outDir, 'grad_precision_fec_sweep.mat');
            save(outFile, 'fecSNR', 'berAll', 'params');
            fprintf('Saved grid sweep to %s\n', outFile);

            testCase.verifyTrue(any(isfinite(fecSNR(:))), ...
                'No grid point ever crossed the FEC limit.');
        end

    end

    %% ================================================================
    %  Helpers
    %% ================================================================
    methods (Static, Access = private)

        % ---- Custom equalize_fxp types table -------------------------
        function T = buildEqTypes(P, gradFL)
            % Only T.grad's numerictype (fraction length) is swept; every
            % other field is high precision so the gradient store is the
            % only quantisation that changes across the sweep.  Integer
            % bits of the gradient type = 16  (WL = 16 + gradFL).  All
            % fields share one fimath (fi binary ops require it).
            F = fimath( ...
                'RoundingMethod',        'Floor', ...
                'OverflowAction',        'Wrap',  ...
                'ProductMode',           'SpecifyPrecision', ...
                'ProductWordLength',      P.DpWL, ...
                'ProductFractionLength',  P.DpFL, ...
                'SumMode',               'SpecifyPrecision', ...
                'SumWordLength',          P.DpWL, ...
                'SumFractionLength',      P.DpFL);

            gWL = P.IntBits + gradFL;
            T.x     = fi([], 1, P.HiWL, P.HiFL, F);
            T.w     = fi([], 1, P.WWL,  P.WFL,  F);
            T.y     = fi([], 1, P.HiWL, P.HiFL, F);
            T.acc   = fi([], 1, P.HiWL, P.HiFL, F);
            T.err   = fi([], 1, P.HiWL, P.HiFL, F);
            T.grad  = fi([], 1, gWL,    gradFL, F);
            T.R_CMA = fi([], 1, P.HiWL, P.HiFL, F);
        end

        % ---- MEX builds ---------------------------------------------
        function buildEqMexInSrc(P, T_eq)
            % codegen equalize_fxp directly with the custom types table
            % (the standard build script only supports uniform configs).
            % SignOnly / UpdateStep are passed as plain (non-constant)
            % scalars so the single MEX handles every runtime value.
            clear mex %#ok<CLMEX>
            srcDir  = fullfile(fileparts(mfilename('fullpath')), ...
                '..', '..', 'src');
            mexFile = fullfile(srcDir, '+adaptive_eq', ...
                ['equalize_fxp_mex.' mexext]);
            if isfile(mexFile), delete(mexFile); end

            x      = fi(complex(0, 0), numerictype(T_eq.x), fimath(T_eq.x));
            InType = coder.typeof(x, [Inf, 2], [true, false]);
            args   = { InType, ...
                double(P.SpS), double(P.NTaps), double(P.Mu), ...
                logical(P.SingleSpike), double(P.N1), double(P.NOut), ...
                logical(false), double(1), T_eq, double(P.PLanes) };

            cfg = coder.config('mex');
            cfg.GenerateReport = false;
            codegen('-config', cfg, 'adaptive_eq.equalize_fxp', ...
                '-args', args, ...
                '-o', fullfile(srcDir, '+adaptive_eq', 'equalize_fxp_mex'));
        end

        function buildVVMexInSrc(P)
            clear mex %#ok<CLMEX>
            B.N_pol          = P.N_pol;
            B.VV_NTaps       = P.VV_NTaps;
            B.PilotLen       = 1;
            B.BlockLen       = P.VV_BlockLen;
            B.StepSize       = P.VV_BlockLen;
            B.FxpConfig_VV   = P.VV_FxpConfig;
            B.PilotThreshold = P.VV_PilotThresh;
            B.CordicIts      = P.CordicIts;

            cfg = coder.config('mex');
            cfg.GenerateReport = false;
            build_carrier_recovery_viterbiViterbi_fxp_mex(B, cfg);
        end

        % ---- One SNR sweep for a given grid point / trial ------------
        function ber = snrSweep(P, T_eq, T_vv, UpdateStep, SignOnly, muVal, ntapsVal, lVal, trial)
            NSN = numel(P.SNR_dB_vec);
            ber = nan(1, NSN);

            % Symbols fixed within a trial; only the noise varies.
            rng(P.SymSeed + trial);
            syms  = grad_precision_fec.randomQPSK(P.Nsym, P.N_pol);
            txSig = repelem(syms, P.SpS, 1);

            for si = 1:NSN
                % --- Constant PMD realisation for every run ---
                rng(P.PmdSeed);
                rxSig = channel.add_pmd(txSig, lVal, P.SpS, ...
                    P.Rs, P.DGDSpec, P.N_pmd);

                % --- Per-trial / per-SNR AWGN ---
                rng(P.NoiseSeed + trial * 131 + si);
                rxSig = channel.add_awgn(rxSig, P.SNR_dB_vec(si));

                % --- Fixed-point adaptive equaliser (MEX) ---
                rxFi  = cast(rxSig, 'like', T_eq.x);
                eqSig = adaptive_eq.equalize_fxp_mex(rxFi, ...
                    double(P.SpS), double(ntapsVal), double(muVal), ...
                    logical(P.SingleSpike), double(P.N1), ...
                    double(P.NOut), logical(SignOnly), ...
                    double(UpdateStep), T_eq, double(P.PLanes));
                eqSig = double(eqSig);

                % --- High-precision Viterbi-Viterbi (MEX) ---
                crSig = grad_precision_fec.runVV(P, T_vv, eqSig, ...
                    P.SNR_dB_vec(si));

                % --- Reference symbols aligned to equaliser output ---
                L   = size(crSig, 1);
                ref = syms(P.NOut + 1 : P.NOut + L, :);

                % Drop VV edge transient before scoring.
                g = P.VV_NTaps;
                ber(si) = grad_precision_fec.berMinAmbiguity( ...
                    crSig(g+1:end-g, :), ref(g+1:end-g, :));
            end
        end

        % ---- High-precision VV (MEX) --------------------------------
        function v = runVV(P, T_vv, x, SNR_dB)
            Lflt = 2 * P.VV_NTaps + 1;
            wVV  = carrier_recovery.genVVFilter(1e3, P.Rs, SNR_dB, ...
                       1, P.N_pol, P.VV_NTaps);
            wVV  = wVV(:);
            if numel(wVV) ~= Lflt
                wVV = ones(Lflt, 1) / Lflt;
            end

            Nsym    = size(x, 1);
            NBlocks = ceil(Nsym / P.VV_BlockLen);

            x_fi   = cast(x,                 'like', T_vv.x);
            w_fi   = cast(wVV,               'like', T_vv.w);
            pil_fi = cast(zeros(NBlocks, P.N_pol), 'like', T_vv.x);

            [v, ~] = carrier_recovery.viterbiViterbi_fxp_mex( ...
                x_fi, double(P.N_pol), double(P.VV_NTaps), w_fi, ...
                pil_fi, double(P.VV_BlockLen), double(P.VV_BlockLen), ...
                double(P.VV_PilotThresh), double(P.CordicIts), T_vv);
            v = double(v);
        end

        % ---- Random unit-energy DP-QPSK -----------------------------
        function s = randomQPSK(Nsym, NPol)
            re = 2 * randi([0 1], Nsym, NPol) - 1;
            im = 2 * randi([0 1], Nsym, NPol) - 1;
            s  = (re + 1j * im) / sqrt(2);   % |s| = 1
        end

        % ---- BER with QPSK pi/2 + polarisation ambiguity ------------
        function ber = berMinAmbiguity(eqSym, refSyms)
            NPol    = size(refSyms, 2);
            totErr  = 0;
            totBits = 0;
            for p = 1:NPol
                refBits = modem.symbolsToBits(refSyms(:, p));
                best    = Inf;
                for q = 1:NPol
                    for k = 0:3
                        rot  = eqSym(:, q) .* exp(-1j * k * pi/2);
                        bits = modem.symbolsToBits( ...
                                   modem.decideSymbols(rot));
                        n    = min(numel(refBits), numel(bits));
                        e    = sum(refBits(1:n) ~= bits(1:n)) / n;
                        if e < best, best = e; end
                    end
                end
                totErr  = totErr  + best * numel(refBits);
                totBits = totBits + numel(refBits);
            end
            ber = totErr / totBits;
        end

        % ---- FEC-SNR crossing (linear interpolation) ----------------
        function xCross = fecCrossing(x, y, yLimit)
            xCross = NaN;
            n = numel(x);
            if n < 2, return; end

            exact = find(y == yLimit, 1, 'first');
            if ~isempty(exact)
                xCross = x(exact);
                return;
            end
            for i = 1:(n - 1)
                if (y(i) - yLimit) * (y(i+1) - yLimit) < 0
                    xCross = x(i) + (yLimit - y(i)) * ...
                        (x(i+1) - x(i)) / (y(i+1) - y(i));
                    return;
                end
            end
        end

    end
end
