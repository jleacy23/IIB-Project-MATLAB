classdef grad_precision_fec < matlab.unittest.TestCase
%GRAD_PRECISION_FEC  FEC-SNR vs gradient-estimation precision sweep.
%
%   Sweeps the fixed-point precision of the CMA gradient estimate
%   (T.grad fractional bits) and the weight-update decimation
%   (UpdateStep) of src/+adaptive_eq/equalize_fxp.m, for both
%   SignOnly = false and SignOnly = true, and records the SNR that
%   achieves the FEC BER limit (2e-2) for every grid point and every
%   Monte-Carlo run.
%
%   Signal chain
%     random unit-energy DP-QPSK  ->  SpS upsample
%       ->  PMD (constant realisation for every run)
%       ->  AWGN (SNR sweep)
%       ->  equalize_fxp  (fixed-point, gradient precision swept)
%       ->  viterbiViterbi_fxp  (double config = high precision)
%       ->  QPSK pi/2 + polarisation ambiguity resolved by min BER
%       ->  BER  ->  FEC-SNR crossing (linear interpolation)
%
%   Grid (integer bits fixed at 16 for the gradient type)
%     GradFL_vec    = [2 4 6 8 10 12 14 16]
%     UpdateStep_vec= [1 2 4 16 32 64]
%     SignOnly      = {false, true}
%
%   Output  ->  grad_precision_fec_sweep.mat
%     fecSNR  - [NGradFL x NUpdateStep x 2 x NTrials] FEC SNR [dB]
%               (3rd dim: 1 = SignOnly false, 2 = SignOnly true).
%               NaN where the BER curve never crosses the FEC limit.
%     berAll  - [NGradFL x NUpdateStep x 2 x NTrials x NSNR] raw BER,
%               kept for further processing.
%     plus the parameter vectors used.
%
%   Run with:
%       runtests('grad_precision_fec')
%       runtests('grad_precision_fec','ProcedureName','test_grid_sweep')
%
%   NOTE  equalize_fxp is executed as MATLAB fixed-point (no MEX), so the
%   full default grid is long-running.  Shrink Nsym / SNR_dB_vec /
%   NTrials for a quick smoke run.

    %% ================================================================
    %  Parameters
    %% ================================================================
    properties (Constant)

        % System
        N_pol   = 2
        Nsym    = 8000          % QPSK symbols per polarisation
        SpS     = 2             % samples per symbol
        Rs      = 32            % [GBd]

        % Monte-Carlo
        NTrials = 3             % runs used to recompute the FEC SNR

        % SNR sweep [dB]
        SNR_dB_vec = 0 : 1 : 30

        % Gradient-precision sweep — integer bits fixed at 16
        IntBits      = 16
        GradFL_vec   = [2, 4, 6, 8, 10, 12, 14, 16]
        UpdateStep_vec = [1, 2, 4, 16, 32, 64]

        % High-precision (non-swept) fixed-point fields.
        %   All fields share ONE fimath (fi binary ops require equal
        %   fimath); only the per-field numerictype WL/FL differs.  The
        %   shared product/sum precision is kept very fine so the only
        %   quantisation that varies across the sweep is T.grad's store.
        DpWL   = 48             % shared product/sum word length
        DpFL   = 32             % shared product/sum fraction length
        HiWL   = 32             % x / y / acc / err / R_CMA numerictype
        HiFL   = 24
        WWL    = 40             % weight numerictype (fine: mu*grad != 0)
        WFL    = 32

        % Adaptive EQ settings
        NTaps       = 3
        Mu          = 1e-3
        SingleSpike = true
        N1          = 1000
        NOut        = 2000

        % Channel
        L_km    = 80            % fibre length [km]
        DGDSpec = 0.5           % PMD coeff [ps/sqrt(km)]
        N_pmd   = 5

        % High-precision Viterbi-Viterbi carrier recovery
        VV_NTaps       = 5
        VV_BlockLen    = 100
        VV_PilotThresh = 1e6    % huge -> pilot correction disabled

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
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src'));
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
            NTR = P.NTrials;
            NSN = numel(P.SNR_dB_vec);

            fecSNR = nan(NFL, NUS, NSO, NTR);
            berAll = nan(NFL, NUS, NSO, NTR, NSN);

            signOnlyVals = [false, true];

            for so = 1:NSO
                SignOnly = signOnlyVals(so);
                for fl = 1:NFL
                    gradFL = P.GradFL_vec(fl);
                    T_eq   = grad_precision_fec.buildEqTypes(P, gradFL);
                    for us = 1:NUS
                        UpdateStep = P.UpdateStep_vec(us);
                        fprintf(['[SignOnly=%d] GradFL=%2d  UpdateStep=%2d ', ...
                                 '(so %d/%d, fl %d/%d, us %d/%d)\n'], ...
                            SignOnly, gradFL, UpdateStep, so, NSO, ...
                            fl, NFL, us, NUS);

                        for tr = 1:NTR
                            ber = grad_precision_fec.snrSweep( ...
                                P, T_eq, UpdateStep, SignOnly, tr);
                            berAll(fl, us, so, tr, :) = ber;
                            fecSNR(fl, us, so, tr) = ...
                                grad_precision_fec.fecCrossing( ...
                                    P.SNR_dB_vec, ber, P.FEC_BER);
                        end
                    end
                end
            end

            params = struct( ...
                'GradFL_vec',     P.GradFL_vec, ...
                'UpdateStep_vec', P.UpdateStep_vec, ...
                'SNR_dB_vec',     P.SNR_dB_vec, ...
                'SignOnly_dim',   {{'false', 'true'}}, ...
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
            % bits of the gradient type = 16  (WL = 16 + gradFL).
            %
            % All fields MUST share one fimath: fi binary operations
            % require equal fimath, and equalize_fxp mixes T.w/T.x in the
            % butterfly and T.x/T.err/T.y in the gradient.  The shared
            % product/sum precision is kept very fine (DpFL) so it never
            % limits the result before the per-field numerictype store.
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

        % ---- One SNR sweep for a given grid point / trial ------------
        function ber = snrSweep(P, T_eq, UpdateStep, SignOnly, trial)
            NSN = numel(P.SNR_dB_vec);
            ber = nan(1, NSN);

            % Symbols are fixed within a trial so only the noise varies
            % across the SNR sweep.
            rng(P.SymSeed + trial);
            syms  = grad_precision_fec.randomQPSK(P.Nsym, P.N_pol);
            txSig = repelem(syms, P.SpS, 1);

            for si = 1:NSN
                % --- Constant PMD realisation for every run ---
                rng(P.PmdSeed);
                rxSig = channel.add_pmd(txSig, P.L_km, P.SpS, ...
                    P.Rs, P.DGDSpec, P.N_pmd);

                % --- Per-trial / per-SNR AWGN ---
                rng(P.NoiseSeed + trial * 131 + si);
                rxSig = channel.add_awgn(rxSig, P.SNR_dB_vec(si));

                % --- Fixed-point adaptive equaliser ---
                rxFi  = cast(rxSig, 'like', T_eq.x);
                eqSig = adaptive_eq.equalize_fxp(rxFi, ...
                    double(P.SpS), double(P.NTaps), double(P.Mu), ...
                    logical(P.SingleSpike), double(P.N1), ...
                    double(P.NOut), logical(SignOnly), ...
                    double(UpdateStep), T_eq);
                eqSig = double(eqSig);

                % --- High-precision Viterbi-Viterbi carrier recovery ---
                crSig = grad_precision_fec.runVV(P, eqSig, ...
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

        % ---- High-precision VV (double config) ----------------------
        function v = runVV(P, x, SNR_dB)
            T_vv  = carrier_recovery.fxp_types('double');
            Lflt  = 2 * P.VV_NTaps + 1;
            wVV   = carrier_recovery.genVVFilter(1e3, P.Rs, SNR_dB, ...
                        1, P.N_pol, P.VV_NTaps);
            wVV   = wVV(:);
            if numel(wVV) ~= Lflt
                wVV = ones(Lflt, 1) / Lflt;
            end

            Nsym    = size(x, 1);
            NBlocks = ceil(Nsym / P.VV_BlockLen);
            pilots  = zeros(NBlocks, P.N_pol);   % unused: correction off

            [v, ~] = carrier_recovery.viterbiViterbi_fxp( ...
                x, P.N_pol, double(P.VV_NTaps), wVV, pilots, ...
                double(P.VV_BlockLen), double(P.VV_BlockLen), ...
                double(P.VV_PilotThresh), 16, T_vv);
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
            % Try all polarisation assignments and all multiples of pi/2;
            % pick the minimum BER per reference polarisation.
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
