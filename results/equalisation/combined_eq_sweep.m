classdef combined_eq_sweep < matlab.unittest.TestCase
%COMBINED_EQ_SWEEP  Grid sweep of CD + adaptive equaliser in cascade.
%
%   Cascades a fixed-point CD equaliser with a fixed-point adaptive
%   (butterfly) equaliser so the two dispersion mechanisms can be
%   compensated end-to-end, and sweeps the variant and precision of each.
%
%   Test parameters
%       1. CD equaliser variant      CDConfigs        {overlap_save,
%                                                      overlap_save_po2}
%                                                     (time_domain skipped:
%                                                      known not viable in
%                                                      cd_eq_precision_sweep)
%       2. CD fxp precision          CD_FL_vec        T.grad/twiddle FL
%       3. AEQ variant               AEQConfigs       (mode, sign_only)
%                                                     {CMA, pilot} x
%                                                     {direct, sign-sign}
%       4. AEQ fxp precision         AEQ_FL_vec       T.grad FL
%       5. AEQ tap count             AEQ_NTaps_vec    butterfly FIR length
%       6. Network configuration     {A, B}           20 km / 80 km
%
%   Channel chain (the receiver under test):
%       RRC pulse-shape (Tx) -> CD -> AWGN -> PMD ->
%       CD equaliser (fxp) -> adaptive equaliser (fxp) -> BER
%   No matched filter is applied between the equalisers, matching the
%   convention adopted in adaptive_eq_grid_sweep (the AEQ sees the raw
%   oversampled signal).
%
%   For each (cd_config, aeq_config, cd_fl, aeq_fl, aeq_ntaps, network)
%   point the test runs NTrials independent realisations across the SNR
%   sweep and stores the raw BER cube.
%
%   Output (saved next to this file):
%       combined_eq_sweep.mat  with variables:
%         tbl     - results table, one row per
%                     (cd_config, aeq_mode, aeq_sign_only,
%                      cd_fl, aeq_fl, aeq_ntaps, network):
%                     cd_config    (string)
%                     aeq_mode     (string)   'CMA' | 'pilot'
%                     aeq_sign_only(logical)  false = direct, true = sign-sign
%                     network      (string)   'A' | 'B'
%                     L_km         (double)
%                     splitting    (string)
%                     n_fft        (double)
%                     ncd_report   (double)
%                     aeq_ntaps    (double)
%                     cd_fl        (double)
%                     cd_wl        (double)
%                     aeq_fl       (double)
%                     aeq_wl       (double)
%                     ber          (cell)     [NTrials x NSNR] raw BER
%         SNR_dB_vec - SNR sweep axis [dB]
%         params     - copy of the run parameters
%
%   Run with:
%       runtests('combined_eq_sweep')

    %% ================================================================
    %  Parameters
    %% ================================================================
    properties (Constant)

        % System (report tab:network_params)
        Rs        = 30.5
        D         = 17
        CWL       = 1550
        N_pol     = 2
        SpS       = 2
        DGDSpec   = 0.1
        N_pmd     = 1

        % Pulse shaping (Nyquist / raised-cosine)
        Rolloff   = 0.25
        Span      = 10

        % Monte-Carlo
        NTrials   = 10
        NSub      = 8

        % SNR sweep [dB]
        SNR_dB_vec = 0 : 2 : 30

        % Bit-width sweeps (integer bits fixed; WL = IntBits + FL)
        IntBits    = 16
        CD_FL_vec  = [4]
        AEQ_FL_vec = [2, 4, 6, 8, 10]

        % CD equaliser variants (time_domain dropped: not viable, see
        % cd_eq_precision_sweep results).
        CDConfigs  = {'overlap_save', 'overlap_save_po2'}

        % Adaptive equaliser variants: rows of {mode, sign_only}.
        AEQConfigs = { ...
            'CMA',   true;}

        % AEQ tap count sweep
        AEQ_NTaps_vec = [1, 3, 5]

        % AEQ step size and re-init iteration per (mode, ntaps, sign_only),
        % taken from adaptive_eq_grid_sweep.Mu_vec/N1_vec.  Fields keyed
        % by mode; each field is [NTaps x NSignOnly] (col 1 = direct,
        % col 2 = sign-sign; rows align with AEQ_NTaps_vec).
        %                                  direct     sign-sign
        Mu = struct( ...
            'CMA',   [9.766e-4,  9.766e-4;     % NTaps = 1
                      6.104e-5,  9.766e-4;     % NTaps = 3
                      6.104e-5,  9.766e-4], ... % NTaps = 5
            'pilot', [9.766e-4,  9.766e-4;     % NTaps = 1
                      9.766e-4,  9.766e-4;     % NTaps = 3
                      9.766e-4,  9.766e-4])    % NTaps = 5
        N1 = struct( ...
            'CMA',   [2000, 1000;
                      2000, 1000;
                      2000, 1000], ...
            'pilot', [2000, 2000;
                      2000, 2000;
                      2000, 2000])
        SingleSpike = true
        NOut        = 4000
        UpdateStep  = 1

    end

    %% ================================================================
    %  Setup
    %% ================================================================
    methods (TestClassSetup)
        function setupPaths(~)
            here = fileparts(mfilename('fullpath'));
            addpath(genpath(fullfile(here, '..', '..', 'src')));
            addpath(fullfile(here, '..', '..', 'build'));
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

        function test_combined_sweep(testCase)
            P    = combined_eq_sweep.extractParams(testCase);
            nets = combined_eq_sweep.networks();
            cdCfgs   = testCase.CDConfigs;
            aeqCfgs  = testCase.AEQConfigs;
            cdFLs    = testCase.CD_FL_vec;
            aeqFLs   = testCase.AEQ_FL_vec;
            aeqTaps  = testCase.AEQ_NTaps_vec;
            mus      = testCase.Mu;
            n1s      = testCase.N1;
            NSNR     = numel(P.SNR_dB_vec);

            NCD_FL   = numel(cdFLs);
            NAEQ_FL  = numel(aeqFLs);
            NCD_CFG  = numel(cdCfgs);
            NAEQ     = size(aeqCfgs, 1);
            NAEQ_TAP = numel(aeqTaps);
            NNET     = numel(nets);
            nRows    = NCD_FL * NAEQ_FL * NCD_CFG * NAEQ * NAEQ_TAP * NNET;

            for qi = 1:NAEQ
                m = aeqCfgs{qi, 1};
                assert(isfield(mus, m) && isfield(n1s, m), ...
                       'Mu/N1 must have a field for AEQ mode %s.', m);
                assert(size(mus.(m), 1) == NAEQ_TAP && ...
                       size(n1s.(m), 1) == NAEQ_TAP, ...
                       'Mu.%s / N1.%s must have NTaps rows.', m, m);
            end

            % Pre-allocate table columns
            cd_config     = strings(nRows, 1);
            aeq_mode      = strings(nRows, 1);
            aeq_sign_only = false(nRows, 1);
            network       = strings(nRows, 1);
            L_km          = nan(nRows, 1);
            splitting     = strings(nRows, 1);
            n_fft         = nan(nRows, 1);
            ncd_report    = nan(nRows, 1);
            aeq_ntaps     = nan(nRows, 1);
            cd_fl         = nan(nRows, 1);
            cd_wl         = nan(nRows, 1);
            aeq_fl        = nan(nRows, 1);
            aeq_wl        = nan(nRows, 1);
            ber           = cell(nRows, 1);

            row = 0;
            for cdi = 1:NCD_FL
                cdFLv  = cdFLs(cdi);
                cdFxp  = struct('WL', P.IntBits + cdFLv, 'FL', cdFLv);
                fprintf('=== CD FL = %2d  (WL = %2d)  [%d/%d] : building CD MEX ===\n', ...
                        cdFLv, P.IntBits + cdFLv, cdi, NCD_FL);
                combined_eq_sweep.buildCDMex(P, cdFxp);
                T_fd = cd_eq.equalize_fxp_types(cdFxp);

                for ai = 1:NAEQ_FL
                    aeqFLv = aeqFLs(ai);
                    aeqFxp = struct('WL', P.IntBits + aeqFLv, 'FL', aeqFLv);
                    fprintf(['--- AEQ FL = %2d  (WL = %2d)  [%d/%d] : ', ...
                             'building AEQ MEX ---\n'], ...
                            aeqFLv, P.IntBits + aeqFLv, ai, NAEQ_FL);
                    combined_eq_sweep.buildAEQMex(P, aeqFxp);
                    T_aeq = adaptive_eq.equalize_fxp_types(aeqFxp);

                    for ni = 1:NNET
                        net = nets(ni);

                        % Accumulate BER over (cd_cfg, aeq_cfg, ntaps) so
                        % each channel realisation drives all combinations.
                        berCells = cell(NCD_CFG, NAEQ, NAEQ_TAP);
                        for ci = 1:NCD_CFG
                            for qi = 1:NAEQ
                                for ti = 1:NAEQ_TAP
                                    berCells{ci, qi, ti} = nan(P.NTrials, NSNR);
                                end
                            end
                        end

                        for tr = 1:P.NTrials
                            fprintf(['    CDFL=%2d AEQFL=%2d net=%s  ', ...
                                     'trial %d/%d\n'], ...
                                    cdFLv, aeqFLv, net.name, tr, P.NTrials);
                            for si = 1:NSNR
                                snr = P.SNR_dB_vec(si);
                                [rxSig, symbols, PilotsAll] = ...
                                    combined_eq_sweep.buildChannel( ...
                                        P, net, snr, tr);

                                for ci = 1:NCD_CFG
                                    cdEq = combined_eq_sweep.runCDEq( ...
                                        cdCfgs{ci}, rxSig, P, net, T_fd);

                                    for qi = 1:NAEQ
                                        mode     = aeqCfgs{qi, 1};
                                        signFlag = aeqCfgs{qi, 2};
                                        col      = 1 + double(signFlag);
                                        for ti = 1:NAEQ_TAP
                                            nt    = aeqTaps(ti);
                                            Pr    = P;
                                            Pr.Mu = mus.(mode)(ti, col);
                                            Pr.N1 = n1s.(mode)(ti, col);
                                            eqSym = combined_eq_sweep.runAEQ( ...
                                                mode, cdEq, PilotsAll, ...
                                                nt, signFlag, Pr, T_aeq);
                                            berCells{ci, qi, ti}(tr, si) = ...
                                                combined_eq_sweep.computeBER( ...
                                                    eqSym, symbols, P.NOut);
                                        end
                                    end
                                end
                            end
                        end

                        for ci = 1:NCD_CFG
                            for qi = 1:NAEQ
                                for ti = 1:NAEQ_TAP
                                    row = row + 1;
                                    cd_config(row)     = string(cdCfgs{ci});
                                    aeq_mode(row)      = string(aeqCfgs{qi, 1});
                                    aeq_sign_only(row) = aeqCfgs{qi, 2};
                                    network(row)       = net.name;
                                    L_km(row)          = net.L;
                                    splitting(row)     = net.split;
                                    n_fft(row)         = net.NFFT;
                                    ncd_report(row)    = net.Ncd;
                                    aeq_ntaps(row)     = aeqTaps(ti);
                                    cd_fl(row)         = cdFLv;
                                    cd_wl(row)         = P.IntBits + cdFLv;
                                    aeq_fl(row)        = aeqFLv;
                                    aeq_wl(row)        = P.IntBits + aeqFLv;
                                    ber{row}           = berCells{ci, qi, ti};
                                end
                            end
                        end
                    end
                end
            end

            tbl = table(cd_config, aeq_mode, aeq_sign_only, network, ...
                        L_km, splitting, n_fft, ncd_report, aeq_ntaps, ...
                        cd_fl, cd_wl, aeq_fl, aeq_wl, ber);

            SNR_dB_vec = P.SNR_dB_vec;  %#ok<NASGU,PROP>
            params     = P;             %#ok<NASGU>

            outFile = fullfile(fileparts(mfilename('fullpath')), ...
                               'combined_eq_sweep.mat');
            save(outFile, 'tbl', 'SNR_dB_vec', 'params');
            fprintf('Saved combined EQ sweep results to %s\n', outFile);
            disp(tbl(:, {'cd_config','aeq_mode','aeq_sign_only', ...
                         'aeq_ntaps','network','cd_fl','aeq_fl'}));

            testCase.verifyEqual(height(tbl), nRows);
        end

    end

    %% ================================================================
    %  Helpers
    %% ================================================================
    methods (Static)

        function P = extractParams(tc)
            P.Rs          = tc.Rs;
            P.D           = tc.D;
            P.CWL         = tc.CWL;
            P.N_pol       = tc.N_pol;
            P.SpS         = tc.SpS;
            P.DGDSpec     = tc.DGDSpec;
            P.N_pmd       = tc.N_pmd;
            P.Rolloff     = tc.Rolloff;
            P.Span        = tc.Span;
            P.NTrials     = tc.NTrials;
            P.NSub        = tc.NSub;
            P.SNR_dB_vec  = tc.SNR_dB_vec;
            P.IntBits     = tc.IntBits;
            P.CD_FL_vec   = tc.CD_FL_vec;
            P.AEQ_FL_vec  = tc.AEQ_FL_vec;
            P.CDConfigs    = tc.CDConfigs;
            P.AEQConfigs   = tc.AEQConfigs;
            P.AEQ_NTaps_vec = tc.AEQ_NTaps_vec;
            % Scalar Mu / N1 are runtime-only build prototypes for codegen;
            % the test loop overrides them per (mode, ntaps, sign_only)
            % from Mu/N1.
            P.Mu          = tc.Mu.CMA(1, 1);
            P.N1          = tc.N1.CMA(1, 1);
            P.SingleSpike = tc.SingleSpike;
            P.NOut        = tc.NOut;
            P.UpdateStep  = tc.UpdateStep;
        end

        function nets = networks()
            % CPON networks (report tab:cd_taps / tab:cd_cost_eval): same
            % L / split as the adaptive sweep with the CD FFT/tap counts
            % attached so the CD MEX has all the info it needs.
            nets = struct( ...
                'name',  {"A", "B"}, ...
                'L',     {20, 80}, ...
                'split', {"1:512", "1:16"}, ...
                'Ncd',   {6, 21}, ...
                'NFFT',  {32, 128});
        end

        function buildCDMex(P, fxp)
            % Build both CD MEX binaries for the given fixed-point config.
            % NFFT, L and po2Twiddle are runtime arguments.
            clear mex %#ok<CLMEX>

            B.FxpConfig_CD = fxp;
            B.N_pol        = P.N_pol;
            B.D            = P.D;
            B.L            = 80;        % prototype (runtime-variable)
            B.CWL          = P.CWL;
            B.Rs           = P.Rs;
            B.SpS          = P.SpS;
            B.NFFT         = 128;       % prototype (runtime-variable)
            B.po2Twiddle   = false;     % prototype (runtime-variable)

            cfg = coder.config('mex');
            cfg.GenerateReport = false;
            build_cd_eq_equalize_fxp_mex(B, cfg);
        end

        function buildAEQMex(P, fxp)
            % Build the adaptive-EQ MEX for the given fixed-point config.
            % NTaps, Mode, Pilots, PLanes and BlockLen are runtime arguments.
            clear mex %#ok<CLMEX>

            B.FxpConfig_AEQ   = fxp;
            B.SpS             = P.SpS;
            B.AEQ_NTaps       = max(P.AEQ_NTaps_vec);   % prototype (runtime-variable)
            B.AEQ_Mu          = P.Mu;
            B.AEQ_SingleSpike = P.SingleSpike;
            B.AEQ_N1          = P.N1;
            B.AEQ_NOut        = P.NOut;
            B.AEQ_SignOnly    = false;
            B.AEQ_UpdateStep  = P.UpdateStep;
            B.AEQ_PLanes      = 32;     % prototype (runtime-variable)
            B.AEQ_Mode        = 0;
            B.AEQ_BlockLen    = 32;     % prototype (runtime-variable)

            cfg = coder.config('mex');
            cfg.GenerateReport = false;
            build_adaptive_eq_equalize_fxp_mex(B, cfg);
        end

        function [rxSig, symbols, PilotsAll] = buildChannel(P, net, SNR_dB, trialSeed)
            % CPON-framed DP-QPSK through CD + AWGN + PMD for one
            % (network, SNR) realisation.  No matched filter (the AEQ
            % sees the raw oversampled signal, matching the convention in
            % adaptive_eq_grid_sweep).
            rng(1000 * trialSeed + round(SNR_dB) + 7 * double(net.L));

            DATA_PER_SUBFRAME = 3586;
            nBits = P.NSub * DATA_PER_SUBFRAME * P.N_pol * 2;
            bits  = randi([0 1], nBits, 1);

            [symbols, pilots, ~, nSub] = modem.modulate(bits);
            PilotsAll = repmat(pilots, nSub, 1);

            txSig = modem.rrcPulse(symbols, P.SpS, P.Rolloff, P.Span);
            rxSig = channel.add_chromatic_dispersion(txSig, net.L, P.SpS, ...
                        P.Rs, P.D, P.CWL);
            rxSig = channel.add_awgn(rxSig, SNR_dB);
            rxSig = channel.add_pmd(rxSig, net.L, P.SpS, P.Rs, ...
                        P.DGDSpec, P.N_pmd);
        end

        function eqSig = runCDEq(cfg, rxSig, P, net, T_fd)
            % Run the selected CD equaliser (MEX) at oversampled rate.
            rx_fi = cast(rxSig, 'like', T_fd.x);
            switch cfg
                case 'overlap_save'
                    eqSig = cd_eq.equalize_fxp_mex(rx_fi, double(P.D), ...
                        double(net.L), double(P.CWL), double(P.Rs), ...
                        double(P.N_pol), double(P.SpS), double(net.NFFT), ...
                        false, T_fd);
                case 'overlap_save_po2'
                    eqSig = cd_eq.equalize_fxp_mex(rx_fi, double(P.D), ...
                        double(net.L), double(P.CWL), double(P.Rs), ...
                        double(P.N_pol), double(P.SpS), double(net.NFFT), ...
                        true, T_fd);
                otherwise
                    error('combined_eq_sweep:badCDConfig', ...
                          'Unknown CD config: %s', cfg);
            end
            eqSig = double(eqSig);
        end

        function eqSym = runAEQ(mode, rxSig, PilotsAll, ntaps, signOnly, P, T)
            % Per-mode PLanes/BlockLen following adaptive_eq_grid_sweep.runAEQ:
            %   CMA   -> parallel-lane CMA, Mode=0, PLanes=32, BlockLen=32
            %   pilot -> data-aided LMS,   Mode=1, PLanes=1,  BlockLen=32
            rx_fi = cast(rxSig, 'like', T.x);
            switch mode
                case 'CMA'
                    pilots_fi = cast(complex(zeros(0, 2)), 'like', T.y);
                    Mode      = 0;
                    pLanes    = 32;
                    blockLen  = 32;
                    subBlocks = 116;
                case 'pilot'
                    pilots_fi = cast(PilotsAll, 'like', T.y);
                    Mode      = 1;
                    pLanes    = 1;
                    blockLen  = 32;
                    subBlocks = 0;
                otherwise
                    error('combined_eq_sweep:badMode', ...
                          'Unknown AEQ mode: %s', mode);
            end

            eq = adaptive_eq.equalize_fxp_mex(rx_fi, double(P.SpS), ...
                double(ntaps), double(P.Mu), logical(P.SingleSpike), ...
                double(P.N1), double(P.NOut), logical(signOnly), ...
                double(P.UpdateStep), T, double(pLanes), ...
                double(Mode), pilots_fi, double(blockLen), double(subBlocks));
            eqSym = double(eq);
        end

        function BER = computeBER(eqSym, symbols, NOut)
            % BER vs TX symbols (offset by NOut), resolving residual phase
            % by a dense rotation search (32 angles of pi/16) plus a
            % polarisation-swap loop -- mirrors adaptive_eq_grid_sweep.
            nlen   = size(eqSym, 1);
            refEnd = min(NOut + nlen, size(symbols, 1));
            ref    = symbols(NOut+1 : refEnd, :);
            m      = min(size(eqSym, 1), size(ref, 1));
            if m < 1, BER = NaN; return; end
            eqSym  = eqSym(1:m, :);
            ref    = ref(1:m, :);

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
