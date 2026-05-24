classdef adaptive_eq_grid_sweep < matlab.unittest.TestCase
%ADAPTIVE_EQ_GRID_SWEEP  Grid sweep of the adaptive (butterfly) equaliser.
%
%   Sweeps, in a full Cartesian grid, the five levers of the fixed-point
%   adaptive equaliser and logs the raw per-iteration BER:
%
%       1. Gradient precision      FL_vec       (T.grad fractional bits,
%                                                WL = IntBits + FL).  The
%                                                rest of the data path
%                                                (T.x/y/acc/err/R_CMA) is
%                                                held at high precision and
%                                                T.w is pinned at 16 FL;
%                                                only the CMA gradient
%                                                store varies (mirrors the
%                                                working grad_precision_fec
%                                                pattern).
%       2. Number of taps          NTaps_vec
%       3. Update mode             {CMA, pilot}  (Mode 0 parallel-lane CMA
%                                                with PLanes=BlockLen=32
%                                                vs Mode 1 pilot-aided LMS
%                                                with PLanes=1, BlockLen=32)
%       4. Update strategy         SignOnly_vec  ({direct, sign-sign} weight
%                                                update — false = direct,
%                                                true = sign-sign)
%       5. Network configuration   {A, B}        (20 km / 80 km, report
%                                                \cref{chapter:energy})
%
%   The channel is CPON-framed DP-QPSK through AWGN + PMD with the
%   RRC + RRC matched-filter Nyquist split.  Both networks share the same
%   fibre PMD coefficient (D_PMD = 0.1 ps/sqrt(km)) but differ in length.
%
%   For each (FL, NTaps, mode, network) point the test runs NTrials
%   independent realisations across the SNR sweep and stores the raw BER.
%   No FEC-SNR threshold is computed — downstream post-processing decides
%   what to derive from the BER cube.
%
%   Output (saved next to this file):
%       adaptive_eq_grid_sweep.mat  with variables:
%         tbl     — results table, one row per
%                     (mode, network, NTaps, FL, sign_only):
%                     mode      (string)   'CMA' | 'pilot'
%                     network   (string)   'A' | 'B'
%                     L_km      (double)
%                     splitting (string)
%                     ntaps     (double)
%                     fl        (double)
%                     wl        (double)
%                     sign_only (logical)  false = direct, true = sign-sign
%                     ber       (cell)     [NTrials x NSNR] raw BER
%         SNR_dB_vec — SNR sweep axis [dB]
%         params     — copy of the run parameters
%
%   Run with:
%       runtests('adaptive_eq_grid_sweep')

    %% ================================================================
    %  Parameters
    %% ================================================================
    properties (Constant)

        % System (report tab:network_params)
        Rs        = 30.5         % symbol / channel bandwidth [GBd]
        CWL       = 1550         % central wavelength [nm]
        N_pol     = 2
        SpS       = 2            % samples per symbol
        DGDSpec   = 0.1          % PMD coefficient [ps/sqrt(km)]
        N_pmd     = 1           % PMD birefringent sections

        % Pulse shaping (Nyquist / raised-cosine)
        Rolloff   = 0.25
        Span      = 10

        % Monte-Carlo
        NTrials   = 20
        NSub      = 8            % CPON subframes per realisation

        % SNR sweep [dB]
        SNR_dB_vec = 0 : 2 : 30

        % Bit-width sweep — integer bits fixed; WL = IntBits + FL
        IntBits   = 16
        FL_vec    = [2, 4, 6, 8, 10, 12]

        % Tap-count sweep
        NTaps_vec = [1, 3, 5]

        % Update modes
        Modes        = {'CMA', 'pilot'}

        % Weight-update strategy:
        %   false = direct (full multiplications)
        %   true  = sign-sign (signum of error and complex-signum of output)
        SignOnly_vec = [false, true]

        % Adaptive-EQ fixed settings
        %   PLanes / BlockLen are set per mode in runAEQ
        %   (CMA -> 32/32 parallel; pilot -> 1/32 LMS).
        %
        % Mu_vec / N1_vec are [NTaps x NSignOnly] matrices so each
        % (tap length, update strategy) pair can use its own step size
        % and re-init iteration.  Columns align with SignOnly_vec
        % (col 1 = direct, col 2 = sign-sign); rows align with NTaps_vec.
        % adaptive_eq_convergence_sweep populates the best (Mu, N1) for
        % each tap length; update these matrices after running that sweep.
        %               direct     sign-sign
        Mu_vec      = [9.766e-4,  9.766e-4;   % NTaps = 1
                       6.104e-5,  9.766e-4;   % NTaps = 3
                       6.104e-5,  9.766e-4]   % NTaps = 5
        N1_vec      = [2000,  1000;
                       2000,  1000;
                       2000,  1000]
        SingleSpike = true
        NOut        = 4000       % discarded transient symbols
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

        function test_grid_sweep(testCase)
            P     = adaptive_eq_grid_sweep.extractParams(testCase);
            nets  = adaptive_eq_grid_sweep.networks();
            modes = testCase.Modes;
            signs = testCase.SignOnly_vec;
            FLs   = testCase.FL_vec;
            taps  = testCase.NTaps_vec;
            mus   = testCase.Mu_vec;
            n1s   = testCase.N1_vec;
            NSNR  = numel(P.SNR_dB_vec);

            assert(size(mus,1) == numel(taps) && size(mus,2) == numel(signs) && ...
                   size(n1s,1) == numel(taps) && size(n1s,2) == numel(signs), ...
                   'Mu_vec and N1_vec must be [NTaps x NSignOnly].');

            NFL   = numel(FLs);
            NNET  = numel(nets);
            NMODE = numel(modes);
            NTAP  = numel(taps);
            NSIGN = numel(signs);
            nRows = NFL * NNET * NMODE * NTAP * NSIGN;

            mode_c    = strings(nRows, 1);
            network   = strings(nRows, 1);
            L_km      = nan(nRows, 1);
            splitting = strings(nRows, 1);
            ntaps_c   = nan(nRows, 1);
            fl        = nan(nRows, 1);
            wl        = nan(nRows, 1);
            sign_only = false(nRows, 1);
            ber       = cell(nRows, 1);

            row = 0;
            for bi = 1:NFL
                flv = FLs(bi);
                fxp = struct('WL', P.IntBits + flv, 'FL', flv);

                fprintf('=== FL = %2d  (WL = %2d)  [%d/%d] : building AEQ MEX ===\n', ...
                        flv, P.IntBits + flv, bi, NFL);
                adaptive_eq_grid_sweep.buildAEQMex(P, fxp);
                T = adaptive_eq.equalize_fxp_types(fxp);

                for ni = 1:NNET
                    net = nets(ni);

                    % BER accumulators for this (FL, net): indexed
                    % {mode, ntaps, sign_only} -> [NTrials x NSNR]
                    berCells = cell(NMODE, NTAP, NSIGN);
                    for mi = 1:NMODE
                        for ti = 1:NTAP
                            for so = 1:NSIGN
                                berCells{mi, ti, so} = nan(P.NTrials, NSNR);
                            end
                        end
                    end

                    for tr = 1:P.NTrials
                        fprintf('    FL=%2d net=%s  trial %d/%d\n', ...
                                flv, net.name, tr, P.NTrials);
                        for si = 1:NSNR
                            snr = P.SNR_dB_vec(si);
                            [rxSig, symbols, PilotsAll] = ...
                                adaptive_eq_grid_sweep.buildChannel(P, net, snr, tr);

                            for mi = 1:NMODE
                                for ti = 1:NTAP
                                    for so = 1:NSIGN
                                        % Per-(tap, sign_only) step size and
                                        % re-init point from Mu_vec/N1_vec.
                                        Pr     = P;
                                        Pr.Mu  = mus(ti, so);
                                        Pr.N1  = n1s(ti, so);
                                        eqSym = adaptive_eq_grid_sweep.runAEQ( ...
                                            modes{mi}, rxSig, PilotsAll, taps(ti), ...
                                            signs(so), Pr, T);
                                        berCells{mi, ti, so}(tr, si) = ...
                                            adaptive_eq_grid_sweep.computeBER( ...
                                                eqSym, symbols, P.NOut);
                                    end
                                end
                            end
                        end
                    end

                    % Emit table rows for this (FL, net)
                    for mi = 1:NMODE
                        for ti = 1:NTAP
                            for so = 1:NSIGN
                                berMat = berCells{mi, ti, so};

                                row = row + 1;
                                mode_c(row)    = string(modes{mi});
                                network(row)   = net.name;
                                L_km(row)      = net.L;
                                splitting(row) = net.split;
                                ntaps_c(row)   = taps(ti);
                                fl(row)        = flv;
                                wl(row)        = P.IntBits + flv;
                                sign_only(row) = signs(so);
                                ber{row}       = berMat;
                            end
                        end
                    end
                end
            end

            tbl = table(mode_c, network, L_km, splitting, ntaps_c, fl, wl, ...
                        sign_only, ber, ...
                        'VariableNames', {'mode','network','L_km','splitting', ...
                                          'ntaps','fl','wl','sign_only','ber'});

            SNR_dB_vec = P.SNR_dB_vec;  %#ok<NASGU,PROP>
            params     = P;             %#ok<NASGU>

            outFile = fullfile(fileparts(mfilename('fullpath')), ...
                               'adaptive_eq_grid_sweep.mat');
            save(outFile, 'tbl', 'SNR_dB_vec', 'params');
            fprintf('Saved adaptive EQ grid sweep results to %s\n', outFile);
            disp(tbl(:, {'mode','network','ntaps','fl','wl'}));

            testCase.verifyEqual(height(tbl), nRows);
        end

    end

    %% ================================================================
    %  Helpers
    %% ================================================================
    methods (Static)

        function P = extractParams(tc)
            P.Rs          = tc.Rs;
            P.CWL         = tc.CWL;
            P.N_pol       = tc.N_pol;
            P.SpS         = tc.SpS;
            P.Rolloff     = tc.Rolloff;
            P.Span        = tc.Span;
            P.DGDSpec     = tc.DGDSpec;
            P.N_pmd       = tc.N_pmd;
            P.NTrials     = tc.NTrials;
            P.NSub        = tc.NSub;
            P.SNR_dB_vec  = tc.SNR_dB_vec;
            P.IntBits     = tc.IntBits;
            P.FL_vec      = tc.FL_vec;
            P.NTaps_vec    = tc.NTaps_vec;
            P.Modes        = tc.Modes;
            P.SignOnly_vec = tc.SignOnly_vec;
            P.Mu_vec      = tc.Mu_vec;
            P.N1_vec      = tc.N1_vec;
            % Scalar Mu / N1 are runtime-only build prototypes for codegen;
            % the test loop overrides them per (tap, sign_only) from Mu_vec/N1_vec.
            P.Mu          = tc.Mu_vec(1, 1);
            P.N1          = tc.N1_vec(1, 1);
            P.SingleSpike = tc.SingleSpike;
            P.NOut        = tc.NOut;
            P.UpdateStep  = tc.UpdateStep;
        end

        function nets = networks()
            % Two CPON network configurations (report \cref{chapter:energy}).
            nets = struct( ...
                'name',  {"A", "B"}, ...
                'L',     {20, 80}, ...
                'split', {"1:512", "1:16"});
        end

        function buildAEQMex(P, fxp)
            % Build the adaptive-EQ MEX for the given fixed-point config.
            % NTaps, Mode, Pilots, PLanes and BlockLen are runtime arguments,
            % so one binary serves the whole (NTaps x mode) sub-grid.
            clear mex %#ok<CLMEX>

            B.FxpConfig_AEQ   = fxp;
            B.SpS             = P.SpS;
            B.AEQ_NTaps       = max(P.NTaps_vec);   % prototype (runtime-variable)
            B.AEQ_Mu          = P.Mu;
            B.AEQ_SingleSpike = P.SingleSpike;
            B.AEQ_N1          = P.N1;
            B.AEQ_NOut        = P.NOut;
            B.AEQ_SignOnly    = false;
            B.AEQ_UpdateStep  = P.UpdateStep;
            B.AEQ_PLanes      = 32;     % prototype value (runtime-variable)
            B.AEQ_Mode        = 0;
            B.AEQ_BlockLen    = 32;     % prototype value (runtime-variable)

            cfg = coder.config('mex');
            cfg.GenerateReport = false;
            build_adaptive_eq_equalize_fxp_mex(B, cfg);
        end

        function [rxSig, symbols, PilotsAll] = buildChannel(P, net, SNR_dB, trialSeed)
            % CPON-framed DP-QPSK through AWGN + PMD for one network/SNR.
            % RRC pulse-shaping is applied at the transmitter; the matched
            % filter is omitted so the adaptive equaliser sees the raw
            % oversampled signal.
            rng(1000 * trialSeed + round(SNR_dB) + 7 * double(net.L));

            DATA_PER_SUBFRAME = 3586;
            nBits = P.NSub * DATA_PER_SUBFRAME * P.N_pol * 2;
            bits  = randi([0 1], nBits, 1);

            [symbols, pilots, ~, nSub] = modem.modulate(bits);
            PilotsAll = repmat(pilots, nSub, 1);     % [nSub*116 x 2]

            txSig = modem.rrcPulse(symbols, P.SpS, P.Rolloff, P.Span);
            rxSig = channel.add_awgn(txSig, SNR_dB);
            rxSig = channel.add_pmd(rxSig, net.L, P.SpS, P.Rs, ...
                        P.DGDSpec, P.N_pmd);
        end

        function eqSym = runAEQ(mode, rxSig, PilotsAll, ntaps, signOnly, P, T)
            % Per-mode PLanes/BlockLen following the validated shapes in
            % test_AdaptiveEqualizer:
            %   CMA   -> parallel-lane CMA: PLanes=32, BlockLen=PLanes=32.
            %            Mode=0, Pilots=[] (the pilot position of each
            %            32-symbol block is still skipped intrinsically).
            %   pilot -> data-aided LMS: PLanes=1, BlockLen=32, Mode=1,
            %            Pilots=PilotsAll.
            % Running CMA with PLanes=1/BlockLen=32 (the pilot shape)
            % causes the gradient to accumulate over 32 samples and be
            % applied as a single ~32*mu step, which is unstable at
            % Mu=1e-3 and diverges to NaN.
            % signOnly selects between the direct (false) and sign-sign
            % (true) variants of the weight update.
            rx_fi = cast(rxSig, 'like', T.x);
            switch mode
                case 'CMA'
                    pilots_fi = cast(complex(zeros(0, 2)), 'like', T.y);
                    Mode      = 0;
                    pLanes    = 32;
                    blockLen  = 32;
                    subBlocks = 116;   % CPON subframe = 116 blocks; skip
                                       % block 1 of each subframe so the
                                       % +/-3+/-3j training symbols don't
                                       % corrupt the CMA gradient.
                case 'pilot'
                    pilots_fi = cast(PilotsAll, 'like', T.y);
                    Mode      = 1;
                    pLanes    = 1;
                    blockLen  = 32;
                    subBlocks = 0;     % pilot mode unaffected
                otherwise
                    error('adaptive_eq_grid_sweep:badMode', ...
                          'Unknown mode: %s', mode);
            end

            eq = adaptive_eq.equalize_fxp_mex(rx_fi, double(P.SpS), ...
                double(ntaps), double(P.Mu), logical(P.SingleSpike), ...
                double(P.N1), double(P.NOut), logical(signOnly), ...
                double(P.UpdateStep), T, double(pLanes), ...
                double(Mode), pilots_fi, double(blockLen), double(subBlocks));
            eqSym = double(eq);
        end

        function BER = computeBER(eqSym, symbols, NOut)
            % BER against TX symbols (offset by NOut), resolving the
            % arbitrary residual phase by a dense rotation search (32
            % angles of pi/16, matching test_AdaptiveEqualizer) plus a
            % polarisation-swap loop.  The earlier 4th-power approach was
            % buggy for diagonal QPSK (it over-estimates the rotation by
            % pi/4, rotating the constellation onto the axes and making
            % decideSymbols tie-break between adjacent grid points -> 50%
            % BER per pol even on a converged equaliser).
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
                for q = 1:size(eqSym, 2)             % polarisation swap
                    for kk = 0:31                    % rotations of pi/16
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
