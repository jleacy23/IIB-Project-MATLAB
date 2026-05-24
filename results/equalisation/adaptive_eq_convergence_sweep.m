classdef adaptive_eq_convergence_sweep < matlab.unittest.TestCase
%ADAPTIVE_EQ_CONVERGENCE_SWEEP  Convergence diagnostic sweep for the AEQ.
%
%   Holds the gradient precision fixed at FL = 16 (high precision) and
%   sweeps the three knobs that govern adaptive-equaliser convergence:
%
%       1. Tap length      NTaps_vec  (1, 3, 5)
%       2. Step size       Mu_vec     (CMA learning rate)
%       3. Re-init symbol  N1_vec     (y-polarisation single-spike re-init
%                                      iteration; "burn-in" length before
%                                      the y-pol is forced to a spike)
%
%   The intent is to find which (Mu, N1) combinations let the longer-tap
%   equalisers (N = 3, 5) converge for the grid sweep — at the default
%   Mu = 1e-3, N1 = 2000 used in adaptive_eq_grid_sweep, only N = 1 taps
%   reach the FEC threshold.  Direct-update CMA only (the case the grid
%   sweep observed diverging) and a single high SNR are used so that the
%   BER signal is dominated by convergence behaviour rather than noise.
%
%   Output (saved next to this file):
%       adaptive_eq_convergence_sweep.mat  with variables:
%         tbl     - results table, one row per
%                     (network, NTaps, Mu, N1):
%                     network   (string)   'A' | 'B'
%                     L_km      (double)
%                     ntaps     (double)
%                     mu        (double)
%                     n1        (double)
%                     ber       (double)   minimum BER across trials
%                     ber_trials(double[]) per-trial BER
%         params  - copy of the run parameters
%
%   Run with:
%       runtests('adaptive_eq_convergence_sweep')

    %% ================================================================
    %  Parameters
    %% ================================================================
    properties (Constant)

        % System (mirrors adaptive_eq_grid_sweep)
        Rs        = 30.5
        CWL       = 1550
        N_pol     = 2
        SpS       = 2
        DGDSpec   = 0.1
        N_pmd     = 1

        % Pulse shaping
        Rolloff   = 0.25
        Span      = 10

        % Monte-Carlo
        NTrials   = 3
        NSub      = 8

        % Single high SNR — convergence, not noise floor, is being probed.
        SNR_dB    = 20

        % Fixed-point: FL = 16, WL = IntBits + FL.
        IntBits   = 16
        FL        = 16

        % Convergence sweep axes
        NTaps_vec = [1, 3, 5]
        Mu_vec    = [2^(-10), 2^(-12), 2^(-14)]
        N1_vec    = [1000, 2000, 4000]

        % Held fixed (direct CMA only — the diverging case).
        Mode_str    = 'CMA'
        SignOnly    = true
        SingleSpike = true
        NOut        = 5000
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

        function test_convergence_sweep(testCase)
            P = adaptive_eq_convergence_sweep.extractParams(testCase);

            % Reuse the network list and helpers from the grid sweep so
            % the channel, RNG seeding and BER calculation match exactly.
            nets  = adaptive_eq_grid_sweep.networks();
            taps  = testCase.NTaps_vec;
            mus   = testCase.Mu_vec;
            n1s   = testCase.N1_vec;

            NNET = numel(nets);
            NTAP = numel(taps);
            NMU  = numel(mus);
            NN1  = numel(n1s);
            nRows = NNET * NTAP * NMU * NN1;

            network   = strings(nRows, 1);
            L_km      = nan(nRows, 1);
            ntaps_c   = nan(nRows, 1);
            mu        = nan(nRows, 1);
            n1        = nan(nRows, 1);
            ber       = nan(nRows, 1);
            ber_trials = cell(nRows, 1);

            % Mu and N1 are runtime arguments to adaptive_eq.equalize_fxp_mex
            % (build script passes them as plain doubles, not coder.Constant),
            % so the MEX only needs to be built once at FL = 16.  NTaps is
            % also a runtime arg.  Channel realisations are shared across the
            % full (NTaps, Mu, N1) sweep so all rows see identical noise.
            fxp = struct('WL', P.IntBits + P.FL, 'FL', P.FL);
            fprintf('=== Building AEQ MEX  FL = %d  (WL = %d) ===\n', ...
                    P.FL, P.IntBits + P.FL);
            adaptive_eq_grid_sweep.buildAEQMex(P, fxp);
            T = adaptive_eq.equalize_fxp_types(fxp);

            row = 0;
            for ni = 1:NNET
                net = nets(ni);

                rxBag    = cell(P.NTrials, 1);
                symBag   = cell(P.NTrials, 1);
                pilotBag = cell(P.NTrials, 1);
                for tr = 1:P.NTrials
                    [rxBag{tr}, symBag{tr}, pilotBag{tr}] = ...
                        adaptive_eq_grid_sweep.buildChannel( ...
                            P, net, P.SNR_dB, tr);
                end

                for ti = 1:NTAP
                    nt = taps(ti);
                    for mi = 1:NMU
                        for n1i = 1:NN1
                            muVal = mus(mi);
                            n1Val = n1s(n1i);

                            % Patch only the two swept fields; runAEQ reads
                            % these directly when calling the MEX.
                            Pr    = P;
                            Pr.Mu = muVal;
                            Pr.N1 = n1Val;

                            berT = nan(P.NTrials, 1);
                            for tr = 1:P.NTrials
                                eqSym = adaptive_eq_grid_sweep.runAEQ( ...
                                    P.Mode_str, rxBag{tr}, pilotBag{tr}, ...
                                    nt, P.SignOnly, Pr, T);
                                berT(tr) = adaptive_eq_grid_sweep.computeBER( ...
                                    eqSym, symBag{tr}, P.NOut);
                            end

                            row = row + 1;
                            network(row)    = net.name;
                            L_km(row)       = net.L;
                            ntaps_c(row)    = nt;
                            mu(row)         = muVal;
                            n1(row)         = n1Val;
                            ber(row)        = min(berT, [], 'omitnan');
                            ber_trials{row} = berT(:).';

                            fprintf(['    net=%s  N=%d  Mu=%.1e  N1=%d  ', ...
                                     '-> minBER = %.3g\n'], ...
                                    net.name, nt, muVal, n1Val, ber(row));
                        end
                    end
                end
            end

            tbl = table(network, L_km, ntaps_c, mu, n1, ber, ber_trials, ...
                        'VariableNames', {'network','L_km','ntaps','mu','n1', ...
                                          'ber','ber_trials'});
            params = P;

            outFile = fullfile(fileparts(mfilename('fullpath')), ...
                               'adaptive_eq_convergence_sweep.mat');
            save(outFile, 'tbl', 'params');
            fprintf('Saved AEQ convergence sweep to %s\n', outFile);
            disp(tbl(:, {'network','ntaps','mu','n1','ber'}));

            % Best (Mu, N1) per (network, NTaps): minimum mean BER, ties
            % broken by the first occurrence.  NaN BERs are ignored.
            fprintf('\n=== Best (Mu, N1) per tap length ===\n');
            netsAll = unique(tbl.network, 'stable');
            for ni = 1:numel(netsAll)
                netName = netsAll(ni);
                for ti = 1:NTAP
                    nt   = taps(ti);
                    rows = tbl(tbl.network == netName & tbl.ntaps == nt, :);
                    if isempty(rows), continue; end
                    [bestBER, idx] = min(rows.ber, [], 'omitnan');
                    if isnan(bestBER)
                        fprintf('  net=%s  N=%d  no converged trials\n', ...
                                netName, nt);
                    else
                        fprintf(['  net=%s  N=%d  Mu=%.3e  N1=%d  ', ...
                                 '-> minBER = %.3g\n'], ...
                                netName, nt, rows.mu(idx), rows.n1(idx), bestBER);
                    end
                end
            end

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
            P.SNR_dB      = tc.SNR_dB;
            P.IntBits     = tc.IntBits;
            P.FL          = tc.FL;
            P.NTaps_vec   = tc.NTaps_vec;
            P.Mu_vec      = tc.Mu_vec;
            P.N1_vec      = tc.N1_vec;
            P.Mode_str    = tc.Mode_str;
            P.SignOnly    = tc.SignOnly;
            P.SingleSpike = tc.SingleSpike;
            P.NOut        = tc.NOut;
            P.UpdateStep  = tc.UpdateStep;

            % buildAEQMex consumes these scalar fields (set per outer-loop
            % iteration before each MEX build).
            P.Mu    = NaN;
            P.N1    = NaN;
            % buildAEQMex prototypes the NTaps dimension at max(P.NTaps_vec).
        end

    end
end
