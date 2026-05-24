classdef cd_eq_precision_sweep < matlab.unittest.TestCase
%CD_EQ_PRECISION_SWEEP  Fixed-point precision sweep of the CD equalisers.
%
%   Sweeps the fixed-point fractional bit width over the three chromatic
%   dispersion (CD) equaliser implementations of report/full/full.tex
%   (tab:cd_cost), for both network configurations (tab:cd_taps), and logs
%   the raw per-iteration BER.  No FEC-SNR threshold is computed —
%   downstream post-processing decides what to derive from the BER cube.
%
%   The three CD equaliser configurations are:
%       'time_domain'       — direct FIR convolution (cd_eq.equalize_td_fxp).
%                             The filter length is supplied by the caller
%                             via cd_eq.computeOverlap (shared with the
%                             overlap-save sizing); the report's N_CD
%                             (tab:cd_taps) is also logged as metadata
%                             (ncd_report) for reference.
%       'overlap_save'      — overlap-save frequency-domain (cd_eq.equalize_fxp,
%                             po2Twiddle = false).
%       'overlap_save_po2'  — overlap-save with power-of-two twiddle factors
%                             (cd_eq.equalize_fxp, po2Twiddle = true).
%
%   The two networks (report \cref{chapter:energy} / tab:cd_taps,
%   tab:cd_cost_eval), with the FFT length N and tap count N_CD taken from
%   the report:
%       Net A : L = 20 km, 1:512 split, N_CD = 6,  N_FFT = 32
%       Net B : L = 80 km, 1:16  split, N_CD = 21, N_FFT = 128
%
%   Channel chain (CD equaliser characterised in isolation):
%       RRC pulse-shape (Tx) -> CD -> AWGN -> CD equaliser (fxp under
%       test) -> RRC matched filter (Rx) -> downsample -> BER.
%   The RRC + RRC matched-filter Nyquist split mirrors the canonical
%   coherent-receiver chain (see results/full_pipeline/run_pipeline.m).
%
%   Each (config, network, FL) point is run over NTrials independent
%   Monte-Carlo realisations and the full SNR sweep.
%
%   Output (saved next to this file):
%       cd_eq_precision_sweep.mat  with variables:
%         tbl     — results table, one row per (config, network, FL):
%                     config     (string)
%                     network    (string)  'A' | 'B'
%                     L_km       (double)
%                     splitting  (string)
%                     n_fft      (double)  report FFT length
%                     ncd_report (double)  report tap count
%                     fl         (double)  fractional bits
%                     wl         (double)  word length
%                     ber        (cell)    [NTrials x NSNR] raw BER
%         SNR_dB_vec — SNR sweep axis [dB]
%         params     — copy of the run parameters
%
%   Run with:
%       runtests('cd_eq_precision_sweep')

    %% ================================================================
    %  Parameters
    %% ================================================================
    properties (Constant)

        % System (report tab:network_params)
        Rs        = 30.5         % symbol / channel bandwidth [GBd]
        D         = 17           % dispersion [ps/(nm*km)]
        CWL       = 1550         % central wavelength [nm]
        N_pol     = 2
        SpS       = 2            % samples per symbol

        % Pulse shaping (Nyquist / raised-cosine)
        Rolloff   = 0.25
        Span      = 10

        % Monte-Carlo
        NTrials   = 10
        Ns        = 8192         % symbols per polarisation per trial

        % SNR sweep [dB]
        SNR_dB_vec = 0 : 2 : 30

        % Bit-width sweep — integer bits fixed; WL = IntBits + FL
        IntBits   = 16
        FL_vec    = [2, 4, 6, 8, 10, 12]

        % CD equaliser configurations
        Configs   = {'time_domain', 'overlap_save', 'overlap_save_po2'}

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

        function test_precision_sweep(testCase)
            P    = cd_eq_precision_sweep.extractParams(testCase);
            nets = cd_eq_precision_sweep.networks();
            cfgs = testCase.Configs;
            FLs  = testCase.FL_vec;
            NSNR = numel(P.SNR_dB_vec);

            NFL  = numel(FLs);
            NNET = numel(nets);
            NCFG = numel(cfgs);
            nRows = NFL * NNET * NCFG;

            % Pre-allocate table columns
            config     = strings(nRows, 1);
            network    = strings(nRows, 1);
            L_km       = nan(nRows, 1);
            splitting  = strings(nRows, 1);
            n_fft      = nan(nRows, 1);
            ncd_report = nan(nRows, 1);
            fl         = nan(nRows, 1);
            wl         = nan(nRows, 1);
            ber        = cell(nRows, 1);

            row = 0;
            for bi = 1:NFL
                flv = FLs(bi);
                fxp = struct('WL', P.IntBits + flv, 'FL', flv);

                fprintf('=== FL = %2d  (WL = %2d)  [%d/%d] : building CD MEX ===\n', ...
                        flv, P.IntBits + flv, bi, NFL);
                cd_eq_precision_sweep.buildCDMex(P, fxp);

                T_fd = cd_eq.equalize_fxp_types(fxp);
                T_td = cd_eq.equalize_td_fxp_types(fxp);

                for ni = 1:NNET
                    net = nets(ni);

                    % Generate each channel realisation once and run all
                    % three configs on it so they see identical inputs.
                    berCells = cell(NCFG, 1);
                    for ci = 1:NCFG
                        berCells{ci} = nan(P.NTrials, NSNR);
                    end

                    for tr = 1:P.NTrials
                        fprintf('    FL=%2d net=%s  trial %d/%d\n', ...
                                flv, net.name, tr, P.NTrials);
                        for si = 1:NSNR
                            snr = P.SNR_dB_vec(si);
                            [rxSig, refSym] = cd_eq_precision_sweep.buildChannel( ...
                                P, net, snr, tr);

                            for ci = 1:NCFG
                                eqSym = cd_eq_precision_sweep.runCDEq( ...
                                    cfgs{ci}, rxSig, P, net, T_fd, T_td);
                                berCells{ci}(tr, si) = ...
                                    cd_eq_precision_sweep.computeBER(eqSym, refSym);
                            end
                        end
                    end

                    for ci = 1:NCFG
                        berMat = berCells{ci};

                        row = row + 1;
                        config(row)     = string(cfgs{ci});
                        network(row)    = net.name;
                        L_km(row)       = net.L;
                        splitting(row)  = net.split;
                        n_fft(row)      = net.NFFT;
                        ncd_report(row) = net.Ncd;
                        fl(row)         = flv;
                        wl(row)         = P.IntBits + flv;
                        ber{row}        = berMat;
                    end
                end
            end

            tbl = table(config, network, L_km, splitting, n_fft, ncd_report, ...
                        fl, wl, ber);

            SNR_dB_vec = P.SNR_dB_vec;  %#ok<NASGU,PROP>
            params     = P;             %#ok<NASGU>

            outFile = fullfile(fileparts(mfilename('fullpath')), ...
                               'cd_eq_precision_sweep.mat');
            save(outFile, 'tbl', 'SNR_dB_vec', 'params');
            fprintf('Saved CD precision sweep results to %s\n', outFile);
            disp(tbl(:, {'config','network','fl','wl'}));

            testCase.verifyEqual(height(tbl), nRows);
        end

    end

    %% ================================================================
    %  Helpers
    %% ================================================================
    methods (Static)

        function P = extractParams(tc)
            P.Rs         = tc.Rs;
            P.D          = tc.D;
            P.CWL        = tc.CWL;
            P.N_pol      = tc.N_pol;
            P.SpS        = tc.SpS;
            P.Rolloff    = tc.Rolloff;
            P.Span       = tc.Span;
            P.NTrials    = tc.NTrials;
            P.Ns         = tc.Ns;
            P.SNR_dB_vec = tc.SNR_dB_vec;
            P.IntBits    = tc.IntBits;
            P.FL_vec     = tc.FL_vec;
            P.Configs    = tc.Configs;
        end

        function nets = networks()
            % Two CPON network configurations from report tab:cd_taps /
            % tab:cd_cost_eval.  N_CD and N_FFT are the report values.
            nets = struct( ...
                'name',  {"A", "B"}, ...
                'L',     {20, 80}, ...
                'split', {"1:512", "1:16"}, ...
                'Ncd',   {6, 21}, ...
                'NFFT',  {32, 128});
        end

        function buildCDMex(P, fxp)
            % Build both CD MEX binaries for the given fixed-point config.
            % NFFT, L and po2Twiddle are runtime arguments, so a single pair
            % of binaries serves both networks and both overlap-save modes.
            clear mex %#ok<CLMEX>

            B.FxpConfig_CD = fxp;
            B.N_pol        = P.N_pol;
            B.D            = P.D;
            B.L            = 80;        % prototype value (runtime-variable)
            B.CWL          = P.CWL;
            B.Rs           = P.Rs;
            B.SpS          = P.SpS;
            B.NFFT         = 128;       % prototype value (runtime-variable)
            B.po2Twiddle   = false;     % prototype value (runtime-variable)

            cfg = coder.config('mex');
            cfg.GenerateReport = false;
            build_cd_eq_equalize_fxp_mex(B, cfg);
            build_cd_eq_equalize_td_fxp_mex(B, cfg);
        end

        function [rxSig, refSym] = buildChannel(P, net, SNR_dB, trialSeed)
            % CD + AWGN channel for one network at one SNR.  Tx uses RRC
            % pulse-shaping; the matching RRC matched filter is applied
            % AFTER the CD equaliser inside runCDEq (so Tx*Rx combines to
            % a Nyquist raised cosine and the equaliser is characterised
            % in a canonical coherent-receiver chain).
            rng(1000 * trialSeed + round(SNR_dB) + 7 * double(net.L));

            txBits  = modem.randomBits(4 * P.Ns);
            symbols = modem.modulate(txBits);
            txSig   = modem.rrcPulse(symbols, P.SpS, P.Rolloff, P.Span);

            rxSig = channel.add_chromatic_dispersion(txSig, net.L, P.SpS, ...
                        P.Rs, P.D, P.CWL);
            rxSig = channel.add_awgn(rxSig, SNR_dB);

            refSym = symbols;
        end

        function eqSym = runCDEq(cfg, rxSig, P, net, T_fd, T_td)
            % Run the selected CD equaliser (MEX), apply the Rx matched
            % filter (RRC, matched to rrcPulse at Tx), and return the
            % symbol-rate output (downsampled by SpS).
            switch cfg
                case 'time_domain'
                    rx_fi = cast(rxSig, 'like', T_td.x);
                    NTap  = cd_eq.computeOverlap(double(P.D), double(net.L), ...
                        double(P.CWL), double(P.Rs), double(P.SpS), ...
                        double(net.NFFT));
                    eqSig = cd_eq.equalize_td_fxp_mex(rx_fi, double(P.D), ...
                        double(net.L), double(P.CWL), double(P.Rs), ...
                        double(P.N_pol), double(P.SpS), double(NTap), T_td);
                case 'overlap_save'
                    rx_fi = cast(rxSig, 'like', T_fd.x);
                    eqSig = cd_eq.equalize_fxp_mex(rx_fi, double(P.D), ...
                        double(net.L), double(P.CWL), double(P.Rs), ...
                        double(P.N_pol), double(P.SpS), double(net.NFFT), ...
                        false, T_fd);
                case 'overlap_save_po2'
                    rx_fi = cast(rxSig, 'like', T_fd.x);
                    eqSig = cd_eq.equalize_fxp_mex(rx_fi, double(P.D), ...
                        double(net.L), double(P.CWL), double(P.Rs), ...
                        double(P.N_pol), double(P.SpS), double(net.NFFT), ...
                        true, T_fd);
                otherwise
                    error('cd_eq_precision_sweep:badConfig', ...
                          'Unknown CD config: %s', cfg);
            end
            % --- Rx matched filter (after the CD equaliser) ---
            mfOut = modem.matched_filter(double(eqSig), P.SpS, 'rrc', ...
                        P.Rolloff, P.Span);
            eqSym = mfOut(1:P.SpS:end, :);
        end

        function BER = computeBER(eqSym, refSym)
            % BER with QPSK pi/2-rotation ambiguity resolved per polarisation.
            n      = min(size(eqSym, 1), size(refSym, 1));
            eqSym  = eqSym(1:n, :);
            refSym = refSym(1:n, :);
            rotations = [1, 1j, -1, -1j];

            totErr = 0; totBits = 0;
            for p = 1:size(refSym, 2)
                refBits = modem.symbolsToBits(refSym(:, p));
                best    = Inf;
                for ri = 1:4
                    dec  = modem.decideSymbols(eqSym(:, p) * rotations(ri));
                    bits = modem.symbolsToBits(dec);
                    m    = min(numel(refBits), numel(bits));
                    e    = sum(refBits(1:m) ~= bits(1:m)) / m;
                    best = min(best, e);
                end
                totErr  = totErr  + best * numel(refBits);
                totBits = totBits + numel(refBits);
            end
            BER = totErr / totBits;
        end

    end
end
