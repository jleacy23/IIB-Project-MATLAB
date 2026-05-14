classdef bit_width_full < matlab.unittest.TestCase
%BIT_WIDTH_FULL  Full bit-width sweep with FEC SNR and energy per bit.
%
%   For every combination of frequency recovery
%       {fft_search (data-aided), differential_kay (data-aided),
%        fft_search_blind}
%   and phase recovery
%       {viterbi-viterbi, pilots-only}
%   the test sweeps the fractional bit width of the frequency-recovery
%   subsystem and the carrier-recovery subsystem separately
%   (FL = 2..16, integer bits = 16) and records:
%
%     * FEC SNR — SNR at which the post-CR BER crosses FEC_BER (linear
%                 interpolation between samples, see fecCrossing).
%     * Energy per bit — operation-count model from
%                 report/frequency_recovery/frequency_recovery.tex
%                 (tab:fft_cost / tab:diffkay_cost / tab:viterbi_cost /
%                 tab:pilot_cost) fed to src/+energy/receiver.m.
%
%   The blind FFT search is additionally swept over observation length
%   BlindD_vec, since this is the dominant accuracy/energy lever for that
%   variant.  All FR MEX binaries are built once with FR_BlindD = max so
%   the runtime D argument can vary without rebuilding.
%
%   test_full_grid_sweep — full Cartesian sweep over (fl_fr, fl_cr) for
%       every (FR, CR) pair, plus BlindD for the blind FFT variant.
%       Output: bit_width_full_grid_sweep.mat
%
%   The output table has columns:
%       fr_algo   (string)  — 'fft_search' | 'differential_kay' | 'fft_search_blind'
%       cr_algo   (string)  — 'viterbi_viterbi' | 'pilots_only'
%       fl_fr     (double)  — FR fractional bits
%       fl_cr     (double)  — CR fractional bits
%       blind_d   (double)  — blind observation length (NaN for non-blind)
%       fec_snr_trials (cell) — per-trial FEC SNR thresholds [dB], NTrials×1 (NaN if curve never crosses)
%       energy_fr_fJ (double) — FR energy per bit [fJ]
%       energy_cr_fJ (double) — CR energy per bit [fJ]
%
%   Run with:
%       runtests('bit_width_full')
%       runtests('bit_width_full', 'ProcedureName', 'test_full_grid_sweep')

    %% ================================================================
    %  Parameters
    %% ================================================================
    properties (Constant)

        % System
        Rs          = 30.5
        N_pol       = 2
        TrainingLen = 11

        % Monte-Carlo
        NTrials     = 100

        % SNR sweep
        SNR_dB_vec  = 0 : 1 : 30

        % Bit width sweep — integer bits fixed at 16; WL = IntBits + FL
        IntBits     = 16
        FL_vec      = [2, 4, 6, 8, 10, 12, 14, 16]
        FL_fixed    = 16

        % Channel conditions
        DeltaF_Hz   = 2e9
        LW_Hz       = 1000e3

        % FFT search parameters
        FR_Nfft       = 512
        FR_Po2Twiddle = false
        MaxFreq       = 0.1

        % Blind FFT observation length — BlindD is the max baked into
        % the MEX type (FR_BlindD); BlindD_vec holds the runtime sweep.
        BlindD     = 512
        BlindD_vec = [32, 64, 128, 256, 512]

        % Phase recovery
        CordicIts      = 16
        BlockLen       = 32
        StepSize       = 32
        PilotThreshold = 5 * pi / 9
        VV_NTaps       = 10

        % FEC threshold
        FEC_BER = 2e-2

        % Energy model — fits to horowitz2014computing
        % E_A(n) = EAdd_fJ * n,  E_M(n) = EMult_fJ * n^2
        EAdd_fJ      = 3.16
        EMult_fJ     = 3.03
        SubframeLen  = 3712
        M            = 4
        Oversampling = 1

    end

    properties
        VVFilters   % {1 x NSNR} Wiener VV filter taps, one per SNR point
    end

    %% ================================================================
    %  Test class setup
    %% ================================================================
    methods (TestClassSetup)

        function setupPath(~)
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src'));
        end

        function seedRng(~)
            rng(42);
        end

        function setupBuildPath(~)
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'build'));
        end

        function calibrateVVFilters(testCase)
            BITS_PER_SF = 3586 * 2 * 2;
            [tmp, ~, ~, ~] = modem.modulate(modem.randomBits(BITS_PER_SF));
            symEnergy = mean(abs(tmp(:)).^2);

            NSNR = length(testCase.SNR_dB_vec);
            testCase.VVFilters = cell(1, NSNR);
            for si = 1:NSNR
                testCase.VVFilters{si} = carrier_recovery.genVVFilter( ...
                    testCase.LW_Hz, testCase.Rs, testCase.SNR_dB_vec(si), ...
                    symEnergy, testCase.N_pol, testCase.VV_NTaps);
            end
        end

    end

    %% ================================================================
    %  Tests
    %% ================================================================
    methods (Test)

        function test_full_grid_sweep(testCase)
            % Full Cartesian grid over (fl_fr, fl_cr) for every FR x CR pair.
            % BlindD is additionally swept for the fft_search_blind variant.
            P        = testCase;
            NFL      = length(P.FL_vec);
            NBD      = length(P.BlindD_vec);
            NTrials  = P.NTrials;
            Prms     = bit_width_full.extractParams(testCase);

            % ---- Phase 1: serial MEX builds — FR variants live in temp dirs
            % (added to the worker's path inside parfor); CR variants are
            % swapped into src/+carrier_recovery serially in the outer loop.
            frDirs = cell(1, NFL);
            for bi = 1:NFL
                fl = P.FL_vec(bi);
                fprintf('[FR FL = %2d] building FR MEX  (%d / %d)\n', fl, bi, NFL);
                frDirs{bi} = bit_width_full.buildFRMex(P, struct('WL', P.IntBits + fl, 'FL', fl));
            end

            % Try to ensure the parpool has at least NFL workers so each
            % parfor iteration runs on a unique worker (no FR MEX cache
            % collisions).  Best effort — if Parallel Computing Toolbox
            % isn't available, parfor falls back to serial execution and
            % MEX caching isn't an issue.
            havePool = bit_width_full.ensureParpool(NFL);

            % ---- Phase 2: outer serial over CR FL; inner parfor over FR FL.
            % Indexing: fec*( fr_idx, cr_idx, cr_col )            for non-blind
            %           fec_blind( fr_idx, cr_idx, bd_idx, cr_col ) for blind
            % cr_col 1 = Viterbi-Viterbi, 2 = Pilots-only.
            fecSNR_fft_DA    = nan(NFL, NFL, NTrials, 2);
            fecSNR_dk_DA     = nan(NFL, NFL, NTrials, 2);
            fecSNR_fft_blind = nan(NFL, NFL, NBD, NTrials, 2);

            FL_vec_b     = P.FL_vec;
            IntBits_b    = P.IntBits;
            SNR_dB_b     = P.SNR_dB_vec;
            FEC_BER_b    = P.FEC_BER;
            BlindD_vec_b = P.BlindD_vec;
            NTrials_b    = NTrials;

            for cri = 1:NFL
                fl_cr = P.FL_vec(cri);
                fprintf('=== CR FL = %2d  (%d / %d) ===\n', fl_cr, cri, NFL);

                fxp_cr = struct('WL', P.IntBits + fl_cr, 'FL', fl_cr);
                T_cr_w = carrier_recovery.fxp_types(fxp_cr);

                % Build CR MEX in src/+carrier_recovery (overwrites). The
                % build script + buildCRMexInSrc both call `clear mex` on the
                % main process, but workers retain their own caches, so
                % invalidate them too before the next parfor.
                bit_width_full.buildCRMexInSrc(P, fxp_cr);
                if havePool
                    wait(parfevalOnAll(gcp, @bit_width_full.clearWorkerMex, 0));
                end

                slice_fft_DA = nan(NFL, NTrials_b, 2);
                slice_dk_DA  = nan(NFL, NTrials_b, 2);
                slice_blind  = nan(NFL, NBD, NTrials_b, 2);

                parfor fri = 1:NFL
                    addpath(frDirs{fri});  %#ok<PFBNS>

                    fl_fr  = FL_vec_b(fri);
                    T_fr_w = freq_recovery.fxp_types(struct('WL', IntBits_b + fl_fr, 'FL', fl_fr));

                    ber = bit_width_full.runSnrSweepStatic(Prms, 'fft_search', 0, T_fr_w, T_cr_w);
                    v_vv = bit_width_full.fecCrossingsAll(SNR_dB_b, ber(:,:,1), FEC_BER_b);
                    v_po = bit_width_full.fecCrossingsAll(SNR_dB_b, ber(:,:,2), FEC_BER_b);
                    slice_fft_DA(fri, :, :) = reshape([v_vv(:), v_po(:)], 1, NTrials_b, 2);

                    ber = bit_width_full.runSnrSweepStatic(Prms, 'differential_kay', 0, T_fr_w, T_cr_w);
                    v_vv = bit_width_full.fecCrossingsAll(SNR_dB_b, ber(:,:,1), FEC_BER_b);
                    v_po = bit_width_full.fecCrossingsAll(SNR_dB_b, ber(:,:,2), FEC_BER_b);
                    slice_dk_DA(fri, :, :) = reshape([v_vv(:), v_po(:)], 1, NTrials_b, 2);

                    row_blind = nan(NBD, NTrials_b, 2);
                    for di = 1:length(BlindD_vec_b)
                        Dval = BlindD_vec_b(di);
                        ber  = bit_width_full.runSnrSweepStatic(Prms, 'fft_search_blind', Dval, T_fr_w, T_cr_w);
                        v_vv = bit_width_full.fecCrossingsAll(SNR_dB_b, ber(:,:,1), FEC_BER_b);
                        v_po = bit_width_full.fecCrossingsAll(SNR_dB_b, ber(:,:,2), FEC_BER_b);
                        row_blind(di, :, :) = reshape([v_vv(:), v_po(:)], 1, NTrials_b, 2);
                    end
                    slice_blind(fri, :, :, :) = row_blind;
                end

                fecSNR_fft_DA(:, cri, :, :)       = slice_fft_DA;
                fecSNR_dk_DA(:, cri, :, :)        = slice_dk_DA;
                fecSNR_fft_blind(:, cri, :, :, :) = slice_blind;

                % Workers may have cached FR MEX bound to the fri they ran.
                % Clear before the next outer iteration so a new fri-to-worker
                % assignment picks up the correct FR binary from its temp dir.
                if havePool
                    wait(parfevalOnAll(gcp, @bit_width_full.clearWorkerMex, 0));
                end
            end

            % ---- Phase 3: energy + table ----
            E_fft_per_fl     = bit_width_full.frEnergy(P, 'fft_search',       P.FL_vec, []);
            E_dk_per_fl      = bit_width_full.frEnergy(P, 'differential_kay', P.FL_vec, []);
            E_fft_blind_grid = bit_width_full.frEnergy(P, 'fft_search_blind', P.FL_vec, P.BlindD_vec);
            E_vv_per_fl      = bit_width_full.crEnergy(P, 'viterbi_viterbi',  P.FL_vec);
            E_po_per_fl      = bit_width_full.crEnergy(P, 'pilots_only',      P.FL_vec);

            tbl = bit_width_full.assembleGridTable(P, ...
                fecSNR_fft_DA, fecSNR_dk_DA, fecSNR_fft_blind, ...
                E_fft_per_fl, E_dk_per_fl, E_fft_blind_grid, ...
                [E_vv_per_fl, E_po_per_fl]);

            outDir  = fileparts(mfilename('fullpath'));
            outFile = fullfile(outDir, 'bit_width_full_grid_sweep.mat');
            save(outFile, 'tbl');
            fprintf('Saved full grid sweep results to %s\n', outFile);
            disp(tbl);

            cellfun(@(d) rmdir(d, 's'), frDirs, 'UniformOutput', false);
        end

    end

    %% ================================================================
    %  Static helpers
    %% ================================================================
    methods (Static, Access = private)

        % ---------------- Table assembly -----------------------------

        function T = assembleGridTable(P, fec_fft, fec_dk, fec_blind, ...
                E_fft_per_fl, E_dk_per_fl, E_fft_blind_grid, E_cr_grid)
            % fec_fft, fec_dk:  NFL x NFL x NTrials x 2  (fri, cri, trial, cr_col) — per-trial FEC SNR
            % fec_blind:        NFL x NFL x NBD x NTrials x 2
            % E_fft_per_fl:     NFL x 1  (FR fft energy per FR FL)
            % E_dk_per_fl:      NFL x 1
            % E_fft_blind_grid: NFL x NBD  (FR FL × BlindD)
            % E_cr_grid:        NFL x 2  (CR FL × {VV, PO})
            NFL    = length(P.FL_vec);
            NBD    = length(P.BlindD_vec);
            crNames = {"viterbi_viterbi", "pilots_only"};

            N_total        = 2 * (2 * NFL * NFL + NFL * NFL * NBD);
            fr_algo        = strings(N_total, 1);
            cr_algo        = strings(N_total, 1);
            fl_fr          = nan(N_total, 1);
            fl_cr          = nan(N_total, 1);
            blind_d        = nan(N_total, 1);
            fec_snr_trials = cell(N_total, 1);
            energy_fr_fJ   = nan(N_total, 1);
            energy_cr_fJ   = nan(N_total, 1);

            row = 0;
            for c = 1:2
                for fri = 1:NFL
                    for cri = 1:NFL
                        row = row + 1;
                        fr_algo(row)        = "fft_search";
                        cr_algo(row)        = crNames{c};
                        fl_fr(row)          = P.FL_vec(fri);
                        fl_cr(row)          = P.FL_vec(cri);
                        blind_d(row)           = NaN;
                        fec_snr_trials{row}    = squeeze(fec_fft(fri, cri, :, c));
                        energy_fr_fJ(row)      = E_fft_per_fl(fri);
                        energy_cr_fJ(row)      = E_cr_grid(cri, c);

                        row = row + 1;
                        fr_algo(row)           = "differential_kay";
                        cr_algo(row)           = crNames{c};
                        fl_fr(row)             = P.FL_vec(fri);
                        fl_cr(row)             = P.FL_vec(cri);
                        blind_d(row)           = NaN;
                        fec_snr_trials{row}    = squeeze(fec_dk(fri, cri, :, c));
                        energy_fr_fJ(row)      = E_dk_per_fl(fri);
                        energy_cr_fJ(row)      = E_cr_grid(cri, c);

                        for di = 1:NBD
                            row = row + 1;
                            fr_algo(row)        = "fft_search_blind";
                            cr_algo(row)        = crNames{c};
                            fl_fr(row)          = P.FL_vec(fri);
                            fl_cr(row)          = P.FL_vec(cri);
                            blind_d(row)        = P.BlindD_vec(di);
                            fec_snr_trials{row} = squeeze(fec_blind(fri, cri, di, :, c));
                            energy_fr_fJ(row)   = E_fft_blind_grid(fri, di);
                            energy_cr_fJ(row)   = E_cr_grid(cri, c);
                        end
                    end
                end
            end

            T = table(fr_algo, cr_algo, fl_fr, fl_cr, blind_d, ...
                fec_snr_trials, energy_fr_fJ, energy_cr_fJ);
        end

    end

    %% ================================================================
    %  Public helper — invoked on workers via parfevalOnAll
    %% ================================================================
    methods (Static)

        function clearWorkerMex()
            % `clear mex` is forbidden inside a parfor body but is fine
            % when invoked via parfevalOnAll, since the worker executes
            % it as a top-level function call between parfor iterations.
            clear mex %#ok<CLMEX>
        end

        function havePool = ensureParpool(N)
            % Best-effort: ensure a parpool of at least N workers exists.
            % Returns true if a pool is available, false if PCT isn't
            % installed or the pool can't be created.
            havePool = false;
            try
                pool = gcp('nocreate');
                if isempty(pool) || pool.NumWorkers < N
                    if ~isempty(pool), delete(pool); end
                    parpool('local', N);
                end
                havePool = true;
            catch ME
                fprintf('Parallel pool unavailable (%s); running serially.\n', ME.identifier);
            end
        end

    end

    methods (Static, Access = private)

        % ---------------- Energy --------------------------------------

        function E = frEnergy(P, algo, fl_vec, blindD_vec)
            % Returns FR energy per bit [fJ].
            %   algo = 'fft_search' / 'differential_kay' / 'fft_search_blind'
            %   fl_vec     — fractional bit widths (vector or scalar)
            %   blindD_vec — only used for 'fft_search_blind'
            % Output shape:
            %   non-blind: column vector length(fl_vec)
            %   blind:     [length(fl_vec) x length(blindD_vec)]
            switch algo
                case 'fft_search'
                    [NM, NA] = bit_width_full.fft_search_counts(P.TrainingLen, P.FR_Nfft, false);
                    NM = NM / P.SubframeLen;  NA = NA / P.SubframeLen;
                    E  = arrayfun(@(n) energy.receiver(NA, NM, P.EAdd_fJ, P.EMult_fJ, ...
                                  P.M, P.Oversampling, n), fl_vec);
                    E  = E(:);
                case 'differential_kay'
                    [NM, NA] = bit_width_full.differential_kay_counts(P.TrainingLen);
                    NM = NM / P.SubframeLen;  NA = NA / P.SubframeLen;
                    E  = arrayfun(@(n) energy.receiver(NA, NM, P.EAdd_fJ, P.EMult_fJ, ...
                                  P.M, P.Oversampling, n), fl_vec);
                    E  = E(:);
                case 'fft_search_blind'
                    E = nan(length(fl_vec), length(blindD_vec));
                    for di = 1:length(blindD_vec)
                        [NM, NA] = bit_width_full.fft_search_counts(blindD_vec(di), P.FR_Nfft, true);
                        NM = NM / P.SubframeLen;  NA = NA / P.SubframeLen;
                        E(:, di) = arrayfun(@(n) energy.receiver(NA, NM, P.EAdd_fJ, P.EMult_fJ, ...
                                            P.M, P.Oversampling, n), fl_vec);
                    end
                otherwise
                    error('bit_width_full:unknownFR', 'Unknown FR algo: %s', algo);
            end
        end

        function E = crEnergy(P, algo, fl_vec)
            switch algo
                case 'viterbi_viterbi'
                    [NM, NA] = bit_width_full.viterbi_counts(P.BlockLen);
                case 'pilots_only'
                    [NM, NA] = bit_width_full.pilots_only_counts();
                otherwise
                    error('bit_width_full:unknownCR', 'Unknown CR algo: %s', algo);
            end
            NM = NM / P.BlockLen;  NA = NA / P.BlockLen;
            E  = arrayfun(@(n) energy.receiver(NA, NM, P.EAdd_fJ, P.EMult_fJ, ...
                          P.M, P.Oversampling, n), fl_vec);
            E  = E(:);
        end

        % ---------------- Operation counts (from .tex tables) --------

        function [NM, NA] = fft_search_counts(L, Nfft, blind)
            if blind
                NM_form = 12 * L;
                NA_form = 9  * L;
            else
                NM_form = 4 * L;
                NA_form = 3 * L;
            end
            NM_fft    = 2 * L * log2(Nfft / L) + 2 * Nfft * log2(L);
            NA_fft    = 3 * L * log2(Nfft / L) + 3 * Nfft * log2(L);
            NM_search = 2 * Nfft;
            NA_search = 2 * Nfft;
            NM_interp = 5;
            NA_interp = 4;
            NM = NM_form + NM_fft + NM_search + NM_interp;
            NA = NA_form + NA_fft + NA_search + NA_interp;
        end

        function [NM, NA] = differential_kay_counts(L)
            NM = 7 * L + 1;
            NA = 4 * L - 2;
        end

        function [NM, NA] = viterbi_counts(N)
            NM = 16 * N + 1;
            NA = 14 * N - 1;
        end

        function [NM, NA] = pilots_only_counts()
            NM = 1;
            NA = 1;
        end

        % ---------------- MEX builds (per-FL temp dirs) ---------------

        function tempDir = buildFRMex(P, fxp_fr)
            clear mex %#ok<CLMEX>

            B.Rs            = P.Rs;
            B.N_pol         = P.N_pol;
            B.TrainingLen   = P.TrainingLen;
            B.FR_Nfft       = P.FR_Nfft;
            B.FR_Po2Twiddle = P.FR_Po2Twiddle;
            B.FR_BlindD     = P.BlindD;
            B.FxpConfig_FR  = fxp_fr;
            B.CordicIts     = P.CordicIts;
            B.MaxFreq       = P.MaxFreq;

            cfg = coder.config('mex');
            cfg.GenerateReport = false;

            build_freq_recovery_fft_search_fxp_mex(B, cfg);
            build_freq_recovery_differential_kay_fxp_mex(B, cfg);

            srcDir  = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src');
            tempDir = tempname;
            dstPkg  = fullfile(tempDir, '+freq_recovery');
            mkdir(dstPkg);
            ext = mexext;
            for f = {'fft_search_fxp_mex', 'differential_kay_fxp_mex'}
                copyfile( ...
                    fullfile(srcDir, '+freq_recovery', [f{1} '.' ext]), ...
                    fullfile(dstPkg,                   [f{1} '.' ext]));
            end
        end

        function buildCRMexInSrc(P, fxp_cr)
            % Build CR MEX directly into src/+carrier_recovery (no temp
            % dir).  Used when only one CR MEX is active at a time and is
            % swapped between outer iterations of the grid sweep.
            clear mex %#ok<CLMEX>

            B.Rs             = P.Rs;
            B.N_pol          = P.N_pol;
            B.FxpConfig_VV   = fxp_cr;
            B.FxpConfig_PO   = fxp_cr;
            B.CordicIts      = P.CordicIts;
            B.VV_NTaps       = P.VV_NTaps;
            B.BlockLen       = P.BlockLen;
            B.StepSize       = P.StepSize;
            B.PilotThreshold = P.PilotThreshold;
            B.PilotLen       = 1;

            cfg = coder.config('mex');
            cfg.GenerateReport = false;

            build_carrier_recovery_viterbiViterbi_fxp_mex(B, cfg);
            build_carrier_recovery_pilots_only_fxp_mex(B, cfg);
        end

        % ---------------- Channel + SNR sweep -------------------------

        function Params = extractParams(testCase)
            Params.Rs             = testCase.Rs;
            Params.N_pol          = testCase.N_pol;
            Params.NTrials        = testCase.NTrials;
            Params.SNR_dB_vec     = testCase.SNR_dB_vec;
            Params.DeltaF_Hz      = testCase.DeltaF_Hz;
            Params.LW_Hz          = testCase.LW_Hz;
            Params.FR_Nfft        = testCase.FR_Nfft;
            Params.FR_Po2Twiddle  = testCase.FR_Po2Twiddle;
            Params.MaxFreq        = testCase.MaxFreq;
            Params.BlindD         = testCase.BlindD;
            Params.CordicIts      = testCase.CordicIts;
            Params.BlockLen       = testCase.BlockLen;
            Params.StepSize       = testCase.StepSize;
            Params.PilotThreshold = testCase.PilotThreshold;
            Params.VV_NTaps       = testCase.VV_NTaps;
            Params.FEC_BER        = testCase.FEC_BER;
            Params.VVFilters      = testCase.VVFilters;
        end

        function ber_all = runSnrSweepStatic(Params, fr_algo, blindD, T_fr, T_cr)
            % Returns per-trial BER: ber_all(tr, si, c) for trial tr, SNR
            % index si and CR column c (1 = VV, 2 = PO).  Callers compute
            % the FEC-SNR mean and std from per-trial crossings.
            NSNR    = length(Params.SNR_dB_vec);
            ber_all = zeros(Params.NTrials, NSNR, 2);

            for tr = 1:Params.NTrials
                for si = 1:NSNR
                    [fr_out, pilots, txRefBits] = bit_width_full.buildChannel( ...
                        Params, Params.SNR_dB_vec(si), fr_algo, blindD, T_fr);

                    fr_fi     = cast(fr_out,                    'like', T_cr.x);
                    pilots_fi = cast(pilots,                     'like', T_cr.x);
                    vvfilt_fi = cast(Params.VVFilters{si},       'like', T_cr.w);

                    [cr_vv, ~] = carrier_recovery.viterbiViterbi_fxp_mex( ...
                        fr_fi, Params.N_pol, Params.VV_NTaps, vvfilt_fi, pilots_fi, ...
                        Params.BlockLen, double(Params.StepSize), Params.PilotThreshold, ...
                        double(Params.CordicIts), T_cr);
                    cr_vv = bit_width_full.resolveAmbiguity(double(cr_vv), txRefBits);
                    ber_all(tr, si, 1) = bit_width_full.computeBER(cr_vv, txRefBits);

                    [cr_po, ~] = carrier_recovery.pilots_only_fxp_mex( ...
                        fr_fi, Params.N_pol, Params.BlockLen, pilots_fi, ...
                        double(Params.CordicIts), T_cr);
                    cr_po = bit_width_full.resolveAmbiguity(double(cr_po), txRefBits);
                    ber_all(tr, si, 2) = bit_width_full.computeBER(cr_po, txRefBits);
                end
            end
        end

        function [fr_out, pilots, txRefBits] = buildChannel(P, SNR_dB, fr_algo, blindD, T_fr)
            CPON_BLOCK_LEN = 32;
            CPON_SF_SYMS   = 3712;
            N_CPON_BLOCKS  = 116;
            BITS_PER_SF    = 3586 * 2 * 2;

            txBits = modem.randomBits(BITS_PER_SF);
            [symbols, pilotSyms, training, ~] = modem.modulate(txBits);

            rx = channel.lo_freq_shift(symbols, P.DeltaF_Hz / 1e6, P.Rs, 1);
            rx = channel.add_awgn(rx, SNR_dB);
            rx = channel.add_phase_noise(rx, P.Rs, P.LW_Hz);

            rx_fi = cast(rx,       'like', T_fr.x);
            tr_fi = cast(training, 'like', T_fr.x);

            switch fr_algo
                case 'fft_search'
                    [fr_fi, ~] = freq_recovery.fft_search_fxp_mex( ...
                        rx_fi, tr_fi, P.Rs, P.FR_Nfft, P.FR_Po2Twiddle, ...
                        P.CordicIts, P.MaxFreq, T_fr, true, 0);
                case 'fft_search_blind'
                    [fr_fi, ~] = freq_recovery.fft_search_fxp_mex( ...
                        rx_fi, tr_fi, P.Rs, P.FR_Nfft, P.FR_Po2Twiddle, ...
                        P.CordicIts, P.MaxFreq, T_fr, false, blindD);
                case 'differential_kay'
                    [fr_fi, ~] = freq_recovery.differential_kay_fxp_mex( ...
                        rx_fi, tr_fi, P.Rs, P.CordicIts, T_fr, true, 0, P.MaxFreq);
                otherwise
                    error('bit_width_full:unknownFR', 'Unknown FR algorithm: %s', fr_algo);
            end

            fr_out = double(fr_fi);

            Nsym    = size(fr_out, 1);
            NBlocks = ceil(Nsym / P.BlockLen);
            pilots  = zeros(NBlocks, P.N_pol);
            for b = 1:NBlocks
                pos     = (b - 1) * P.BlockLen + 1;
                posInSf = mod(pos - 1, CPON_SF_SYMS) + 1;
                blk     = min(floor((posInSf - 1) / CPON_BLOCK_LEN) + 1, N_CPON_BLOCKS);
                pilots(b, :) = pilotSyms(blk, :);
            end

            txRefBits = modem.symbolsToBits(symbols);
        end

        function BER = computeBER(crSym, refBits)
            dec  = modem.decideSymbols(crSym);
            bits = modem.symbolsToBits(dec);
            n    = min(length(refBits), length(bits));
            BER  = sum(refBits(1:n) ~= bits(1:n)) / n;
        end

        function best = resolveAmbiguity(crSym, txRefBits)
            bestBER = Inf;
            best    = crSym;
            for k = 0:3
                rot  = crSym .* exp(-1j * k * pi/2);
                dec  = modem.decideSymbols(rot);
                bits = modem.symbolsToBits(dec);
                n    = min(length(txRefBits), length(bits));
                ber  = sum(txRefBits(1:n) ~= bits(1:n)) / n;
                if ber < bestBER
                    bestBER = ber;
                    best    = rot;
                end
            end
        end

        function snrs = fecCrossingsAll(snrDb, berPerTrial, fecBer)
            % Per-trial FEC SNR crossings. berPerTrial is [NTrials × NSNR].
            % Returns NTrials×1 vector; NaN where the curve never crosses fecBer.
            NTrials = size(berPerTrial, 1);
            snrs    = nan(NTrials, 1);
            for t = 1:NTrials
                snrs(t) = bit_width_full.fecCrossing(snrDb, berPerTrial(t, :), fecBer);
            end
        end

        function xCross = fecCrossing(x, y, yLimit)
            xCross = NaN;
            n      = length(x);
            if n < 2, return; end

            exact = find(y == yLimit, 1, 'first');
            if ~isempty(exact)
                xCross = x(exact);
                return;
            end

            for i = 1:(n - 1)
                if (y(i) - yLimit) * (y(i+1) - yLimit) < 0
                    xCross = x(i) + (yLimit - y(i)) * (x(i+1) - x(i)) / (y(i+1) - y(i));
                    return;
                end
            end
        end

    end
end
