classdef subframe_averaging < matlab.unittest.TestCase
%SUBFRAME_AVERAGING  Subframes-to-average for a target CFO estimate, + op counts.
%
%   Compares three high-precision fixed-point frequency-recovery configurations
%       1. differential_kay  (data-aided, L = 11 training symbols)
%       2. fft_search        (blind, D = N_FFT = 2048 — largest FFT in a subframe)
%       3. fft_search        (data-aided, L = 11, zero-padded to N_FFT = 2048)
%
%   For each, the CFO is estimated once per subframe and the estimates are
%   averaged.  The test finds the smallest number of subframes N that must be
%   averaged so the RMS CFO-estimation error falls below TargetErr_Hz — the
%   residual CFO the downstream phase recovery can absorb (10 MHz here; the
%   channel still carries 1 MHz Wiener phase noise).  It then reports the
%   operation counts of each algorithm (from the report's cost tables
%   tab:fft_cost / tab:diffkay_cost) multiplied by that N — i.e. the total
%   arithmetic to reach the target estimate.
%
%   Why the three behave differently
%     - Both data-aided variants use only the 11 fixed CPON training symbols, so
%       their per-subframe variance is large.  fft_search additionally carries a
%       systematic FFT leakage / interpolation bias that is identical every
%       subframe (the training never changes) and therefore does NOT average
%       away, so its RMS error can floor above the target (N -> infeasible).
%     - Blind fft_search uses D = 2048 random data symbols, giving a far smaller
%       per-subframe variance (and a data-dependent bias that does average out),
%       reaching the target in very few subframes but at a high per-estimate cost.
%
%   Channel: CFO (random per trial) + 1 MHz Wiener phase noise + AWGN, then
%   modem.normalise.  Phase recovery / decoding are not involved — only the
%   frequency estimate is assessed.
%
%   Output: subframe_averaging.mat (table 'tbl') + a printed summary.
%
%   Run:  runtests('subframe_averaging')

    %% ================================================================
    %  Parameters
    %% ================================================================
    properties (Constant)
        % ---- System -------------------------------------------------
        Rs          = 30.5            % symbol rate [GBd]
        N_pol       = 2
        TrainingLen = 11              % CPON training symbols (data-aided L)
        SubframeLen = 3712

        % ---- Channel ------------------------------------------------
        SNR_dB        = 20            % operating SNR for the estimate
        LW_Hz         = 1e6           % 1 MHz laser linewidth (phase noise)
        DeltaF_GHz    = [0.5, 3.0]    % CFO drawn uniformly in this range per trial
        NormPct       = 99.9          % modem.normalise percentile

        % ---- Target -------------------------------------------------
        TargetErr_Hz  = 10e6          % required RMS CFO error (residual the phase recovery can absorb)

        % ---- FFT / blind --------------------------------------------
        FR_Nfft       = 2048          % largest power-of-2 FFT inside a subframe
        FR_Po2Twiddle = false
        BlindD        = 2048          % blind observation length (= N_FFT)

        % ---- Fixed-point (high precision) ---------------------------
        FxpConfig = 'fixed32'         % 32-bit, 16 fractional bits
        CordicIts = 16
        MaxFreq   = 1

        % ---- Monte-Carlo --------------------------------------------
        NTrials   = 30                % independent CFO/noise realisations
        % Max subframes searched per algorithm (data-aided are cheap to run for
        % many subframes; the 2048-pt FFT variants are capped lower).
        Nmax_diffkay   = 4096
        Nmax_fft_data  = 256
        Nmax_fft_blind = 128

        % ---- Execution / build --------------------------------------
        UseMex  = true
        Rebuild = false
    end

    properties
        Tfr        % FR fixed-point type table
        Training   % [11 x 2] CPON training symbols
    end

    %% ================================================================
    %  Setup
    %% ================================================================
    methods (TestClassSetup)

        function setupPath(~)
            here = fileparts(mfilename('fullpath'));
            addpath(genpath(fullfile(here, '..', '..', 'src')));
            addpath(fullfile(here, '..', '..', 'build'));
        end

        function seedRng(~)
            rng(20260530);
        end

        function setupTypesAndTraining(testCase)
            testCase.Tfr = freq_recovery.fxp_types(testCase.FxpConfig);
            % CPON training preamble (identical every subframe).
            [~, ~, training, ~] = modem.modulate(modem.randomBits(3586 * 2 * 2));
            testCase.Training = training;
        end

        function buildMex(testCase)
            if ~testCase.UseMex || ~testCase.Rebuild
                return;
            end
            cfg = coder.config('mex');
            cfg.GenerateReport       = false;
            cfg.IntegrityChecks      = false;
            cfg.ResponsivenessChecks = false;

            P.N_pol         = testCase.N_pol;
            P.TrainingLen   = testCase.TrainingLen;
            P.Rs            = testCase.Rs;
            P.FR_Nfft       = testCase.FR_Nfft;
            P.FR_Po2Twiddle = testCase.FR_Po2Twiddle;
            P.FR_BlindD     = testCase.BlindD;
            P.FxpConfig_FR  = testCase.FxpConfig;
            P.CordicIts     = testCase.CordicIts;
            P.MaxFreq       = testCase.MaxFreq;

            fprintf('  Building fft_search_fxp_mex...\n');
            build_freq_recovery_fft_search_fxp_mex(P, cfg);
            fprintf('  Building differential_kay_fxp_mex...\n');
            build_freq_recovery_differential_kay_fxp_mex(P, cfg);
        end

    end

    %% ================================================================
    %  Test
    %% ================================================================
    methods (Test)

        function test_subframe_averaging(testCase)
            L    = testCase.TrainingLen;
            Nfft = testCase.FR_Nfft;
            D    = testCase.BlindD;

            % --- Algorithm definitions ----------------------------------
            % mode      : runner switch
            % Nmax      : max subframes to search
            % [NM, NA]  : per-subframe op counts from the report cost tables
            algos = { ...
                struct('name', "diffkay_data", 'mode', 'dk_data', ...
                       'Nmax', testCase.Nmax_diffkay, ...
                       'counts', {bit_width_full_counts('diffkay', L, Nfft, D)}); ...
                struct('name', "fft_blind",    'mode', 'fft_blind', ...
                       'Nmax', testCase.Nmax_fft_blind, ...
                       'counts', {bit_width_full_counts('fft_blind', L, Nfft, D)}); ...
                struct('name', "fft_data",     'mode', 'fft_data', ...
                       'Nmax', testCase.Nmax_fft_data, ...
                       'counts', {bit_width_full_counts('fft_data', L, Nfft, D)}) };

            nAlgo  = numel(algos);
            name      = strings(nAlgo, 1);
            N_needed  = nan(nAlgo, 1);
            achievable = false(nAlgo, 1);
            rmse_floor_MHz = nan(nAlgo, 1);   % RMS error at Nmax (residual after averaging)
            RM_sub   = nan(nAlgo, 1);
            RA_sub   = nan(nAlgo, 1);
            RM_total = nan(nAlgo, 1);
            RA_total = nan(nAlgo, 1);

            for a = 1:nAlgo
                A = algos{a};
                fprintf('=== %s ===\n', A.name);
                [Nreq, rmseVec] = testCase.subframesToTarget(A.mode, A.Nmax);

                name(a)           = A.name;
                rmse_floor_MHz(a) = rmseVec(end) / 1e6;
                RM_sub(a)         = A.counts(1);
                RA_sub(a)         = A.counts(2);

                if isfinite(Nreq)
                    achievable(a) = true;
                    N_needed(a)   = Nreq;
                    RM_total(a)   = A.counts(1) * Nreq;
                    RA_total(a)   = A.counts(2) * Nreq;
                    fprintf('  reaches %.3f MHz error in %d subframes\n', ...
                        testCase.TargetErr_Hz / 1e6, Nreq);
                else
                    fprintf('  does NOT reach %.3f MHz within %d subframes (floor %.3f MHz)\n', ...
                        testCase.TargetErr_Hz / 1e6, A.Nmax, rmse_floor_MHz(a));
                end
            end

            % --- Operation-count table ----------------------------------
            tbl = table(name, N_needed, achievable, rmse_floor_MHz, ...
                        RM_sub, RA_sub, RM_total, RA_total);
            outDir  = fileparts(mfilename('fullpath'));
            outFile = fullfile(outDir, 'subframe_averaging.mat');
            save(outFile, 'tbl');
            fprintf('\nSaved results to %s\n\n', outFile);
            disp(tbl);

            % --- Sanity checks ------------------------------------------
            % The variance-limited estimators must reach the target by averaging.
            testCase.verifyTrue(achievable(strcmp(name, "diffkay_data")), ...
                'differential_kay should reach the target by averaging.');
            testCase.verifyTrue(achievable(strcmp(name, "fft_blind")), ...
                'blind fft_search should reach the target by averaging.');
        end

    end

    %% ================================================================
    %  Private helpers
    %% ================================================================
    methods (Access = private)

        function [N_needed, rmse] = subframesToTarget(testCase, mode, Nmax)
            % Monte-Carlo: for each trial fix a random CFO, estimate it once per
            % subframe (independent noise), cumulatively average the estimates,
            % and record the error vs the number of averaged subframes.  RMSE(N)
            % is taken over trials; N_needed is the first N with RMSE <= target.
            NT   = testCase.NTrials;
            errN = nan(NT, Nmax);             % cumulative-average error per (trial, N)

            for t = 1:NT
                fLo = testCase.DeltaF_GHz(1) * 1e9;
                fHi = testCase.DeltaF_GHz(2) * 1e9;
                f_true = fLo + (fHi - fLo) * rand;

                est = nan(Nmax, 1);
                for sf = 1:Nmax
                    [rx, tr] = testCase.buildObservation(mode, f_true);
                    est(sf)  = testCase.estimate(mode, rx, tr);
                end
                cumAvg     = cumsum(est) ./ (1:Nmax).';
                errN(t, :) = (cumAvg - f_true).';
            end

            rmse     = sqrt(mean(errN.^2, 1)).';          % [Nmax x 1], Hz
            idx      = find(rmse <= testCase.TargetErr_Hz, 1, 'first');
            if isempty(idx)
                N_needed = Inf;
            else
                N_needed = idx;
            end
        end

        function [rx, training] = buildObservation(testCase, mode, f_true_Hz)
            % Build one subframe's observation window for the estimator.
            %   data-aided : the 11 CPON training symbols (identical each subframe)
            %   blind      : training + BlindD random QPSK data symbols
            training = testCase.Training;

            switch mode
                case {'dk_data', 'fft_data'}
                    syms = training;                       % [L x 2]
                case 'fft_blind'
                    D    = testCase.BlindD;
                    data = (2*randi([0 1], D, testCase.N_pol) - 1) + ...
                       1j*(2*randi([0 1], D, testCase.N_pol) - 1);
                    syms = [training; data];               % [(L+D) x 2]
                otherwise
                    error('subframe_averaging:badMode', 'Unknown mode %s', mode);
            end

            rx = channel.lo_freq_shift(syms, f_true_Hz / 1e6, testCase.Rs, 1);
            rx = channel.add_phase_noise(rx, testCase.Rs, testCase.LW_Hz);
            rx = channel.add_awgn(rx, testCase.SNR_dB);
            rx = modem.normalise(rx, testCase.NormPct);
        end

        function dF = estimate(testCase, mode, rx, training)
            T     = testCase.Tfr;
            rx_fi = cast(rx,       'like', T.x);
            tr_fi = cast(training, 'like', T.x);
            rsGBd = testCase.Rs;
            ci    = double(testCase.CordicIts);
            mf    = double(testCase.MaxFreq);
            Nfft  = double(testCase.FR_Nfft);
            po2   = logical(testCase.FR_Po2Twiddle);
            D     = double(testCase.BlindD);

            if testCase.UseMex
                switch mode
                    case 'dk_data'
                        [~, dF] = freq_recovery.differential_kay_fxp_mex( ...
                            rx_fi, tr_fi, rsGBd, ci, T, true, 0, mf);
                    case 'fft_data'
                        [~, dF] = freq_recovery.fft_search_fxp_mex( ...
                            rx_fi, tr_fi, rsGBd, Nfft, po2, ci, mf, T, true, 0);
                    case 'fft_blind'
                        [~, dF] = freq_recovery.fft_search_fxp_mex( ...
                            rx_fi, tr_fi, rsGBd, Nfft, po2, ci, mf, T, false, D);
                end
            else
                switch mode
                    case 'dk_data'
                        [~, dF] = freq_recovery.differential_kay_fxp( ...
                            rx_fi, tr_fi, rsGBd, ci, T, true, 0, mf);
                    case 'fft_data'
                        [~, dF] = freq_recovery.fft_search_fxp( ...
                            rx_fi, tr_fi, rsGBd, Nfft, po2, ci, mf, T, true, 0);
                    case 'fft_blind'
                        [~, dF] = freq_recovery.fft_search_fxp( ...
                            rx_fi, tr_fi, rsGBd, Nfft, po2, ci, mf, T, false, D);
                end
            end
        end

    end
end

%% ====================================================================
%  Local: per-subframe operation counts from the report cost tables
%  (tab:fft_cost / tab:diffkay_cost).  Returns [NM, NA] = [mults, adds].
%% ====================================================================
function c = bit_width_full_counts(which, L, Nfft, D)
    switch which
        case 'diffkay'      % tab:diffkay_cost, data-aided L
            NM = 7*L + 1;
            NA = 4*L - 2;
        case 'fft_data'     % tab:fft_cost, data-aided (L training, padded to Nfft)
            [NM, NA] = fft_counts(L, Nfft, false);
        case 'fft_blind'    % tab:fft_cost, blind (observation length D, FFT Nfft)
            [NM, NA] = fft_counts(D, Nfft, true);
        otherwise
            error('subframe_averaging:badCounts', 'Unknown algo %s', which);
    end
    c = [NM, NA];
end

function [NM, NA] = fft_counts(Lobs, Nfft, blind)
    if blind
        NM_form = 12 * Lobs;   NA_form = 9 * Lobs;
    else
        NM_form = 4 * Lobs;    NA_form = 3 * Lobs;
    end
    NM_fft    = 2*Lobs*log2(Nfft/Lobs) + 2*Nfft*log2(Lobs);
    NA_fft    = 3*Lobs*log2(Nfft/Lobs) + 3*Nfft*log2(Lobs);
    NM_search = 2*Nfft;        NA_search = 2*Nfft;
    NM_interp = 5;             NA_interp = 4;
    NM = NM_form + NM_fft + NM_search + NM_interp;
    NA = NA_form + NA_fft + NA_search + NA_interp;
end
