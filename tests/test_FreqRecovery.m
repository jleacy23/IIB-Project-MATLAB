classdef test_FreqRecovery < matlab.unittest.TestCase
%TEST_FREQRECOVERY  Fixed-point MEX frequency-recovery tests.
%
%   MEX tests for the CORDIC-based fixed-point frequency-recovery estimators
%       freq_recovery.fft_search_fxp        (FFT peak search)
%       freq_recovery.differential_kay_fxp  (differential + Tretter/Kay)
%
%       modem.modulate  ->  CFO (lo_freq_shift)  ->  AWGN  ->  modem.normalise
%                       ->  fixed-point FR (CORDIC)  ->  estimated offset
%
%   A known carrier-frequency offset plus AWGN are applied (no phase noise,
%   to isolate the frequency estimator).  The estimator is run once per
%   subframe (each subframe carries its own CPON training preamble) and the
%   per-subframe frequency estimates are averaged over NumSubframes; the
%   averaged estimate is checked against the applied offset within FreqTol_MHz.
%
%   Fixed-point configuration
%     The FR type table pins T.theta / T.acc wide (>= 16 fraction bits) and
%     sweeps only the signal type T.x, so the CORDIC iteration count is held
%     at CordicIts (matched to the wide angle precision), independent of T.x.
%
%   Execution mode (UseMex)
%     UseMex = true  -> builds and calls the compiled *_fxp_mex binaries
%                       (requires MATLAB Coder + Fixed-Point Designer).
%     UseMex = false -> calls the *_fxp functions directly (interpreted fi;
%                       identical arithmetic, no toolbox/build needed).
%
%   Prerequisites (UseMex = true)
%     - MATLAB Coder and Fixed-Point Designer toolboxes licensed.
%     - build_freq_recovery_fft_search_fxp_mex.m and
%       build_freq_recovery_differential_kay_fxp_mex.m on the path.

    properties (Constant)
        % ---- Signal -------------------------------------------------
        N_pol        = 2
        Rs           = 30.5          % symbol rate [GBd]
        TrainingLen  = 11            % CPON training symbols per subframe
        SubframeLen  = 3712          % symbols per CPON subframe / pol
        NumSubframes = 10             % subframes to estimate over and average

        % ---- Channel ------------------------------------------------
        SNR_dB     = 20            % good SNR to isolate FR errors
        DeltaF_MHz = 2000          % applied carrier-frequency offset [MHz]
        NormPct    = 99.9          % modem.normalise percentile

        % ---- Fixed-point config -------------------------------------
        FxpConfig     = 'fixed32'    % 'fixed16' | 'fixed32' | struct('WL',wl,'FL',fl)
        CordicIts     = 16           % CORDIC iterations (matched to wide angle FL)
        FR_Nfft       = 2048         % FFT size (power of 2 >= TrainingLen)
        FR_Po2Twiddle = false
        FR_BlindD     = 64           % blind data length baked into the MEX type
        MaxFreq       = 1            % phase-scaling factor

        % ---- Pass / fail --------------------------------------------
        % Data-aided fft_search estimates from only the 11 CPON training
        % symbols, so it is inherently coarse (tens of MHz); the tolerance is
        % sized for that.  differential_kay is far more accurate (~10 MHz).
        FreqTol_MHz = 200            % allowed |estimate - applied| [MHz]

        % ---- Execution / build --------------------------------------
        UseMex  = true
        Rebuild = true
    end

    properties
        Tfr    % FR fixed-point type table
    end

    %% ================================================================
    %  One-time setup
    %% ================================================================
    methods (TestClassSetup)

        function setupPath(~)
            here = fileparts(mfilename('fullpath'));
            addpath(genpath(fullfile(here, '..', 'src')));
            addpath(fullfile(here, '..', 'build'));
        end

        function seedRng(~)
            rng(20260530);
        end

        function setupTypes(testCase)
            testCase.Tfr = freq_recovery.fxp_types(testCase.FxpConfig);
        end

        function buildMex(testCase)
            % Compile both MEX binaries so their baked-in fixed-point type and
            % CORDIC-iteration constant match this test's configuration exactly.
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
            P.FR_BlindD     = testCase.FR_BlindD;
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
    %  Tests
    %% ================================================================
    methods (Test)

        function testFFTSearch_Fxp(testCase)
            [rxSym, training] = buildChannel(testCase);
            [dF, dFs] = averageOverSubframes(testCase, ...
                @(rx, tr) runFFTSearch(testCase, rx, tr), rxSym, training);
            reportEstimate(testCase, 'FFT-search', dF, dFs);
        end

        function testDifferentialKay_Fxp(testCase)
            [rxSym, training] = buildChannel(testCase);
            [dF, dFs] = averageOverSubframes(testCase, ...
                @(rx, tr) runDifferentialKay(testCase, rx, tr), rxSym, training);
            reportEstimate(testCase, 'Differential-Kay', dF, dFs);
        end

    end

    %% ================================================================
    %  Private helpers
    %% ================================================================
    methods (Access = private)

        % ---- Channel: modulate -> CFO -> AWGN -----------------------
        function [rxSym, training] = buildChannel(testCase)
            BITS_PER_SF = 3586 * 2 * 2;
            txBits = modem.randomBits(testCase.NumSubframes * BITS_PER_SF);
            [symbols, ~, training, ~] = modem.modulate(txBits);

            % Frequency offset then AWGN (no phase noise, to isolate FR).
            rxSym = channel.lo_freq_shift(symbols, testCase.DeltaF_MHz, testCase.Rs, 1);
            rxSym = channel.add_awgn(rxSym, testCase.SNR_dB);

            % Normalise into the unit box before the fixed-point cast.  The
            % angle-domain estimators are scale-invariant, so the ±1±1j
            % training reference need not be rescaled.
            rxSym = modem.normalise(rxSym, testCase.NormPct);
        end

        % ---- Per-subframe estimate, averaged over NumSubframes ------
        function [dFmean, dFs] = averageOverSubframes(testCase, estimatorFn, rxSym, training)
            % Run the estimator independently on each subframe (each starts
            % with its own CPON training preamble) and average the estimates.
            SF  = testCase.SubframeLen;
            nSF = min(testCase.NumSubframes, floor(size(rxSym, 1) / SF));
            dFs = zeros(nSF, 1);
            for sf = 1:nSF
                idx      = (sf - 1) * SF + (1:SF);
                dFs(sf)  = estimatorFn(rxSym(idx, :), training);
            end
            dFmean = mean(dFs);
        end

        % ---- fft_search runner (MEX or interpreted) -----------------
        function dF = runFFTSearch(testCase, rxSym, training)
            T = testCase.Tfr;
            rx_fi = cast(rxSym,    'like', T.x);
            tr_fi = cast(training, 'like', T.x);

            if testCase.UseMex
                [~, dF] = freq_recovery.fft_search_fxp_mex( ...
                    rx_fi, tr_fi, testCase.Rs, double(testCase.FR_Nfft), ...
                    logical(testCase.FR_Po2Twiddle), double(testCase.CordicIts), ...
                    double(testCase.MaxFreq), T, true, 0);
            else
                [~, dF] = freq_recovery.fft_search_fxp( ...
                    rx_fi, tr_fi, testCase.Rs, double(testCase.FR_Nfft), ...
                    logical(testCase.FR_Po2Twiddle), double(testCase.CordicIts), ...
                    double(testCase.MaxFreq), T, true, 0);
            end
        end

        % ---- differential_kay runner (MEX or interpreted) -----------
        function dF = runDifferentialKay(testCase, rxSym, training)
            T = testCase.Tfr;
            rx_fi = cast(rxSym,    'like', T.x);
            tr_fi = cast(training, 'like', T.x);

            if testCase.UseMex
                [~, dF] = freq_recovery.differential_kay_fxp_mex( ...
                    rx_fi, tr_fi, testCase.Rs, double(testCase.CordicIts), ...
                    T, true, 0, double(testCase.MaxFreq));
            else
                [~, dF] = freq_recovery.differential_kay_fxp( ...
                    rx_fi, tr_fi, testCase.Rs, double(testCase.CordicIts), ...
                    T, true, 0, double(testCase.MaxFreq));
            end
        end

        % ---- Report + check the (averaged) estimated frequency ------
        function reportEstimate(testCase, algoName, dF, dFs)
            dF_MHz  = dF / 1e6;
            err_MHz = dF_MHz - testCase.DeltaF_MHz;
            std_MHz = std(dFs(:) / 1e6);
            fprintf(['%s fxp (%s): estimated %.3f MHz (avg over %d subframes, ' ...
                     'std %.3f MHz) | applied %.3f MHz | error %+.3f MHz\n'], ...
                algoName, cfgName(testCase), dF_MHz, numel(dFs), std_MHz, ...
                testCase.DeltaF_MHz, err_MHz);

            testCase.verifyTrue(all(isfinite(dFs)), ...
                sprintf('%s fxp produced a non-finite frequency estimate.', algoName));
            testCase.verifyLessThanOrEqual(abs(err_MHz), testCase.FreqTol_MHz, ...
                sprintf('%s fxp averaged estimate off by %.3f MHz (tol %.3f MHz).', ...
                    algoName, err_MHz, testCase.FreqTol_MHz));
        end

    end
end

%% ====================================================================
%  Local function
%% ====================================================================
function s = cfgName(testCase)
    if ischar(testCase.FxpConfig) || isstring(testCase.FxpConfig)
        s = char(testCase.FxpConfig);
    else
        s = sprintf('WL=%d,FL=%d', testCase.FxpConfig.WL, testCase.FxpConfig.FL);
    end
end
