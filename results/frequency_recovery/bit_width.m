classdef bit_width < matlab.unittest.TestCase
%BIT_WIDTH  SNR threshold at FEC limit vs fractional bit width.
%
%   Two test methods, each producing three figures (one per FR variant):
%
%   test_snr_vs_fr_bitwidth  — sweeps the FR fractional bit width (FL)
%       while the CR bit width is held at FL_fixed.  Three FR variants:
%       data-aided FFT search, data-aided differential Kay, and blind FFT
%       search (D = BlindD, Nfft = FR_Nfft).  Each figure shows two lines,
%       one per phase recovery algorithm (Viterbi-Viterbi / pilots-only).
%
%   test_snr_vs_cr_bitwidth  — same structure but sweeps CR FL with FR
%       held fixed.
%
%   Parallelism
%     Each FL iteration is compiled into a private temp directory that has
%     the same +freq_recovery / +carrier_recovery package structure as src/.
%     A parfor loop then runs the SNR sweeps concurrently: each worker calls
%     addpath on its own temp directory, so it uses the MEX built for that
%     specific FL value without interfering with other workers.
%
%   Run with:
%       runtests('bit_width')
%       runtests('bit_width', 'ProcedureName', 'test_snr_vs_fr_bitwidth')

    %% ================================================================
    %  Parameters
    %% ================================================================
    properties (Constant)

        % System
        Rs          = 30.5          % symbol rate [GBd]
        N_pol       = 2
        TrainingLen = 11

        % Monte-Carlo
        NTrials     = 100           % trials per (SNR, FL) point

        % SNR sweep
        SNR_dB_vec  = 0 : 1 : 30   % [dB]

        % Bit width sweep — word length is fixed; fractional length varies
        WL          = 32
        FL_vec      = [4, 6, 8, 10, 12, 14, 16]  % fractional bit widths
        FL_fixed    = 16            % held constant for the non-swept subsystem

        % Channel conditions
        DeltaF_Hz   = 1.5e9           % frequency offset [Hz]
        LW_Hz       = 1000e3        % laser linewidth  [Hz]

        % FFT search parameters — Nfft constant for all runs
        FR_Nfft       = 512
        FR_Po2Twiddle = false
        MaxFreq       = 0.1

        % Blind FFT search observation length
        BlindD = 512

        % Phase recovery
        CordicIts      = 16
        BlockLen       = 32
        StepSize       = 32
        PilotThreshold = 5 * pi / 9
        VV_NTaps       = 10

        % FEC threshold
        FEC_BER = 2e-2

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

        function test_snr_vs_fr_bitwidth(testCase)
            % Sweep FR fractional bit width; CR held at FL_fixed.
            P    = testCase;
            NFL  = length(P.FL_vec);
            Prms = bit_width.extractParams(testCase);

            fxp_cr_fixed = struct('WL', P.WL, 'FL', P.FL_fixed);
            T_cr = carrier_recovery.fxp_types(fxp_cr_fixed);

            % ---- Phase 1: Serial MEX builds into per-FL temp dirs -------
            fprintf('Building CR MEX (FL = %d, fixed)...\n', P.FL_fixed);
            crDir = bit_width.buildCRMex(P, fxp_cr_fixed);

            frDirs = cell(1, NFL);
            for bi = 1:NFL
                fl = P.FL_vec(bi);
                fprintf('[FR FL = %2d] building FR MEX  (%d / %d)\n', fl, bi, NFL);
                frDirs{bi} = bit_width.buildFRMex(P, struct('WL', P.WL, 'FL', fl));
            end

            % ---- Phase 2: Parallel SNR sweeps ---------------------------
            fecSNR_fft_DA    = nan(NFL, 2);
            fecSNR_dk_DA     = nan(NFL, 2);
            fecSNR_fft_blind = nan(NFL, 2);

            FL_vec_b   = P.FL_vec;
            WL_b       = P.WL;
            SNR_dB_b   = P.SNR_dB_vec;
            FEC_BER_b  = P.FEC_BER;
            BlindD_b   = P.BlindD;

            parfor pi = 1:NFL
                % Each worker prepends its own MEX dirs to its local path.
                % Workers are separate processes, so addpath is local.
                addpath(frDirs{pi});  %#ok<PFBNS>
                addpath(crDir);

                fl     = FL_vec_b(pi);
                T_fr_w = freq_recovery.fxp_types(struct('WL', WL_b, 'FL', fl));

                ber = bit_width.runSnrSweepStatic(Prms, 'fft_search', 0, T_fr_w, T_cr);
                fecSNR_fft_DA(pi, :) = [ ...
                    bit_width.fecCrossing(SNR_dB_b, ber(:,1), FEC_BER_b), ...
                    bit_width.fecCrossing(SNR_dB_b, ber(:,2), FEC_BER_b)];

                ber = bit_width.runSnrSweepStatic(Prms, 'differential_kay', 0, T_fr_w, T_cr);
                fecSNR_dk_DA(pi, :) = [ ...
                    bit_width.fecCrossing(SNR_dB_b, ber(:,1), FEC_BER_b), ...
                    bit_width.fecCrossing(SNR_dB_b, ber(:,2), FEC_BER_b)];

                ber = bit_width.runSnrSweepStatic(Prms, 'fft_search_blind', BlindD_b, T_fr_w, T_cr);
                fecSNR_fft_blind(pi, :) = [ ...
                    bit_width.fecCrossing(SNR_dB_b, ber(:,1), FEC_BER_b), ...
                    bit_width.fecCrossing(SNR_dB_b, ber(:,2), FEC_BER_b)];
            end

            bit_width.plotFRSweep(P, P.FL_vec, fecSNR_fft_DA, fecSNR_dk_DA, fecSNR_fft_blind);

            cellfun(@(d) rmdir(d, 's'), frDirs, 'UniformOutput', false);
            rmdir(crDir, 's');
        end

        function test_snr_vs_cr_bitwidth(testCase)
            % Sweep CR fractional bit width; FR held at FL_fixed.
            P    = testCase;
            NFL  = length(P.FL_vec);
            Prms = bit_width.extractParams(testCase);

            fxp_fr_fixed = struct('WL', P.WL, 'FL', P.FL_fixed);
            T_fr = freq_recovery.fxp_types(fxp_fr_fixed);

            % ---- Phase 1: Serial MEX builds into per-FL temp dirs -------
            fprintf('Building FR MEX (FL = %d, fixed)...\n', P.FL_fixed);
            frDir = bit_width.buildFRMex(P, fxp_fr_fixed);

            crDirs = cell(1, NFL);
            for bi = 1:NFL
                fl = P.FL_vec(bi);
                fprintf('[CR FL = %2d] building CR MEX  (%d / %d)\n', fl, bi, NFL);
                crDirs{bi} = bit_width.buildCRMex(P, struct('WL', P.WL, 'FL', fl));
            end

            % ---- Phase 2: Parallel SNR sweeps ---------------------------
            fecSNR_fft_DA    = nan(NFL, 2);
            fecSNR_dk_DA     = nan(NFL, 2);
            fecSNR_fft_blind = nan(NFL, 2);

            FL_vec_b   = P.FL_vec;
            WL_b       = P.WL;
            SNR_dB_b   = P.SNR_dB_vec;
            FEC_BER_b  = P.FEC_BER;
            BlindD_b   = P.BlindD;

            parfor pi = 1:NFL
                addpath(crDirs{pi});  %#ok<PFBNS>
                addpath(frDir);

                fl     = FL_vec_b(pi);
                T_cr_w = carrier_recovery.fxp_types(struct('WL', WL_b, 'FL', fl));

                ber = bit_width.runSnrSweepStatic(Prms, 'fft_search', 0, T_fr, T_cr_w);
                fecSNR_fft_DA(pi, :) = [ ...
                    bit_width.fecCrossing(SNR_dB_b, ber(:,1), FEC_BER_b), ...
                    bit_width.fecCrossing(SNR_dB_b, ber(:,2), FEC_BER_b)];

                ber = bit_width.runSnrSweepStatic(Prms, 'differential_kay', 0, T_fr, T_cr_w);
                fecSNR_dk_DA(pi, :) = [ ...
                    bit_width.fecCrossing(SNR_dB_b, ber(:,1), FEC_BER_b), ...
                    bit_width.fecCrossing(SNR_dB_b, ber(:,2), FEC_BER_b)];

                ber = bit_width.runSnrSweepStatic(Prms, 'fft_search_blind', BlindD_b, T_fr, T_cr_w);
                fecSNR_fft_blind(pi, :) = [ ...
                    bit_width.fecCrossing(SNR_dB_b, ber(:,1), FEC_BER_b), ...
                    bit_width.fecCrossing(SNR_dB_b, ber(:,2), FEC_BER_b)];
            end

            bit_width.plotCRSweep(P, P.FL_vec, fecSNR_fft_DA, fecSNR_dk_DA, fecSNR_fft_blind);

            cellfun(@(d) rmdir(d, 's'), crDirs, 'UniformOutput', false);
            rmdir(frDir, 's');
        end

    end

    %% ================================================================
    %  Static helpers — public interface for build/sweep/plot
    %% ================================================================
    methods (Static, Access = private)

        function tempDir = buildFRMex(P, fxp_fr)
            % Clears all loaded MEX, builds FR MEX to the default src/
            % location, then copies the result into a fresh temp directory
            % that mirrors the +freq_recovery package structure.  Returns
            % the temp directory path so callers can addpath it.
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

        function tempDir = buildCRMex(P, fxp_cr)
            % Same pattern as buildFRMex but for carrier recovery.
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

            srcDir  = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src');
            tempDir = tempname;
            dstPkg  = fullfile(tempDir, '+carrier_recovery');
            mkdir(dstPkg);
            ext = mexext;
            for f = {'viterbiViterbi_fxp_mex', 'pilots_only_fxp_mex'}
                copyfile( ...
                    fullfile(srcDir, '+carrier_recovery', [f{1} '.' ext]), ...
                    fullfile(dstPkg,                      [f{1} '.' ext]));
            end
        end

        function Params = extractParams(testCase)
            % Convert TestCase properties to a plain struct so it can be
            % broadcast into a parfor body without serialising the full
            % TestCase object.
            Params.Rs            = testCase.Rs;
            Params.N_pol         = testCase.N_pol;
            Params.NTrials       = testCase.NTrials;
            Params.SNR_dB_vec    = testCase.SNR_dB_vec;
            Params.DeltaF_Hz     = testCase.DeltaF_Hz;
            Params.LW_Hz         = testCase.LW_Hz;
            Params.FR_Nfft       = testCase.FR_Nfft;
            Params.FR_Po2Twiddle = testCase.FR_Po2Twiddle;
            Params.MaxFreq       = testCase.MaxFreq;
            Params.BlindD        = testCase.BlindD;
            Params.CordicIts     = testCase.CordicIts;
            Params.BlockLen      = testCase.BlockLen;
            Params.StepSize      = testCase.StepSize;
            Params.PilotThreshold = testCase.PilotThreshold;
            Params.VV_NTaps      = testCase.VV_NTaps;
            Params.FEC_BER       = testCase.FEC_BER;
            Params.VVFilters     = testCase.VVFilters;
        end

        function ber_avg = runSnrSweepStatic(Params, fr_algo, blindD, T_fr, T_cr)
            % Static version used inside parfor.  Accepts a plain Params
            % struct (from extractParams) rather than a TestCase instance.
            NSNR    = length(Params.SNR_dB_vec);
            ber_all = zeros(Params.NTrials, NSNR, 2);

            for tr = 1:Params.NTrials
                for si = 1:NSNR
                    [fr_out, pilots, txRefBits] = bit_width.buildChannel( ...
                        Params, Params.SNR_dB_vec(si), fr_algo, blindD, T_fr);

                    fr_fi     = cast(fr_out,                    'like', T_cr.x);
                    pilots_fi = cast(pilots,                     'like', T_cr.x);
                    vvfilt_fi = cast(Params.VVFilters{si},       'like', T_cr.w);

                    [cr_vv, ~] = carrier_recovery.viterbiViterbi_fxp_mex( ...
                        fr_fi, Params.N_pol, Params.VV_NTaps, vvfilt_fi, pilots_fi, ...
                        Params.BlockLen, double(Params.StepSize), Params.PilotThreshold, ...
                        double(Params.CordicIts), T_cr);
                    cr_vv = bit_width.resolveAmbiguity(double(cr_vv), txRefBits);
                    ber_all(tr, si, 1) = bit_width.computeBER(cr_vv, txRefBits);

                    [cr_po, ~] = carrier_recovery.pilots_only_fxp_mex( ...
                        fr_fi, Params.N_pol, Params.BlockLen, pilots_fi, ...
                        double(Params.CordicIts), T_cr);
                    cr_po = bit_width.resolveAmbiguity(double(cr_po), txRefBits);
                    ber_all(tr, si, 2) = bit_width.computeBER(cr_po, txRefBits);
                end
            end

            ber_avg = reshape(mean(ber_all, 1), [NSNR, 2]);

            nBits    = 3712 * Params.N_pol * 2;
            berFloor = 1 / (Params.NTrials * nBits);
            ber_avg(ber_avg == 0) = berFloor;
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
                    error('bit_width:unknownFR', 'Unknown FR algorithm: %s', fr_algo);
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

        function plotFRSweep(P, fl_vec, fecSNR_fft_DA, fecSNR_dk_DA, fecSNR_fft_blind)
            % 2 figures, one per CR method (VV / PO).
            % Each figure: 3 lines, one per FR variant.
            frLabels = {'FFT Search (data-aided)', 'Diff. Kay (data-aided)', ...
                        sprintf('FFT Search (blind D=%d)', P.BlindD)};
            crTitles = {'Viterbi-Viterbi', 'Pilots-only'};
            fecAll   = {fecSNR_fft_DA, fecSNR_dk_DA, fecSNR_fft_blind};
            colors   = lines(3);
            markers  = {'o', 's', '^'};

            for c = 1:2
                figure('Name', sprintf('FR bit width  |  %s', crTitles{c}), ...
                    'Position', [100 + (c-1)*80, 100, 750, 500], 'Color', 'w');
                ax = axes;
                hold(ax, 'on');
                grid(ax, 'on');

                for fr = 1:3
                    valid = ~isnan(fecAll{fr}(:, c));
                    if any(valid)
                        plot(ax, fl_vec(valid), fecAll{fr}(valid, c), ...
                            'LineStyle', '-', 'Marker', markers{fr}, ...
                            'MarkerSize', 6, 'LineWidth', 1.8, ...
                            'Color', colors(fr, :), ...
                            'DisplayName', frLabels{fr});
                    end
                end

                set(ax, 'FontSize', 11, 'Box', 'on', ...
                    'XTick', fl_vec, 'XLim', [fl_vec(1) - 1, fl_vec(end) + 1]);
                xlabel(ax, 'FR fractional bit width  [bits]', 'FontSize', 12);
                ylabel(ax, sprintf('SNR at BER = %.0e  [dB]', P.FEC_BER), 'FontSize', 12);
                title(ax, sprintf('%s  |  CR FL = %d fixed\nWL = %d,  \\DeltaF = %.0f MHz,  LW = %.0f kHz', ...
                    crTitles{c}, P.FL_fixed, P.WL, P.DeltaF_Hz/1e6, P.LW_Hz/1e3), 'FontSize', 11);
                legend(ax, 'Location', 'northeast', 'FontSize', 10);
            end
        end

        function plotCRSweep(P, fl_vec, fecSNR_fft_DA, fecSNR_dk_DA, fecSNR_fft_blind)
            % 3 figures, one per FR method.
            % Each figure: 2 lines, one per CR method (VV / PO).
            frTitles = {'FFT Search (data-aided)', 'Differential Kay (data-aided)', ...
                        sprintf('FFT Search (blind D=%d)', P.BlindD)};
            crLabels = {'Viterbi-Viterbi', 'Pilots-only'};
            fecAll   = {fecSNR_fft_DA, fecSNR_dk_DA, fecSNR_fft_blind};
            colors   = lines(2);
            markers  = {'o', 's'};

            for fr = 1:3
                figure('Name', sprintf('CR bit width  |  %s', frTitles{fr}), ...
                    'Position', [100 + (fr-1)*80, 100, 750, 500], 'Color', 'w');
                ax = axes;
                hold(ax, 'on');
                grid(ax, 'on');

                for c = 1:2
                    valid = ~isnan(fecAll{fr}(:, c));
                    if any(valid)
                        plot(ax, fl_vec(valid), fecAll{fr}(valid, c), ...
                            'LineStyle', '-', 'Marker', markers{c}, ...
                            'MarkerSize', 6, 'LineWidth', 1.8, ...
                            'Color', colors(c, :), ...
                            'DisplayName', crLabels{c});
                    end
                end

                set(ax, 'FontSize', 11, 'Box', 'on', ...
                    'XTick', fl_vec, 'XLim', [fl_vec(1) - 1, fl_vec(end) + 1]);
                xlabel(ax, 'CR fractional bit width  [bits]', 'FontSize', 12);
                ylabel(ax, sprintf('SNR at BER = %.0e  [dB]', P.FEC_BER), 'FontSize', 12);
                title(ax, sprintf('%s  |  FR FL = %d fixed\nWL = %d,  \\DeltaF = %.0f MHz,  LW = %.0f kHz', ...
                    frTitles{fr}, P.FL_fixed, P.WL, P.DeltaF_Hz/1e6, P.LW_Hz/1e3), 'FontSize', 11);
                legend(ax, 'Location', 'northeast', 'FontSize', 10);
            end
        end

    end
end
