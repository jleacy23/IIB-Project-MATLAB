classdef full_cr < matlab.unittest.TestCase
%FULL_CR  End-to-end carrier-recovery pipeline benchmark.
%
%   Tests every combination of:
%     Frequency Recovery : fft_search, fft_search_blind, differential_kay_data_aided, differential_kay_blind
%     Phase Recovery     : Viterbi-Viterbi, pilots_only
%   across a grid of SNRs, frequency offsets and laser linewidths.
%
%   Teardown outputs:
%     1) text summary table (absolute SNR at FEC);
%
%   Uses the true CPON symbol rate (30.5 GBd), frequency offsets in Hz,
%   and linewidths in Hz.
%
%   Run with:
%       runtests('full_cr')

    %% ================================================================
    %  Parameters
    %% ================================================================
    properties (Constant)

        % System
        Rs          = 30.5              % symbol rate [GBd]
        N_pol       = 2
        TrainingLen = 11                % training symbols per subframe
        NTrials     = 100               % independent channel realisations per point

        % Sweep grids
        SNR_dB_vec    = 0:0.5:20   % [dB]
        DeltaF_Hz_vec = [3e9]               % frequency offset [Hz]
        LW_Hz_vec     = [1000e3]            % laser linewidth  [Hz]

        % Frequency recovery — fixed-point settings
        FxpConfig_FR  = 'fixed32'       % 'fixed16' | 'fixed32'
        FR_Nfft       = 512             % FFT size for fft_search_fxp (power of 2)
        FR_Po2Twiddle = false           % round FFT twiddles to powers of 2
        FR_BlindD_vec = [64, 256, 512] % blind data lengths [symbols] (FFT + DiffKay)
        MaxFreq = 0.1

        % Phase recovery — shared settings
        BlockLen       = 32             % CPON block length [symbols]
        StepSize       = 32
        PilotThreshold = 5 * pi / 9           % cycle-slip detection threshold [rad]

        % Viterbi-Viterbi
        VV_NTaps = 10

        % Enable/disable figure output
        Plot = true

        % Fixed-point configuration (CR)
        FxpConfig = 'fixed16'           % 'fixed16' | 'fixed32'
        CordicIts = 16                  % CORDIC iterations (shared FR + CR)

        % Enable/disable MEX rebuild
        Rebuild = true

    end

    %% ================================================================
    %  Pre-computed resources
    %% ================================================================
    properties
        VVFilters   % {NSNR x NLW} cell of Wiener VV filters
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
            full_cr.fecSummaryStore('reset');
            full_cr.berSummaryStore('reset');
        end

        function buildMex(testCase)
            buildDir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'build');
            addpath(buildDir);

            P.Rs             = testCase.Rs;
            P.N_pol          = testCase.N_pol;
            P.TrainingLen    = testCase.TrainingLen;
            P.FR_Nfft        = testCase.FR_Nfft;
            P.FR_Po2Twiddle  = testCase.FR_Po2Twiddle;
            P.FR_BlindD      = testCase.FR_BlindD_vec(1);
            P.FxpConfig_FR   = testCase.FxpConfig_FR;
            P.FxpConfig_VV   = testCase.FxpConfig;
            P.FxpConfig_PO   = testCase.FxpConfig;
            P.CordicIts      = testCase.CordicIts;
            P.VV_NTaps       = testCase.VV_NTaps;
            P.BlockLen       = testCase.BlockLen;
            P.StepSize       = testCase.StepSize;
            P.PilotThreshold = testCase.PilotThreshold;
            P.PilotLen       = 1;
            P.MaxFreq        = testCase.MaxFreq;

            cfg = coder.config('mex');
            cfg.GenerateReport = false;

            if testCase.Rebuild
                fprintf('Building MEX objects...\n');
                build_freq_recovery_fft_search_fxp_mex(P, cfg);
                build_freq_recovery_differential_kay_fxp_mex(P, cfg);
                build_carrier_recovery_viterbiViterbi_fxp_mex(P, cfg);
                build_carrier_recovery_pilots_only_fxp_mex(P, cfg);
                fprintf('All MEX objects built.\n');
            end
        end

        function precomputeVVFilters(testCase)
            % Estimate symbol energy from one real CPON subframe so the
            % Wiener filter is calibrated to the true mixed pilot/data energy.
            BITS_PER_SF = 3586 * 2 * 2;
            [tmp_syms, ~, ~, ~] = modem.modulate(modem.randomBits(BITS_PER_SF));
            symEnergy = mean(abs(tmp_syms(:)).^2);

            NSNR = length(testCase.SNR_dB_vec);
            NLW  = length(testCase.LW_Hz_vec);
            testCase.VVFilters = cell(NSNR, NLW);

            for si = 1:NSNR
                for li = 1:NLW
                    testCase.VVFilters{si, li} = carrier_recovery.genVVFilter( ...
                        testCase.LW_Hz_vec(li), testCase.Rs,       ...
                        testCase.SNR_dB_vec(si), symEnergy,        ...
                        testCase.N_pol, testCase.VV_NTaps);
                end
            end
        end

    end

    %% ================================================================
    %  Test class teardown
    %% ================================================================
    methods (TestClassTeardown)

        function printFecSummaryAtEnd(testCase)
            [absSNR, rowLabels, colLabels] = full_cr.buildAbsoluteSnrTable(testCase);
            full_cr.writeAbsoluteSnrSummaryFile(testCase, absSNR, rowLabels, colLabels);
            full_cr.fecSummaryStore('reset');
            full_cr.berSummaryStore('reset');
        end

    end

    %% ================================================================
    %  Tests — one per frequency-recovery algorithm
    %% ================================================================
    methods (Test)

        function test_fft_search(testCase)
            testCase.runPipeline('fft_search', 0);
        end

        function test_fft_search_blind(testCase)
            for d = testCase.FR_BlindD_vec
                testCase.runPipeline('fft_search_blind', d);
            end
        end

        function test_differential_kay_blind(testCase)
            for d = testCase.FR_BlindD_vec
                testCase.runPipeline('differential_kay_blind', d);
            end
        end

        function test_differential_kay_data_aided(testCase)
            testCase.runPipeline('differential_kay_data_aided', 0);
        end

    end

    %% ================================================================
    %  Private pipeline implementation
    %% ================================================================
    methods (Access = private)

        function runPipeline(testCase, fr_algo, blindD)
            P    = testCase;
            NSNR = length(P.SNR_dB_vec);
            NFO  = length(P.DeltaF_Hz_vec);
            NLW  = length(P.LW_Hz_vec);

            if nargin < 3
                blindD = 0;
            end

            % BER storage: (trial, SNR, DeltaF, LW, PR)
            %   PR index:  1=ViterbiViterbi  2=PilotsOnly
            BER_all = zeros(P.NTrials, NSNR, NFO, NLW, 2);

            for tr = 1:P.NTrials
                fprintf('[%s] Trial %d / %d\n', fr_algo, tr, P.NTrials);

                for fi = 1:NFO
                    DeltaF_Hz = P.DeltaF_Hz_vec(fi);

                    for li = 1:NLW
                        LW = P.LW_Hz_vec(li);

                        T_cr = carrier_recovery.fxp_types(P.FxpConfig);

                        for si = 1:NSNR
                            SNR_dB   = P.SNR_dB_vec(si);
                            VVFilter = P.VVFilters{si, li};

                            [fr_out, pilots, txRefBits, ~] = ...
                                full_cr.buildChannel(P, SNR_dB, DeltaF_Hz, LW, fr_algo, blindD);

                            % fprintf('    SNR=%2ddB  freq_offset_est = %+.3f MHz\n', ...
                            %     SNR_dB, freq_offset/1e6);

                            fr_out_fi   = cast(fr_out,   'like', T_cr.x);
                            pilots_fi   = cast(pilots,   'like', T_cr.x);
                            vvfilter_fi = cast(VVFilter, 'like', T_cr.w);

                            %-- Viterbi-Viterbi (fxp MEX) --
                            [cr_vv, ~] = carrier_recovery.viterbiViterbi_fxp_mex( ...
                                fr_out_fi, P.N_pol, P.VV_NTaps, vvfilter_fi, ...
                                pilots_fi, P.BlockLen, double(P.StepSize), P.PilotThreshold, ...
                                double(P.CordicIts), T_cr);
                            cr_vv = full_cr.resolvePhaseAmbiguity(double(cr_vv), txRefBits);
                            BER_all(tr, si, fi, li, 1) = full_cr.computeBER(cr_vv, txRefBits);

                            %-- Pilots only (fxp MEX) --
                            [cr_po, ~] = carrier_recovery.pilots_only_fxp_mex( ...
                                fr_out_fi, P.N_pol, P.BlockLen, pilots_fi, double(P.CordicIts), T_cr);
                            cr_po = full_cr.resolvePhaseAmbiguity(double(cr_po), txRefBits);
                            BER_all(tr, si, fi, li, 2) = full_cr.computeBER(cr_po, txRefBits);

                        end % SNR
                    end % LW
                end % DeltaF
            end % trial

            % Average over trials → [NSNR x NFO x NLW x 2]
            BER = reshape(mean(BER_all, 1), [NSNR, NFO, NLW, 2]);

            % BER = 0 cannot be plotted on a log scale.  Replace with the
            % minimum observable BER given NTrials * bitsPerTrial total bits.
            % CPON subframe: 3712 symbols x N_pol x 2 bits/QPSK symbol.
            nBits    = 3712 * P.N_pol * 2;
            berFloor = 1 / (P.NTrials * nBits);
            BER(BER == 0) = berFloor;

            fecSNR = full_cr.computeFecCrossingSNR(P, BER, 2e-2);
            full_cr.fecSummaryStore('set', fr_algo, blindD, fecSNR);
            full_cr.berSummaryStore('set', fr_algo, blindD, BER, berFloor);

            full_cr.printResults(P, BER, fr_algo);
        end

    end

    %% ================================================================
    %  Static helpers
    %% ================================================================
    methods (Static, Access = private)

        function storeOut = fecSummaryStore(action, fr_algo, blindD, fecSNR)
            persistent S
            if isempty(S)
                S = struct();
            end

            if nargin < 2
                fr_algo = '';
            end
            if nargin < 3
                blindD = 0;
            end
            if nargin < 4
                fecSNR = [];
            end

            switch action
                case 'reset'
                    S = struct();
                case 'set'
                    key = full_cr.makeScenarioKey(fr_algo, blindD);
                    S.(key) = fecSNR;
                case 'get'
                    % no-op
                otherwise
                    error('full_cr:invalidSummaryAction', 'Unknown action: %s', action);
            end
            storeOut = S;
        end

        function key = makeScenarioKey(fr_algo, blindD)
            key = sprintf('%s_D%d', fr_algo, round(blindD));
            key = strrep(key, '+', 'p');
            key = strrep(key, '-', 'm');
            key = matlab.lang.makeValidName(key);
        end

        function storeOut = berSummaryStore(action, fr_algo, blindD, BER, berFloor)
            persistent S
            if isempty(S)
                S = struct();
            end

            if nargin < 2
                fr_algo = '';
            end
            if nargin < 3
                blindD = 0;
            end
            if nargin < 4
                BER = [];
            end
            if nargin < 5
                berFloor = NaN;
            end

            switch action
                case 'reset'
                    S = struct();
                case 'set'
                    key = full_cr.makeScenarioKey(fr_algo, blindD);
                    S.(key).BER = BER;
                    S.(key).berFloor = berFloor;
                case 'get'
                    % no-op
                otherwise
                    error('full_cr:invalidBerSummaryAction', 'Unknown action: %s', action);
            end
            storeOut = S;
        end

        function fecSNR = computeFecCrossingSNR(P, BER, fecLimit)
            % BER: [NSNR x NFO x NLW x 2]
            NSNR = length(P.SNR_dB_vec);
            NFO  = length(P.DeltaF_Hz_vec);
            NLW  = length(P.LW_Hz_vec);
            fecSNR = nan(NFO, NLW, 2);

            x = P.SNR_dB_vec(:);
            for fi = 1:NFO
                for li = 1:NLW
                    for pr = 1:2
                        y = squeeze(BER(:, fi, li, pr));
                        if length(y) ~= NSNR
                            continue;
                        end
                        fecSNR(fi, li, pr) = full_cr.interpolateFecCrossing(x, y, fecLimit);
                    end
                end
            end
        end

        function xCross = interpolateFecCrossing(x, y, yLimit)
            % Returns first SNR crossing (in ascending SNR order) where BER
            % reaches yLimit using linear interpolation between bracketing points.
            xCross = NaN;
            n = length(x);
            if n < 2
                return;
            end

            % Exact hit takes precedence.
            idxExact = find(y == yLimit, 1, 'first');
            if ~isempty(idxExact)
                xCross = x(idxExact);
                return;
            end

            for i = 1:(n - 1)
                y1 = y(i);
                y2 = y(i + 1);
                if (y1 - yLimit) * (y2 - yLimit) < 0
                    x1 = x(i);
                    x2 = x(i + 1);
                    xCross = x1 + (yLimit - y1) * (x2 - x1) / (y2 - y1);
                    return;
                end
            end
        end

        function [absSNR, rowLabels, colLabels] = buildAbsoluteSnrTable(P)
            S = full_cr.fecSummaryStore('get');
            NFO = length(P.DeltaF_Hz_vec);
            NLW = length(P.LW_Hz_vec);

            colLabels = {'FFT|VV', 'FFT|PO', 'DK|VV', 'DK|PO'};
            rowLabels = cell(1, 1 + length(P.FR_BlindD_vec));
            rowLabels{1} = 'DA';
            for r = 2:length(rowLabels)
                rowLabels{r} = sprintf('BL-D%d', P.FR_BlindD_vec(r - 1));
            end

            nRows = length(rowLabels);
            nCols = length(colLabels);
            absSNR = nan(nRows, nCols, NFO, NLW);

            key_fft_da = full_cr.makeScenarioKey('fft_search', 0);
            key_dk_da  = full_cr.makeScenarioKey('differential_kay_data_aided', 0);
            if isfield(S, key_fft_da)
                absSNR(1, 1, :, :) = S.(key_fft_da)(:, :, 1);
                absSNR(1, 2, :, :) = S.(key_fft_da)(:, :, 2);
            end
            if isfield(S, key_dk_da)
                absSNR(1, 3, :, :) = S.(key_dk_da)(:, :, 1);
                absSNR(1, 4, :, :) = S.(key_dk_da)(:, :, 2);
            end

            for di = 1:length(P.FR_BlindD_vec)
                D = P.FR_BlindD_vec(di);
                r = di + 1;
                key_fft_bl = full_cr.makeScenarioKey('fft_search_blind', D);
                key_dk_bl  = full_cr.makeScenarioKey('differential_kay_blind', D);
                if isfield(S, key_fft_bl)
                    absSNR(r, 1, :, :) = S.(key_fft_bl)(:, :, 1);
                    absSNR(r, 2, :, :) = S.(key_fft_bl)(:, :, 2);
                end
                if isfield(S, key_dk_bl)
                    absSNR(r, 3, :, :) = S.(key_dk_bl)(:, :, 1);
                    absSNR(r, 4, :, :) = S.(key_dk_bl)(:, :, 2);
                end
            end

        end

        function writeAbsoluteSnrSummaryFile(P, absSNR, rowLabels, colLabels)
            NFO = length(P.DeltaF_Hz_vec);
            NLW = length(P.LW_Hz_vec);
            outPath = fullfile(fileparts(mfilename('fullpath')), 'fec_absolute_snr_summary.txt');
            fid = fopen(outPath, 'w');
            if fid < 0
                error('full_cr:summaryWriteFailed', 'Failed to open summary file: %s', outPath);
            end

            fprintf(fid, 'Absolute SNR at FEC summary\n');
            fprintf(fid, 'FEC limit = 2e-2\n\n');

            for fi = 1:NFO
                for li = 1:NLW
                    fprintf(fid, 'DeltaF = %.0f MHz, LW = %.0f kHz\n', ...
                        P.DeltaF_Hz_vec(fi)/1e6, P.LW_Hz_vec(li)/1e3);
                    fprintf(fid, '%-10s', 'Row');
                    for c = 1:length(colLabels)
                        fprintf(fid, '  %-12s', colLabels{c});
                    end
                    fprintf(fid, '\n');

                    for r = 1:length(rowLabels)
                        fprintf(fid, '%-10s', rowLabels{r});
                        for c = 1:length(colLabels)
                            val = absSNR(r, c, fi, li);
                            if isnan(val)
                                fprintf(fid, '  %-12s', 'N/A');
                            else
                                fprintf(fid, '  %-12.3f', val);
                            end
                        end
                        fprintf(fid, '\n');
                    end
                    fprintf(fid, '\n');
                end
            end
            fclose(fid);
            fprintf('Wrote summary table to: %s\n', outPath);
        end

        function printDataAidedBerSummary(P, colLabels)
            S = full_cr.berSummaryStore('get');
            key_fft_da = full_cr.makeScenarioKey('fft_search', 0);
            key_dk_da  = full_cr.makeScenarioKey('differential_kay_data_aided', 0);
            if ~isfield(S, key_fft_da) || ~isfield(S, key_dk_da)
                fprintf('Data-aided BER summary skipped: missing DA BER results.\n');
                return;
            end

            ber_fft = S.(key_fft_da).BER;
            ber_dk  = S.(key_dk_da).BER;

            if nargin < 2 || isempty(colLabels)
                colLabels = {'FFT|VV', 'FFT|PO', 'DK|VV', 'DK|PO'};
            end

            fprintf('\n==============================================================\n');
            fprintf('Data-aided BER vs SNR summary\n');
            fprintf('Columns: %s | %s | %s | %s\n', ...
                colLabels{1}, colLabels{2}, colLabels{3}, colLabels{4});
            fprintf('==============================================================\n');

            for fi = 1:length(P.DeltaF_Hz_vec)
                for li = 1:length(P.LW_Hz_vec)
                    fprintf('\nDeltaF = %.0f MHz, LW = %.0f kHz\n', ...
                        P.DeltaF_Hz_vec(fi)/1e6, P.LW_Hz_vec(li)/1e3);
                    fprintf('%8s  %12s  %12s  %12s  %12s\n', ...
                        'SNR[dB]', colLabels{1}, colLabels{2}, colLabels{3}, colLabels{4});

                    for si = 1:length(P.SNR_dB_vec)
                        fprintf('%8.1f  %12.5e  %12.5e  %12.5e  %12.5e\n', ...
                            P.SNR_dB_vec(si), ...
                            ber_fft(si, fi, li, 1), ...
                            ber_fft(si, fi, li, 2), ...
                            ber_dk(si, fi, li, 1), ...
                            ber_dk(si, fi, li, 2));
                    end
                end
            end
        end

        function plotDataAidedBerCurves(P, colLabels)
            S = full_cr.berSummaryStore('get');
            key_fft_da = full_cr.makeScenarioKey('fft_search', 0);
            key_dk_da  = full_cr.makeScenarioKey('differential_kay_data_aided', 0);
            if ~isfield(S, key_fft_da) || ~isfield(S, key_dk_da)
                return;
            end

            if nargin < 2 || isempty(colLabels)
                colLabels = {'FFT|VV', 'FFT|PO', 'DK|VV', 'DK|PO'};
            end

            ber_fft = S.(key_fft_da).BER;
            ber_dk  = S.(key_dk_da).BER;
            berFloor = min(S.(key_fft_da).berFloor, S.(key_dk_da).berFloor);

            NFO = length(P.DeltaF_Hz_vec);
            NLW = length(P.LW_Hz_vec);
            curveColors = lines(4);
            curveMarkers = {'o', 's', 'd', '^'};

            fig_w = max(900, 420 * NLW);
            fig_h = max(600, 360 * NFO);
            figure('Name', 'Full CR  |  BER vs SNR (Data-Aided Combinations)', ...
                'Position', [80, 80, fig_w, fig_h], 'Color', 'w');

            for fi = 1:NFO
                for li = 1:NLW
                    ax = subplot(NFO, NLW, (fi - 1) * NLW + li);
                    set(ax, 'YScale', 'log', 'FontSize', 11, 'Box', 'on', 'Color', 'w');
                    hold(ax, 'on');

                    yData = [ ...
                        ber_fft(:, fi, li, 1); ...
                        ber_fft(:, fi, li, 2); ...
                        ber_dk(:, fi, li, 1); ...
                        ber_dk(:, fi, li, 2)  ...
                    ];

                    semilogy(ax, P.SNR_dB_vec, ber_fft(:, fi, li, 1), ...
                        'LineStyle', '-', 'Marker', curveMarkers{1}, 'MarkerSize', 4, 'LineWidth', 1.8, ...
                        'Color', curveColors(1, :), 'DisplayName', colLabels{1});
                    semilogy(ax, P.SNR_dB_vec, ber_fft(:, fi, li, 2), ...
                        'LineStyle', '-', 'Marker', curveMarkers{2}, 'MarkerSize', 4, 'LineWidth', 1.8, ...
                        'Color', curveColors(2, :), 'DisplayName', colLabels{2});
                    semilogy(ax, P.SNR_dB_vec, ber_dk(:, fi, li, 1), ...
                        'LineStyle', '-', 'Marker', curveMarkers{3}, 'MarkerSize', 4, 'LineWidth', 1.8, ...
                        'Color', curveColors(3, :), 'DisplayName', colLabels{3});
                    semilogy(ax, P.SNR_dB_vec, ber_dk(:, fi, li, 2), ...
                        'LineStyle', '-', 'Marker', curveMarkers{4}, 'MarkerSize', 4, 'LineWidth', 1.8, ...
                        'Color', curveColors(4, :), 'DisplayName', colLabels{4});

                    yline(ax, 2e-2, 'k--', 'LineWidth', 1.2, 'DisplayName', 'FEC limit (2\times10^{-2})');
                    yline(ax, berFloor, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1.0, ...
                        'DisplayName', sprintf('Zero-error floor (%.2g)', berFloor));

                    full_cr.setBerYAxisToData(ax, yData, berFloor);

                    grid(ax, 'on');
                    xlabel(ax, 'SNR [dB]', 'FontSize', 11);
                    ylabel(ax, 'BER', 'FontSize', 11);
                    title(ax, sprintf('\\DeltaF = %.0f MHz, LW = %.0f kHz', ...
                        P.DeltaF_Hz_vec(fi)/1e6, P.LW_Hz_vec(li)/1e3), 'FontSize', 11);
                    legend(ax, 'Location', 'southwest', 'FontSize', 8);
                end
            end

            sgtitle('BER vs SNR  |  Data-Aided FR+PR Combinations', ...
                'FontSize', 14, 'FontWeight', 'bold');
        end

        function plotBerColumnsDaVsBlind(P, colLabels)
            S = full_cr.berSummaryStore('get');
            key_fft_da = full_cr.makeScenarioKey('fft_search', 0);
            key_dk_da  = full_cr.makeScenarioKey('differential_kay_data_aided', 0);
            if ~isfield(S, key_fft_da) || ~isfield(S, key_dk_da)
                return;
            end

            berFftDa = S.(key_fft_da).BER;
            berDkDa  = S.(key_dk_da).BER;
            berFloor = min(S.(key_fft_da).berFloor, S.(key_dk_da).berFloor);

            nBlind = length(P.FR_BlindD_vec);
            berFftBlind = cell(1, nBlind);
            berDkBlind  = cell(1, nBlind);
            for di = 1:nBlind
                D = P.FR_BlindD_vec(di);
                key_fft_bl = full_cr.makeScenarioKey('fft_search_blind', D);
                key_dk_bl  = full_cr.makeScenarioKey('differential_kay_blind', D);
                if isfield(S, key_fft_bl)
                    berFftBlind{di} = S.(key_fft_bl).BER;
                    berFloor = min(berFloor, S.(key_fft_bl).berFloor);
                end
                if isfield(S, key_dk_bl)
                    berDkBlind{di} = S.(key_dk_bl).BER;
                    berFloor = min(berFloor, S.(key_dk_bl).berFloor);
                end
            end

            NFO = length(P.DeltaF_Hz_vec);
            NLW = length(P.LW_Hz_vec);

            comboNames = {'FFT|VV', 'FFT|PO', 'DK|VV', 'DK|PO'};
            lineColors = lines(1 + nBlind);      % DA + blind-D curves
            blindMarkers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', 'x', '+'};

            for ci = 1:length(comboNames)
                fig_w = max(900, 420 * NLW);
                fig_h = max(600, 360 * NFO);
                figure('Name', sprintf('Full CR  |  BER vs SNR (%s)', comboNames{ci}), ...
                    'Position', [90, 90, fig_w, fig_h], ...
                    'Color', 'w');

                for fi = 1:NFO
                    for li = 1:NLW
                        ax = subplot(NFO, NLW, (fi - 1) * NLW + li);
                        set(ax, 'YScale', 'log', 'FontSize', 11, 'Box', 'on', 'Color', 'w');
                        hold(ax, 'on');

                        [berDa, berBlindSet] = full_cr.getBerForColumn(ci, berFftDa, berDkDa, berFftBlind, berDkBlind);
                        yData = berDa(:, fi, li);

                        semilogy(ax, P.SNR_dB_vec, berDa(:, fi, li), ...
                            'LineStyle', '-', 'Marker', 'o', 'MarkerSize', 4, 'LineWidth', 1.8, ...
                            'Color', lineColors(1, :), 'DisplayName', 'DA');

                        for di = 1:nBlind
                            if isempty(berBlindSet{di})
                                continue;
                            end
                            markerIdx = mod(di - 1, length(blindMarkers)) + 1;
                            semilogy(ax, P.SNR_dB_vec, berBlindSet{di}(:, fi, li), ...
                                'LineStyle', '--', 'Marker', blindMarkers{markerIdx}, ...
                                'MarkerSize', 3.5, 'LineWidth', 1.6, ...
                                'Color', lineColors(di + 1, :), ...
                                'DisplayName', sprintf('BL-D%d', P.FR_BlindD_vec(di)));
                            yData = [yData; berBlindSet{di}(:, fi, li)];
                        end

                        yline(ax, 2e-2, 'k--', 'LineWidth', 1.2, ...
                            'DisplayName', 'FEC limit (2\times10^{-2})');
                        yline(ax, berFloor, 'Color', [0.5 0.5 0.5], ...
                            'LineStyle', '--', 'LineWidth', 1.0, ...
                            'DisplayName', sprintf('Zero-error floor (%.2g)', berFloor));

                        full_cr.setBerYAxisToData(ax, yData, berFloor);

                        grid(ax, 'on');
                        xlabel(ax, 'SNR [dB]', 'FontSize', 11);
                        ylabel(ax, 'BER', 'FontSize', 11);
                        title(ax, sprintf('\\DeltaF = %.0f MHz, LW = %.0f kHz', ...
                            P.DeltaF_Hz_vec(fi)/1e6, P.LW_Hz_vec(li)/1e3), ...
                            'FontSize', 11);
                        legend(ax, 'Location', 'southwest', 'FontSize', 8);
                    end
                end

                if nargin >= 2 && ~isempty(colLabels) && ci <= length(colLabels)
                    plotTitle = colLabels{ci};
                else
                    plotTitle = comboNames{ci};
                end
                sgtitle(sprintf('BER vs SNR  |  %s  (DA solid, blind dashed)', plotTitle), ...
                    'FontSize', 14, 'FontWeight', 'bold');
            end
        end

        function [berDa, berBlindSet] = getBerForColumn(colIdx, berFftDa, berDkDa, berFftBlind, berDkBlind)
            switch colIdx
                case 1
                    berDa = berFftDa(:, :, :, 1);
                    berBlindSet = cell(size(berFftBlind));
                    for i = 1:length(berFftBlind)
                        if ~isempty(berFftBlind{i})
                            berBlindSet{i} = berFftBlind{i}(:, :, :, 1);
                        else
                            berBlindSet{i} = [];
                        end
                    end
                case 2
                    berDa = berFftDa(:, :, :, 2);
                    berBlindSet = cell(size(berFftBlind));
                    for i = 1:length(berFftBlind)
                        if ~isempty(berFftBlind{i})
                            berBlindSet{i} = berFftBlind{i}(:, :, :, 2);
                        else
                            berBlindSet{i} = [];
                        end
                    end
                case 3
                    berDa = berDkDa(:, :, :, 1);
                    berBlindSet = cell(size(berDkBlind));
                    for i = 1:length(berDkBlind)
                        if ~isempty(berDkBlind{i})
                            berBlindSet{i} = berDkBlind{i}(:, :, :, 1);
                        else
                            berBlindSet{i} = [];
                        end
                    end
                case 4
                    berDa = berDkDa(:, :, :, 2);
                    berBlindSet = cell(size(berDkBlind));
                    for i = 1:length(berDkBlind)
                        if ~isempty(berDkBlind{i})
                            berBlindSet{i} = berDkBlind{i}(:, :, :, 2);
                        else
                            berBlindSet{i} = [];
                        end
                    end
                otherwise
                    error('full_cr:invalidColumn', 'Unknown column index: %d', colIdx);
            end
        end

        function setBerYAxisToData(ax, yData, berFloor)
            if nargin < 3
                berFloor = NaN;
            end

            y = yData(isfinite(yData) & yData > 0);
            if isempty(y)
                return;
            end

            % Ignore floor-clipped points when selecting limits so the axis
            % follows the visible BER trend instead of collapsing to berFloor.
            if isfinite(berFloor) && berFloor > 0
                yNoFloor = y(y > 1.05 * berFloor);
                if ~isempty(yNoFloor)
                    y = yNoFloor;
                end
            end

            yMin = min(y) / 1.5;
            yMax = max(y) * 1.5;
            yMin = max(yMin, 1e-8);
            yMax = min(yMax, 1);
            if yMax <= yMin
                yMax = yMin * 10;
            end

            ylim(ax, [yMin, yMax]);
        end

        function [fr_out, pilots, txRefBits, freq_offset] = buildChannel( ...
            P, SNR_dB, DeltaF_Hz, LW, fr_algo, blindD)
            % Generate one CPON subframe, apply channel impairments and
            % frequency recovery, then build the per-CR-block pilot matrix.
            CPON_BLOCK_LEN = 32;
            CPON_SF_SYMS   = 3712;
            N_CPON_BLOCKS  = 116;
            BITS_PER_SF    = 3586 * 2 * 2;   % data_syms x N_pol x bits/QPSK_sym

            txBits = modem.randomBits(BITS_PER_SF);
            [symbols, pilotSyms, training, ~] = modem.modulate(txBits);

            % Channel: frequency offset -> AWGN -> phase noise
            rx = channel.lo_freq_shift(symbols, DeltaF_Hz / 1e6, P.Rs, 1);
            rx = channel.add_awgn(rx, SNR_dB);
            rx = channel.add_phase_noise(rx, P.Rs, LW);

            % Frequency recovery (fixed-point MEX)
            T_fr   = freq_recovery.fxp_types(P.FxpConfig_FR);
            rx_fi  = cast(rx,       'like', T_fr.x);
            tr_fi  = cast(training, 'like', T_fr.x);
            switch fr_algo
                case 'fft_search'
                    [fr_out, freq_offset] = freq_recovery.fft_search_fxp_mex( ...
                        rx_fi, tr_fi, P.Rs, P.FR_Nfft, P.FR_Po2Twiddle, P.CordicIts, P.MaxFreq, T_fr, true, 0);
                case 'fft_search_blind'
                    [fr_out, freq_offset] = freq_recovery.fft_search_fxp_mex( ...
                        rx_fi, tr_fi, P.Rs, P.FR_Nfft, P.FR_Po2Twiddle, P.CordicIts, P.MaxFreq, T_fr, false, blindD);
                case 'differential_kay_data_aided'
                    [fr_out, freq_offset] = freq_recovery.differential_kay_fxp_mex( ...
                        rx_fi, tr_fi, P.Rs, P.CordicIts, T_fr, true, 0, P.MaxFreq);
                case {'differential_kay_blind', 'differential_kay'}
                    [fr_out, freq_offset] = freq_recovery.differential_kay_fxp_mex( ...
                        rx_fi, tr_fi, P.Rs, P.CordicIts, T_fr, false, blindD, P.MaxFreq);
                otherwise
                    error('full_cr:unknownFR', 'Unknown FR algorithm: %s', fr_algo);
            end
            fr_out = double(fr_out);

            % Build per-CR-block pilot matrix [NBlocks x N_pol].
            % For any BlockLen, CR block b starts at signal position
            % (b-1)*BlockLen+1; map through the CPON subframe layout to
            % find the corresponding CPON pilot symbol.
            Nsym    = size(fr_out, 1);
            NBlocks = ceil(Nsym / P.BlockLen);
            pilots  = zeros(NBlocks, P.N_pol);
            for b = 1:NBlocks
                pos       = (b - 1) * P.BlockLen + 1;
                posInSf   = mod(pos - 1, CPON_SF_SYMS) + 1;
                cponBlock = floor((posInSf - 1) / CPON_BLOCK_LEN) + 1;
                cponBlock = min(cponBlock, N_CPON_BLOCKS);
                pilots(b, :) = pilotSyms(cponBlock, :);
            end

            txRefBits = modem.symbolsToBits(symbols);
        end

        function BER = computeBER(crSym, refBits)
            decidedSyms = modem.decideSymbols(crSym);
            rxBits      = modem.symbolsToBits(decidedSyms);
            nBits       = min(length(refBits), length(rxBits));
            BER         = sum(refBits(1:nBits) ~= rxBits(1:nBits)) / nBits;
        end

        function bestSym = resolvePhaseAmbiguity(crSym, txRefBits)
            bestBER = Inf;
            bestSym = crSym;
            for k = 0:3
                rotated     = crSym .* exp(-1j * k * pi/2);
                decidedSyms = modem.decideSymbols(rotated);
                rxBits      = modem.symbolsToBits(decidedSyms);
                nBits       = min(length(txRefBits), length(rxBits));
                thisBER     = sum(txRefBits(1:nBits) ~= rxBits(1:nBits)) / nBits;
                if thisBER < bestBER
                    bestBER = thisBER;
                    bestSym = rotated;
                end
            end
        end

        function plotConstellation(sym, titleStr)
            N_pol = size(sym, 2);
            figure('Name', titleStr, 'Position', [100 100 900 400], 'Color', 'w');
            for p = 1:N_pol
                subplot(1, N_pol, p);
                plot(real(sym(:,p)), imag(sym(:,p)), '.', 'MarkerSize', 2);
                grid on; axis equal;
                title(sprintf('Pol %d', p), 'FontSize', 11);
                xlabel('In-Phase'); ylabel('Quadrature');
            end
            sgtitle(titleStr, 'FontSize', 12, 'Interpreter', 'none');
        end

        function printResults(P, BER, fr_algo)
            % BER: [NSNR x NFO x NLW x 2]
            pr_names = {'VV', 'PilotsOnly'};
            for fi = 1:length(P.DeltaF_Hz_vec)
                for li = 1:length(P.LW_Hz_vec)
                    fprintf('\n[%s]  DeltaF = %.0f MHz  |  LW = %.0f kHz\n', ...
                        fr_algo, P.DeltaF_Hz_vec(fi)/1e6, P.LW_Hz_vec(li)/1e3);
                    fprintf('  SNR [dB]    : ');
                    fprintf('%7.1f  ', P.SNR_dB_vec);
                    fprintf('\n');
                    for pr = 1:2
                        fprintf('  %-12s: ', pr_names{pr});
                        fprintf('%7.5f  ', BER(:, fi, li, pr)');
                        fprintf('\n');
                    end
                end
            end
        end

    end
end
