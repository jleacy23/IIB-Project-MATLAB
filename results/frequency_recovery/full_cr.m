classdef full_cr < matlab.unittest.TestCase
%FULL_CR  End-to-end carrier-recovery pipeline benchmark.
%
%   Tests every combination of:
%     Frequency Recovery : fft_search, differential_kay
%     Phase Recovery     : BPS, Viterbi-Viterbi, pilots_only
%   across a grid of SNRs, frequency offsets and laser linewidths.
%
%   One figure is produced per frequency-recovery algorithm.  Each figure
%   contains NFO x NLW subplots (frequency offsets as rows, linewidths as
%   columns); within each subplot BER vs SNR is plotted with one line per
%   phase-recovery algorithm.
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
        NTrials     = 100                % independent channel realisations per point

        % Sweep grids
        SNR_dB_vec    = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]   % [dB]
        DeltaF_Hz_vec = [3e9]               % frequency offset [Hz]
        LW_Hz_vec     = [1000e3]            % laser linewidth  [Hz]

        % FFT pruning factor
        FFT_K = 47

        % Phase recovery — shared settings
        BlockLen       = 32             % CPON block length [symbols]
        StepSize       = 1
        PilotThreshold = 5 * pi / 9           % cycle-slip detection threshold [rad]

        % Viterbi-Viterbi
        VV_NTaps = 5

        % BPS
        BPS_N = 5
        BPS_B = 64
        BPS_M = 4                       % QPSK

        % Enable/disable figure output
        Plot = true

        % Fixed-point configuration (CR only)
        FxpConfig = 'fixed16'           % 'fixed16' | 'fixed32'
        CordicIts = 16                  % CORDIC iterations

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
        end

        function buildMex(testCase)
            buildDir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'build');
            addpath(buildDir);

            P.Rs             = testCase.Rs;
            P.N_pol          = testCase.N_pol;
            P.FxpConfig_BPS  = testCase.FxpConfig;
            P.FxpConfig_VV   = testCase.FxpConfig;
            P.CordicIts      = testCase.CordicIts;
            P.BPS_N          = testCase.BPS_N;
            P.BPS_B          = testCase.BPS_B;
            P.M              = testCase.BPS_M;
            P.VV_NTaps       = testCase.VV_NTaps;
            P.BlockLen       = testCase.BlockLen;
            P.StepSize       = testCase.StepSize;
            P.PilotThreshold = testCase.PilotThreshold;
            P.PilotLen       = 1;

            cfg = coder.config('mex');
            cfg.GenerateReport = false;

            if testCase.Rebuild
                fprintf('Building MEX objects...\n');
                build_carrier_recovery_bps_fxp_mex(P, cfg);
                build_carrier_recovery_viterbiViterbi_fxp_mex(P, cfg);
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
    %  Tests — one per frequency-recovery algorithm
    %% ================================================================
    methods (Test)

        function test_fft_search(testCase)
            testCase.runPipeline('fft_search');
        end

        function test_differential_kay(testCase)
            testCase.runPipeline('differential_kay');
        end

    end

    %% ================================================================
    %  Private pipeline implementation
    %% ================================================================
    methods (Access = private)

        function runPipeline(testCase, fr_algo)
            P    = testCase;
            NSNR = length(P.SNR_dB_vec);
            NFO  = length(P.DeltaF_Hz_vec);
            NLW  = length(P.LW_Hz_vec);

            % BER storage: (trial, SNR, DeltaF, LW, PR)
            %   PR index:  1=BPS  2=ViterbiViterbi  3=PilotsOnly
            BER_all = zeros(P.NTrials, NSNR, NFO, NLW, 3);

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

                            [fr_out, pilots, txRefBits, rx_preFR, freq_offset] = ...
                                full_cr.buildChannel(P, SNR_dB, DeltaF_Hz, LW, fr_algo);

                            fprintf('    SNR=%2ddB  freq_offset_est = %+.3f MHz\n', ...
                                SNR_dB, freq_offset/1e6);

                            fr_out_fi   = cast(fr_out,   'like', T_cr.x);
                            pilots_fi   = cast(pilots,   'like', T_cr.x);
                            vvfilter_fi = cast(VVFilter, 'like', T_cr.w);

                            % Diagnostic constellation plots: first trial, first FO/LW, highest SNR
                            doPlot = P.Plot && (tr == 1) && (fi == 1) && (li == 1) && (si == NSNR);
                            % if doPlot
                            %     full_cr.plotConstellation(rx_preFR, ...
                            %         sprintf('Before FR | %s | SNR=%ddB', strrep(fr_algo,'_',' '), SNR_dB));
                            %     full_cr.plotConstellation(fr_out, ...
                            %         sprintf('After FR | %s | SNR=%ddB', strrep(fr_algo,'_',' '), SNR_dB));
                            % end

                            %-- BPS (fxp MEX) --
                            [cr_bps, ~] = carrier_recovery.bps_fxp_mex( ...
                                fr_out_fi, P.BPS_N, P.N_pol, P.BPS_M, P.BPS_B, ...
                                P.BlockLen, double(P.StepSize), pilots_fi, P.PilotThreshold, ...
                                double(P.CordicIts), T_cr);
                            cr_bps = full_cr.resolvePhaseAmbiguity(double(cr_bps), txRefBits);
                            BER_all(tr, si, fi, li, 1) = full_cr.computeBER(cr_bps, txRefBits);
                            % if doPlot
                            %     full_cr.plotConstellation(cr_bps, ...
                            %         sprintf('After BPS fxp | %s | SNR=%ddB', strrep(fr_algo,'_',' '), SNR_dB));
                            % end

                            %-- Viterbi-Viterbi (fxp MEX) --
                            [cr_vv, ~] = carrier_recovery.viterbiViterbi_fxp_mex( ...
                                fr_out_fi, P.N_pol, P.VV_NTaps, vvfilter_fi, ...
                                pilots_fi, P.BlockLen, double(P.StepSize), P.PilotThreshold, ...
                                double(P.CordicIts), T_cr);
                            cr_vv = full_cr.resolvePhaseAmbiguity(double(cr_vv), txRefBits);
                            BER_all(tr, si, fi, li, 2) = full_cr.computeBER(cr_vv, txRefBits);
                            % if doPlot
                            %     full_cr.plotConstellation(cr_vv, ...
                            %         sprintf('After VV fxp | %s | SNR=%ddB', strrep(fr_algo,'_',' '), SNR_dB));
                            % end

                            %-- Pilots only (float) --
                            [cr_po, ~] = carrier_recovery.pilots_only( ...
                                fr_out, P.N_pol, P.BlockLen, pilots);
                            cr_po = full_cr.resolvePhaseAmbiguity(cr_po, txRefBits);
                            BER_all(tr, si, fi, li, 3) = full_cr.computeBER(cr_po, txRefBits);
                            % if doPlot
                            %     full_cr.plotConstellation(cr_po, ...
                            %         sprintf('After PilotsOnly | %s | SNR=%ddB', strrep(fr_algo,'_',' '), SNR_dB));
                            % end

                        end % SNR
                    end % LW
                end % DeltaF
            end % trial

            % Average over trials → [NSNR x NFO x NLW x 3]
            BER = reshape(mean(BER_all, 1), [NSNR, NFO, NLW, 3]);

            % BER = 0 cannot be plotted on a log scale.  Replace with the
            % minimum observable BER given NTrials * bitsPerTrial total bits.
            % CPON subframe: 3712 symbols x N_pol x 2 bits/QPSK symbol.
            nBits    = 3712 * P.N_pol * 2;
            berFloor = 1 / (P.NTrials * nBits);
            BER(BER == 0) = berFloor;

            full_cr.printResults(P, BER, fr_algo);

            if P.Plot
                full_cr.plotResults(P, BER, fr_algo, berFloor);
            end
        end

    end

    %% ================================================================
    %  Static helpers
    %% ================================================================
    methods (Static, Access = private)

        function [fr_out, pilots, txRefBits, rx_preFR, freq_offset] = buildChannel( ...
                P, SNR_dB, DeltaF_Hz, LW, fr_algo)
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

            % Frequency recovery (floating-point)
            rx_preFR = rx;
            switch fr_algo
                case 'fft_search'
                    [fr_out, freq_offset] = freq_recovery.fft_search(rx, training, P.Rs, P.FFT_K);
                case 'differential_kay'
                    [fr_out, freq_offset] = freq_recovery.differential_kay(rx, training, P.Rs);
                otherwise
                    error('full_cr:unknownFR', 'Unknown FR algorithm: %s', fr_algo);
            end

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
            % BER: [NSNR x NFO x NLW x 3]
            pr_names = {'BPS', 'VV', 'PilotsOnly'};
            for fi = 1:length(P.DeltaF_Hz_vec)
                for li = 1:length(P.LW_Hz_vec)
                    fprintf('\n[%s]  DeltaF = %.0f MHz  |  LW = %.0f kHz\n', ...
                        fr_algo, P.DeltaF_Hz_vec(fi)/1e6, P.LW_Hz_vec(li)/1e3);
                    fprintf('  SNR [dB]    : ');
                    fprintf('%7.1f  ', P.SNR_dB_vec);
                    fprintf('\n');
                    for pr = 1:3
                        fprintf('  %-12s: ', pr_names{pr});
                        fprintf('%7.5f  ', BER(:, fi, li, pr)');
                        fprintf('\n');
                    end
                end
            end
        end

        function plotResults(P, BER, fr_algo, berFloor)
            % BER: [NSNR x NFO x NLW x 3]
            NFO = length(P.DeltaF_Hz_vec);
            NLW = length(P.LW_Hz_vec);

            PR_names  = {'BPS', 'Viterbi-Viterbi', 'Pilots Only'};
            PR_styles = {'-s', '--o', ':^'};
            PR_colors = lines(3);

            switch fr_algo
                case 'fft_search',       fr_title = 'FFT Search';
                case 'differential_kay', fr_title = 'Differential + Kay';
                otherwise,               fr_title = strrep(fr_algo, '_', ' ');
            end

            fig_w = max(900,  420 * NLW);
            fig_h = max(600,  360 * NFO);
            figure('Name',     sprintf('Full CR  |  FR: %s', fr_title), ...
                   'Position', [80, 80, fig_w, fig_h],                  ...
                   'Color',    'w');

            for fi = 1:NFO
                for li = 1:NLW
                    ax = subplot(NFO, NLW, (fi - 1) * NLW + li);
                    set(ax, 'YScale', 'log', 'FontSize', 11, 'Box', 'on', 'Color', 'w');
                    hold(ax, 'on');

                    for pr = 1:3
                        semilogy(ax, P.SNR_dB_vec, BER(:, fi, li, pr), ...
                            PR_styles{pr}, 'LineWidth', 1.8,            ...
                            'Color',       PR_colors(pr, :),            ...
                            'DisplayName', PR_names{pr});
                    end

                    % FEC threshold
                    yline(ax, 2e-2, 'k--', 'LineWidth', 1.2, ...
                        'DisplayName', 'FEC limit (2\times10^{-2})');

                    % Zero-BER floor (minimum observable BER)
                    yline(ax, berFloor, 'Color', [0.5 0.5 0.5], ...
                        'LineStyle', '--', 'LineWidth', 1.0, ...
                        'DisplayName', sprintf('Zero-error floor (%.2g)', berFloor));

                    grid(ax, 'on');
                    xlabel(ax, 'SNR [dB]', 'FontSize', 11);
                    ylabel(ax, 'BER',      'FontSize', 11);
                    title(ax,  sprintf('\\DeltaF = %.0f MHz,  LW = %.0f kHz', ...
                        P.DeltaF_Hz_vec(fi)/1e6, P.LW_Hz_vec(li)/1e3), ...
                        'FontSize', 11);
                    legend(ax, 'Location', 'northeast', 'FontSize', 9);
                end
            end

            sgtitle(sprintf('BER vs SNR  |  FR: %s  |  %d trials per point', ...
                fr_title, P.NTrials), 'FontSize', 14, 'FontWeight', 'bold');
        end

    end
end
