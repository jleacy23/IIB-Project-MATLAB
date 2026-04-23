classdef blind_data < matlab.unittest.TestCase
%BLIND_DATA  SNR threshold at FEC limit vs blind observation length.
%
%   Sweeps the blind data-observation length D and plots, for each FR
%   algorithm separately, the SNR required to achieve BER = 2e-2 with
%   both Viterbi-Viterbi and pilots-only phase recovery.  Data-aided
%   results are shown as horizontal dashed reference lines.
%
%   FR algorithms:  fft_search (Nfft fixed at FR_Nfft for all D),
%                   differential_kay  — both as blind and data-aided
%   PR algorithms:  viterbiViterbi, pilots_only  (fixed-point MEX)
%   Output:         two figures, one per FR algorithm
%
%   Run with:
%       runtests('blind_data')

    %% ================================================================
    %  Parameters
    %% ================================================================
    properties (Constant)

        % System
        Rs          = 30.5              % symbol rate [GBd]
        N_pol       = 2
        TrainingLen = 11

        % Monte-Carlo
        NTrials     = 30                % trials per (SNR, D) point

        % Sweep grids
        SNR_dB_vec = 0 : 1 : 20        % [dB]
        BlindD_vec = [16, 32, 64, 128, 256, 512]  % observation lengths [symbols]

        % Channel conditions (single operating point)
        DeltaF_Hz      = 3e9            % frequency offset [Hz]
        LW_Hz          = 1000e3         % laser linewidth  [Hz]

        % Fixed-point configuration (word length / fraction length)
        FxpConfig_FR   = struct('WL', 32, 'FL', 16)   % frequency recovery
        FxpConfig_CR   = struct('WL', 32, 'FL', 16)  % carrier recovery (VV + PO)
        CordicIts      = 16

        % FFT search parameters — Nfft is constant across all D values
        FR_Nfft        = 512
        FR_Po2Twiddle  = false
        MaxFreq        = 0.1

        % Phase recovery
        BlockLen       = 32
        StepSize       = 32
        PilotThreshold = 5 * pi / 9    % cycle-slip detection threshold [rad]
        VV_NTaps       = 10

        % Build control
        Rebuild        = true

        % FEC threshold
        FEC_BER        = 2e-2

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

        function buildMex(testCase)
            buildDir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'build');
            addpath(buildDir);

            B.Rs             = testCase.Rs;
            B.N_pol          = testCase.N_pol;
            B.TrainingLen    = testCase.TrainingLen;
            B.FR_Nfft        = testCase.FR_Nfft;
            B.FR_Po2Twiddle  = testCase.FR_Po2Twiddle;
            B.FR_BlindD      = max(testCase.BlindD_vec);
            B.FxpConfig_FR   = testCase.FxpConfig_FR;
            B.FxpConfig_VV   = testCase.FxpConfig_CR;
            B.FxpConfig_PO   = testCase.FxpConfig_CR;
            B.CordicIts      = testCase.CordicIts;
            B.VV_NTaps       = testCase.VV_NTaps;
            B.BlockLen       = testCase.BlockLen;
            B.StepSize       = testCase.StepSize;
            B.PilotThreshold = testCase.PilotThreshold;
            B.PilotLen       = 1;
            B.MaxFreq        = testCase.MaxFreq;

            cfg = coder.config('mex');
            cfg.GenerateReport = false;

            if testCase.Rebuild
                fprintf('Building MEX objects...\n');
                build_freq_recovery_fft_search_fxp_mex(B, cfg);
                build_freq_recovery_differential_kay_fxp_mex(B, cfg);
                build_carrier_recovery_viterbiViterbi_fxp_mex(B, cfg);
                build_carrier_recovery_pilots_only_fxp_mex(B, cfg);
                fprintf('All MEX objects built.\n');
            end
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

        function test_snr_vs_blind_length(testCase)
            P  = testCase;
            ND = length(P.BlindD_vec);

            % fecSNR: [ND x 2]  columns = VV, PO
            fecSNR_fft_blind = nan(ND, 2);
            fecSNR_dk_blind  = nan(ND, 2);
            fecSNR_fft_DA    = nan(1, 2);
            fecSNR_dk_DA     = nan(1, 2);

            % ---- Data-aided references ----------------------------------
            fprintf('Data-aided reference (FFT search)...\n');
            ber = testCase.runSnrSweep('fft_search', 0);
            fecSNR_fft_DA(1) = blind_data.fecCrossing(P.SNR_dB_vec, ber(:,1), P.FEC_BER);
            fecSNR_fft_DA(2) = blind_data.fecCrossing(P.SNR_dB_vec, ber(:,2), P.FEC_BER);

            fprintf('Data-aided reference (differential Kay)...\n');
            ber = testCase.runSnrSweep('differential_kay', 0);
            fecSNR_dk_DA(1) = blind_data.fecCrossing(P.SNR_dB_vec, ber(:,1), P.FEC_BER);
            fecSNR_dk_DA(2) = blind_data.fecCrossing(P.SNR_dB_vec, ber(:,2), P.FEC_BER);

            % ---- Blind sweep over D -------------------------------------
            for di = 1:ND
                D = P.BlindD_vec(di);

                fprintf('[FFT search]       blind D = %4d  (%d/%d)\n', D, di, ND);
                ber = testCase.runSnrSweep('fft_search_blind', D);
                fecSNR_fft_blind(di,1) = blind_data.fecCrossing(P.SNR_dB_vec, ber(:,1), P.FEC_BER);
                fecSNR_fft_blind(di,2) = blind_data.fecCrossing(P.SNR_dB_vec, ber(:,2), P.FEC_BER);

                fprintf('[Differential Kay] blind D = %4d  (%d/%d)\n', D, di, ND);
                ber = testCase.runSnrSweep('differential_kay_blind', D);
                fecSNR_dk_blind(di,1) = blind_data.fecCrossing(P.SNR_dB_vec, ber(:,1), P.FEC_BER);
                fecSNR_dk_blind(di,2) = blind_data.fecCrossing(P.SNR_dB_vec, ber(:,2), P.FEC_BER);
            end

            % ---- Plot ---------------------------------------------------
            blind_data.plotAlgoResults(P, ...
                sprintf('FFT Search  (N_{fft} = %d,  fixed)', P.FR_Nfft), ...
                fecSNR_fft_blind, fecSNR_fft_DA);

            blind_data.plotAlgoResults(P, ...
                'Differential Kay', ...
                fecSNR_dk_blind, fecSNR_dk_DA);
        end

    end

    %% ================================================================
    %  Private instance helpers
    %% ================================================================
    methods (Access = private)

        function ber_avg = runSnrSweep(testCase, fr_algo, blindD)
            % Returns BER averaged over NTrials: [NSNR x 2] (col1=VV, col2=PO).
            P       = testCase;
            NSNR    = length(P.SNR_dB_vec);
            ber_all = zeros(P.NTrials, NSNR, 2);

            T_fr = freq_recovery.fxp_types(P.FxpConfig_FR);
            T_cr = carrier_recovery.fxp_types(P.FxpConfig_CR);

            for tr = 1:P.NTrials
                for si = 1:NSNR
                    [fr_out, pilots, txRefBits] = blind_data.buildChannel( ...
                        P, P.SNR_dB_vec(si), fr_algo, blindD, T_fr);

                    fr_fi     = cast(fr_out,              'like', T_cr.x);
                    pilots_fi = cast(pilots,               'like', T_cr.x);
                    vvfilt_fi = cast(testCase.VVFilters{si}, 'like', T_cr.w);

                    % Viterbi-Viterbi
                    [cr_vv, ~] = carrier_recovery.viterbiViterbi_fxp_mex( ...
                        fr_fi, P.N_pol, P.VV_NTaps, vvfilt_fi, pilots_fi, ...
                        P.BlockLen, double(P.StepSize), P.PilotThreshold, ...
                        double(P.CordicIts), T_cr);
                    cr_vv = blind_data.resolveAmbiguity(double(cr_vv), txRefBits);
                    ber_all(tr, si, 1) = blind_data.computeBER(cr_vv, txRefBits);

                    % Pilots-only
                    [cr_po, ~] = carrier_recovery.pilots_only_fxp_mex( ...
                        fr_fi, P.N_pol, P.BlockLen, pilots_fi, ...
                        double(P.CordicIts), T_cr);
                    cr_po = blind_data.resolveAmbiguity(double(cr_po), txRefBits);
                    ber_all(tr, si, 2) = blind_data.computeBER(cr_po, txRefBits);
                end
            end

            ber_avg = reshape(mean(ber_all, 1), [NSNR, 2]);

            nBits    = 3712 * P.N_pol * 2;
            berFloor = 1 / (P.NTrials * nBits);
            ber_avg(ber_avg == 0) = berFloor;
        end

    end

    %% ================================================================
    %  Static helpers
    %% ================================================================
    methods (Static, Access = private)

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
                case 'differential_kay_blind'
                    [fr_fi, ~] = freq_recovery.differential_kay_fxp_mex( ...
                        rx_fi, tr_fi, P.Rs, P.CordicIts, T_fr, false, blindD, P.MaxFreq);
                otherwise
                    error('blind_data:unknownFR', 'Unknown FR algorithm: %s', fr_algo);
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

        function plotAlgoResults(P, algoTitle, fecSNR_blind, fecSNR_DA)
            % fecSNR_blind : [ND x 2]  columns = VV, PO
            % fecSNR_DA    : [1  x 2]  same order — plotted as horizontal lines

            prNames = {'Viterbi-Viterbi', 'Pilots-only'};
            colors  = lines(2);
            markers = {'o', 's'};

            figure('Name', sprintf('SNR at FEC vs Blind D  |  %s', algoTitle), ...
                'Position', [100 100 780 520], 'Color', 'w');
            ax = axes;
            hold(ax, 'on');
            grid(ax, 'on');

            % Blind curves
            for c = 1:2
                valid = ~isnan(fecSNR_blind(:, c));
                if any(valid)
                    plot(ax, P.BlindD_vec(valid), fecSNR_blind(valid, c), ...
                        'LineStyle', '-', 'Marker', markers{c}, ...
                        'MarkerSize', 6, 'LineWidth', 1.8, ...
                        'Color', colors(c, :), ...
                        'DisplayName', sprintf('%s (blind)', prNames{c}));
                end
            end

            % Data-aided horizontal reference lines
            for c = 1:2
                if ~isnan(fecSNR_DA(c))
                    yline(ax, fecSNR_DA(c), ...
                        'LineStyle', '--', 'LineWidth', 1.2, ...
                        'Color', colors(c, :), ...
                        'HandleVisibility', 'off');
                    plot(ax, NaN, NaN, ...
                        'LineStyle', '--', 'LineWidth', 1.2, ...
                        'Color', colors(c, :), ...
                        'DisplayName', sprintf('%s (data-aided)', prNames{c}));
                end
            end

            set(ax, 'XScale', 'log', 'FontSize', 11, 'Box', 'on');
            xlim(ax, [P.BlindD_vec(1) * 0.7, P.BlindD_vec(end) * 1.4]);

            xlabel(ax, 'Blind observation length  D  [symbols]', 'FontSize', 12);
            ylabel(ax, sprintf('SNR at BER = %.0e  [dB]', P.FEC_BER), 'FontSize', 12);
            title(ax, sprintf('%s\n\\DeltaF = %.0f MHz,  LW = %.0f kHz', ...
                algoTitle, P.DeltaF_Hz/1e6, P.LW_Hz/1e3), 'FontSize', 12);

            lgd = legend(ax, 'Location', 'northeast', 'FontSize', 10);
            lgd.Title.String = 'Phase recovery  (solid = blind, dashed = DA)';
        end

    end
end
