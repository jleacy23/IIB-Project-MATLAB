classdef data < matlab.unittest.TestCase
%data  Compare Tretter-Kay, FFT-search and LRP
%   frequency estimators over a range of offsets and SNRs.
%
%   For each selected frequency offset one figure is produced showing the
%   normalised RMSE (RMSE / Rs) versus SNR for all three algorithms.
%   Each point is averaged over NTrials independent AWGN realisations.
%
%   Run with:
%       results = runtests('data');

    % ================================================================
    %  Parameters
    % ================================================================
    properties (Constant)
        Rs        = 1e-9         % symbol rate normalized to 1 symbol/s
        NTrials   = 1000                % independent noise trials per point
        SNR_dB_vec  = 0:1:20             % SNR sweep [dB]
        % Frequency offset sweep [normalised to Rs = 1] 
        DeltaF_vec     = [0.3]
        % Subset of offsets for which individual SNR-sweep figures are produced
        PlotDeltaF_vec = [0.3]

        % Zero-padding factor K for the floating-point fft_search
        % K * TrainingLen(11) ≈ 1024
        FR_FFT_K      = 40
    end

    methods (TestMethodSetup)
        function seedRng(~)
            rng(42);
        end
    end

    % ================================================================
    %  Test
    % ================================================================
    methods (Test)

        function testMSEComparison(testCase)
            %TESTMSECOMPARISON
            %   Sweeps SNR × DeltaF, runs NTrials per point, plots MSE.

            addpath(fullfile(fileparts(mfilename('fullpath')), ...
                            '..', '..', 'src'));

            Rs_      = testCase.Rs;
            DeltaF_v = testCase.DeltaF_vec;
            SNR_v    = testCase.SNR_dB_vec;
            NT       = testCase.NTrials;
            NF       = numel(DeltaF_v);
            NSNR     = numel(SNR_v);

            % Results: NMSE_alg(SNR_idx, DeltaF_idx)
            NMSE_fft = zeros(NSNR, NF);
            NMSE_dk  = zeros(NSNR, NF);

            % Modified Cramér-Rao bound (MCRB) for frequency estimation
            N_train   = 11;   % TrainingLen
            SNR_lin   = 10.^(SNR_v(:)' ./ 10);
            NMSE_MCRB = 3 * 0.5 ./ (2 * pi^2 * N_train^3 .* SNR_lin); % 0.5 for 2-pol average

            % Pre-compute indices of the plot frequencies inside DeltaF_v
            PlotDeltaF_v    = testCase.PlotDeltaF_vec;
            NPlot           = numel(PlotDeltaF_v);
            plot_df_indices = zeros(1, NPlot);
            for pi_ = 1:NPlot
                [~, plot_df_indices(pi_)] = min(abs(DeltaF_v - PlotDeltaF_v(pi_)));
            end

            for si = 1:NSNR
                SNR_dB = SNR_v(si);
                fprintf('SNR = %d dB\n', SNR_dB);

                for fi = 1:NF
                    df = DeltaF_v(fi);

                    se_fft = zeros(NT, 1);
                    se_dk  = zeros(NT, 1);

                    for tr = 1:NT
                        % Random ±3±3j training symbols (16-QAM corners), 11 x 2
                        training = (2*randi([0 1], 11, 2) - 1)*3 + 1j*(2*randi([0 1], 11, 2) - 1)*3;

                        % Apply LO shift and AWGN to training block only
                        % lo_freq_shift expects DeltaF in MHz; df is in Hz (Rs=1)
                        rx_shifted = channel.lo_freq_shift(training, df * 1e-6, Rs_, 1);
                        rx         = channel.add_awgn(rx_shifted, SNR_dB);

                        [~, est_fft] = freq_recovery.fft_search( ...
                            rx, training, Rs_, testCase.FR_FFT_K);
                        [~, est_dk]  = freq_recovery.differential_kay( ...
                            rx, training, Rs_);

                        se_fft(tr) = (est_fft - df)^2;
                        se_dk(tr)  = (est_dk  - df)^2;
                    end

                    NMSE_fft(si, fi) = mean(se_fft);
                    NMSE_dk(si,  fi) = mean(se_dk);
                end
            end

            % ============================================================
            %  Plot — one figure per selected frequency, SNR on x-axis
            % ============================================================
            colors   = lines(2);
            algNames = {'FFT search (float)', 'Diff+Kay (float)'};

            for pi_ = 1:NPlot
                df_idx    = plot_df_indices(pi_);
                df_actual = DeltaF_v(df_idx);

                figure('Name', sprintf('Freq Recovery NMSE | df = %g MHz', df_actual), ...
                       'Position', [60 + (pi_-1)*40, 60 + (pi_-1)*40, 820, 520], ...
                       'Color', 'w');

                semilogy(SNR_v, NMSE_fft(:, df_idx), '-',  'Color', colors(1,:), 'LineWidth', 1.8, 'DisplayName', algNames{1});
                hold on;
                semilogy(SNR_v, NMSE_dk(:,  df_idx), '--', 'Color', colors(2,:), 'LineWidth', 1.8, 'DisplayName', algNames{2});
                semilogy(SNR_v, NMSE_MCRB,            'k-', 'LineWidth', 2.0,     'DisplayName', sprintf('MCRB (N=%d)', N_train));
                hold off;

                grid on;
                set(gca, 'FontSize', 13, 'LineWidth', 1, 'Box', 'on');
                xlabel('SNR [dB]', 'FontSize', 14);
                ylabel('Normalised RMSE  (RMSE / R_s)  [-]', 'FontSize', 14);
                legend('Location', 'best', 'FontSize', 12);
                title(sprintf('Normalised Frequency Estimation MSE  |  \Deltaf = %g MHz  |  %d trials', ...
                              df_actual, NT));
            end

            % Print summary table
            fprintf('\n%-14s', 'DeltaF [MHz]');
            for si = 1:NSNR
                fprintf('  SNR=%ddB FFT      DiffKay', SNR_v(si));
            end
            fprintf('\n');
            for fi = 1:NF
                fprintf('%-14.0f', DeltaF_v(fi));
                for si = 1:NSNR
                    fprintf('  %8.2e %8.2e', NMSE_fft(si,fi), NMSE_dk(si,fi));
                end
                fprintf('\n');
            end
        end

    end

    % ================================================================
    %  Private static helpers
    % ================================================================
    methods (Static, Access = private)

        function [symbols, training] = generateSubframe()
            %GENERATESUBFRAME  Produce one CPON subframe and its training sequence.
            DATA_PER_SUBFRAME = 3586;
            Nbits = DATA_PER_SUBFRAME * 2 * 2;  % 2 pol, 2 bits/sym/pol

            bits = modem.randomBits(Nbits);
            [symbols_full, ~, training, ~] = modem.modulate(bits);

            SUBFRAME_SYMS = 3712;
            symbols = symbols_full(1:SUBFRAME_SYMS, :);  % [3712 x 2]
            % training is [11 x 2], same for every subframe
        end

    end
end
