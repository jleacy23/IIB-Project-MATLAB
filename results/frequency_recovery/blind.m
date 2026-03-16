classdef blind < matlab.unittest.TestCase
%blind  Blind frequency estimator MSE comparison.
%
%   Runs Tretter-Kay, FFT-search and LRP in blind (4th-power) mode over a
%   range of data observation lengths D and SNRs.  For each selected
%   frequency offset one figure is produced showing the MSE vs SNR for all
%   three algorithms at every value of D, together with the MCRB for each D.
%
%   Algorithms use the floating-point implementations directly (no MEX
%   required).  The 4th-power operation removes unknown data modulation;
%   the estimate is divided by 4 inside each estimator.
%
%   Run with:
%       results = runtests('blind');

    % ================================================================
    %  Parameters
    % ================================================================
    properties (Constant)
        Rs             = 1e-9          % symbol rate (normalised, same as training-aided test)
        NTrials        = 1000          % independent noise trials per point
        SNR_dB_vec     = 0:2:24        % SNR sweep [dB]
        DeltaF_vec     = [0.00]         % frequency offset(s) to sweep
        PlotDeltaF_vec = [0.00]         % subset of offsets for which figures are produced
        D_vec          = [16, 32, 64, 128]   % blind observation lengths [symbols]
        TrainingLen    = 11            % training block length (sets start of data block)
        FR_FFT_K       = 46             % zero-padding factor for fft_search
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

        function testBlindMSE(testCase)
            %TESTBLINDMSE  Sweep SNR × DeltaF × D in blind mode, plot MSE.

            addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src'));

            Rs_      = testCase.Rs;
            DeltaF_v = testCase.DeltaF_vec;
            SNR_v    = testCase.SNR_dB_vec;
            NT       = testCase.NTrials;
            L        = testCase.TrainingLen;
            D_v      = testCase.D_vec;
            ND       = numel(D_v);
            NF       = numel(DeltaF_v);
            NSNR     = numel(SNR_v);
            D_max    = max(D_v);

            % MSE arrays [NSNR x NF x ND]
            NMSE_fft = zeros(NSNR, NF, ND);
            NMSE_dk  = zeros(NSNR, NF, ND);

            % Modified CRB for frequency estimation from D constant-envelope
            % symbols via 4th-power method (reduces to same formula as
            % training-aided CRB with N->D, after the 1/4 scaling cancels):
            %   MCRB(f0) = 3 / (2*pi^2 * D^3 * SNR)
            % The 0.5 factor accounts for averaging over 2 polarisations.
            SNR_lin   = 10.^(SNR_v(:)' ./ 10);
            NMSE_MCRB = zeros(ND, NSNR);
            for di = 1:ND
                NMSE_MCRB(di, :) = 3 * 0.5 ./ (2 * pi^2 * D_v(di)^3 .* SNR_lin);
            end

            % Pre-compute indices of plot frequencies inside DeltaF_v
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

                    se_fft = zeros(NT, ND);
                    se_dk  = zeros(NT, ND);

                    for tr = 1:NT
                        % Random ±1±1j symbols (constant envelope) for training
                        % and data blocks
                        training  = (2*randi([0 1], L,     2) - 1)*1 + 1j*(2*randi([0 1], L,     2) - 1)*1;
                        data_syms = (2*randi([0 1], D_max, 2) - 1)*1 + 1j*(2*randi([0 1], D_max, 2) - 1)*1;
                        x_clean   = [training; data_syms];   % (L+D_max) x 2

                        % Apply LO shift and AWGN to the full block once
                        % lo_freq_shift expects DeltaF in MHz; df is in Hz (Rs=1)
                        rx_shifted = channel.lo_freq_shift(x_clean, df * 1e-6, Rs_, 1);
                        rx         = channel.add_awgn(rx_shifted, SNR_dB);

                        for di = 1:ND
                            D_val   = D_v(di);
                            x_block = rx(1:L+D_val, :);   % (L+D_val) x 2

                            % Blind mode: data_aided=false, D=D_val
                            % training is only used to determine L (offset to data block)
                            [~, est_fft] = freq_recovery.fft_search( ...
                                x_block, training, Rs_, testCase.FR_FFT_K, false, D_val);
                            [~, est_dk]  = freq_recovery.differential_kay( ...
                                x_block, training, Rs_, false, D_val);

                            se_fft(tr, di) = (est_fft - df)^2;
                            se_dk(tr,  di) = (est_dk  - df)^2;
                        end
                    end

                    for di = 1:ND
                        NMSE_fft(si, fi, di) = mean(se_fft(:, di));
                        NMSE_dk(si,  fi, di) = mean(se_dk(:,  di));
                    end
                end
            end

            % ============================================================
            %  Plot — one figure per algorithm per selected frequency
            %  Each figure shows curves for all D values + matching MCRBs
            %  Colours  → D value
            %  Solid    → algorithm estimate,  dashed → MCRB
            % ============================================================
            alg_names = {'FFT search', 'Diff+Kay'};
            D_colors  = lines(ND);        % one colour per D value

            for pi_ = 1:NPlot
                df_idx    = plot_df_indices(pi_);
                df_actual = DeltaF_v(df_idx);

                for ai = 1:2   % algorithm index
                    figure('Name', sprintf('Blind %s MSE | df = %g', alg_names{ai}, df_actual), ...
                           'Position', [60 + (ai-1)*30 + (pi_-1)*20, ...
                                        60 + (ai-1)*30 + (pi_-1)*20, 820, 520], ...
                           'Color', 'w');
                    hold on;

                    for di = 1:ND
                        D_val = D_v(di);
                        c     = D_colors(di, :);

                        switch ai
                            case 1; data = NMSE_fft(:, df_idx, di);
                            case 2; data = NMSE_dk(:,  df_idx, di);
                        end

                        semilogy(SNR_v, data,              '-',  'Color', c, 'LineWidth', 1.8, ...
                                 'DisplayName', sprintf('D=%d',      D_val));
                        semilogy(SNR_v, NMSE_MCRB(di, :), '--', 'Color', c, 'LineWidth', 1.2, ...
                                 'DisplayName', sprintf('MCRB D=%d', D_val));
                    end

                    hold off;
                    grid on;
                    set(gca, 'FontSize', 13, 'LineWidth', 1, 'Box', 'on', 'YScale', 'log');
                    xlabel('SNR [dB]', 'FontSize', 14);
                    ylabel('Normalised MSE', 'FontSize', 14);
                    legend('Location', 'best', 'FontSize', 11);
                    title(sprintf('%s', ...
                                  alg_names{ai}));
                end
            end

            % ============================================================
            %  Summary table
            % ============================================================
            fprintf('\n%-10s %-6s', 'DeltaF', 'D');
            for si = 1:NSNR
                fprintf('  SNR=%ddB   MCRB      FFT       DiffKay', SNR_v(si));
            end
            fprintf('\n');
            for fi = 1:NF
                for di = 1:ND
                    fprintf('%-10.4f %-6d', DeltaF_v(fi), D_v(di));
                    for si = 1:NSNR
                        fprintf('  %8.2e  %8.2e  %8.2e', ...
                            NMSE_MCRB(di,si), NMSE_fft(si,fi,di), NMSE_dk(si,fi,di));
                    end
                    fprintf('\n');
                end
            end
        end

    end

end
