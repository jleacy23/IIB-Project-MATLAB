classdef test_phase_recovery_pilotthreshold < matlab.unittest.TestCase
%TEST_PHASE_RECOVERY_PILOTTHRESHOLD
%   PilotThreshold sweep for floating-point VV and BPS carrier recovery.
%   For each SNR, BER vs PilotThreshold is plotted with one curve per
%   linewidth.  All (SNR_i, LW_j) combinations are evaluated.
%
%   Run with:
%       runtests('test_phase_recovery_pilotthreshold')

    %% ================================================================
    %  Constant Parameters
    %% ================================================================
    properties (Constant)

        % Modulation / system
        M          = 4
        N_pol      = 2
        Ns         = 2^11          % symbols per polarisation per trial
        Rs         = 30.5            % symbol rate [GBd]
        NTrials    = 20

        % SNR and linewidth grids
        SNR_dB_vec = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]                      % [dB]
        LW_vec     = [100e3, 1000e3, 10000e3]          % [Hz]

        % PilotThreshold sweep
        PilotThresholds = linspace(pi/8, 2*pi, 16)          % [rad]

        % Block / step
        BlockLen   = 32
        StepSize   = 1

        % VV
        VV_NTaps   = 5

        % BPS
        BPS_N      = 5
        BPS_B      = 64

        % Plot
        Plot       = true
    end

    %% ================================================================
    %  Shared Resources (pre-computed VV filters)
    %% ================================================================
    properties
        VVFilters   % cell {NSNR x NLW}, one VVFilter per (SNR, LW) pair
    end

    %% ================================================================
    %  Test Class Setup
    %% ================================================================
    methods (TestClassSetup)

        function seedRng(~)
            rng(42);
        end

        function precomputeVVFilters(testCase)
            % Compute VV Wiener filter for every (SNR, LW) combination.
            % symEnergy is estimated from one sample CPON subframe so that
            % the filter design matches the actual mixed pilot/data energy.
            tmp_bits  = modem.randomBits(4 * 3712);
            [tmp_syms, ~, ~, ~] = modem.modulate(tmp_bits);
            symEnergy = mean(abs(tmp_syms(:)).^2);

            NSNR = length(testCase.SNR_dB_vec);
            NLW  = length(testCase.LW_vec);
            testCase.VVFilters = cell(NSNR, NLW);

            for si = 1:NSNR
                for li = 1:NLW
                    testCase.VVFilters{si, li} = carrier_recovery.genVVFilter( ...
                        testCase.LW_vec(li), testCase.Rs,         ...
                        testCase.SNR_dB_vec(si), symEnergy,       ...
                        testCase.N_pol, testCase.VV_NTaps);
                end
            end
        end

    end

    %% ================================================================
    %  Test
    %% ================================================================
    methods (Test)

        function test_pilotthreshold_sweep(testCase)

            P          = testCase;
            Thresholds = P.PilotThresholds;
            SNR_dB_vec = P.SNR_dB_vec;
            LW_vec     = P.LW_vec;
            NT         = length(Thresholds);
            NSNR       = length(SNR_dB_vec);
            NLW        = length(LW_vec);

            % Results: (trial, threshold, SNR, LW)
            BER_VV_all  = zeros(P.NTrials, NT, NSNR, NLW);
            BER_BPS_all = zeros(P.NTrials, NT, NSNR, NLW);

            for tr = 1:P.NTrials
                fprintf('--- Trial %d / %d ---\n', tr, P.NTrials);

                for si = 1:NSNR
                    for li = 1:NLW

                        SNR_dB   = SNR_dB_vec(si);
                        LW       = LW_vec(li);
                        VVFilter = testCase.VVFilters{si, li};

                        % One channel realisation shared across threshold sweep
                        [~, pilots, txRefBits, rxSym] = ...
                            test_phase_recovery_pilotthreshold.buildChannel( ...
                                P.Ns, SNR_dB, LW, P.Rs, P.BlockLen, P.N_pol);

                        for ti = 1:NT
                            thresh = Thresholds(ti);

                            % fprintf('  SNR=%.1f dB | LW=%.0f kHz | Thresh=%.4f rad\n', ...
                            %     SNR_dB, LW/1e3, thresh);

                            %% VV (floating-point)
                            [cr_vv, ~] = carrier_recovery.viterbiViterbi( ...
                                rxSym, P.N_pol, VVFilter, ...
                                P.BlockLen, P.StepSize, pilots, thresh);

                            BER_VV_all(tr, ti, si, li) = ...
                                test_phase_recovery_pilotthreshold.computeBER( ...
                                    cr_vv, txRefBits);

                            %% BPS (floating-point)
                            [cr_bps, ~] = carrier_recovery.bps( ...
                                rxSym, P.BPS_N, P.N_pol, ...
                                P.M, P.BPS_B, P.BlockLen, P.StepSize, ...
                                pilots, thresh);

                            BER_BPS_all(tr, ti, si, li) = ...
                                test_phase_recovery_pilotthreshold.computeBER( ...
                                    cr_bps, txRefBits);

                        end  % threshold loop
                    end  % LW loop
                end  % SNR loop
            end  % trial loop

            %% Average over trials — always [NT x NSNR x NLW]
            BER_VV  = reshape(mean(BER_VV_all,  1), NT, NSNR, NLW);
            BER_BPS = reshape(mean(BER_BPS_all, 1), NT, NSNR, NLW);

            %% Minimum BER over pilot thresholds → [NSNR x NLW]
            BER_VV_min  = reshape(min(BER_VV,  [], 1), NSNR, NLW);
            BER_BPS_min = reshape(min(BER_BPS, [], 1), NSNR, NLW);

            %% Print table of optimal pilot thresholds
            thresh_deg = Thresholds * 180 / pi;
            sep = repmat('=', 1, 72);
            fprintf('\n%s\n', sep);
            fprintf('Optimal Pilot Thresholds (all thresholds achieving minimum BER)\n');
            fprintf('%s\n', sep);

            for si = 1:NSNR
                for li = 1:NLW
                    fprintf('\nSNR = %5.1f dB  |  LW = %7.0f kHz\n', ...
                        SNR_dB_vec(si), LW_vec(li)/1e3);

                    vv_col  = BER_VV(:, si, li);
                    opt_vv  = thresh_deg(vv_col == BER_VV_min(si, li));
                    fprintf('  VV  | min BER = %.5f | optimal thresholds [deg]:', ...
                        BER_VV_min(si, li));
                    fprintf('  %.1f', opt_vv); fprintf('\n');

                    bps_col = BER_BPS(:, si, li);
                    opt_bps = thresh_deg(bps_col == BER_BPS_min(si, li));
                    fprintf('  BPS | min BER = %.5f | optimal thresholds [deg]:', ...
                        BER_BPS_min(si, li));
                    fprintf('  %.1f', opt_bps); fprintf('\n');
                end
            end
            fprintf('%s\n', sep);

            %% Plot — single figure: min-BER vs SNR, one line per LW
            if P.Plot

                LW_colors = lines(NLW);

                figure('Name',     'CR Min-BER vs SNR | Optimal Pilot Threshold', ...
                       'Position', [100, 100, 900, 600], ...
                       'Color',    'w');
                ax = axes;
                set(ax, 'FontSize', 14, 'LineWidth', 1, 'Box', 'on', ...
                    'Color', 'w', 'YScale', 'log');
                hold on;

                for li = 1:NLW
                    semilogy(ax, SNR_dB_vec, BER_VV_min(:, li), 's-', ...
                        'LineWidth', 2, 'Color', LW_colors(li, :), ...
                        'DisplayName', sprintf('VV  LW=%.0f kHz', LW_vec(li)/1e3));
                    semilogy(ax, SNR_dB_vec, BER_BPS_min(:, li), 'd--', ...
                        'LineWidth', 2, 'Color', LW_colors(li, :), ...
                        'DisplayName', sprintf('BPS LW=%.0f kHz', LW_vec(li)/1e3));
                end

                grid on;
                xlabel('SNR [dB]',                          'Interpreter', 'latex', 'FontSize', 16);
                ylabel('Min BER over pilot thresholds',     'Interpreter', 'latex', 'FontSize', 16);
                legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);
                title(sprintf('Min-BER vs SNR | QPSK | %d trials  (solid=VV, dashed=BPS)', ...
                    P.NTrials));

            end  % Plot
        end

    end

    %% ================================================================
    %  Private Static Helpers
    %% ================================================================
    methods (Static, Access = private)

        function [symbols, pilots, txRefBits, rxSym] = buildChannel( ...
                Ns, SNR_dB, LW, Rs, BlockLen, N_pol) %#ok<INUSD>
            % Generate a CPON signal and build a correct per-CR-block pilot matrix.
            %
            % With BlockLen=32 = CPON_BLOCK_LEN each CR block maps 1:1 to a
            % CPON block.  For any other BlockLen, the pilot at CR block b is
            % the CPON pilot at the signal position (b-1)*BlockLen+1.
            CPON_BLOCK_LEN = 32;
            CPON_SF_SYMS   = 3712;
            N_CPON_BLOCKS  = 116;

            Nbits  = 4 * Ns;
            txBits = modem.randomBits(Nbits);
            [symbols, pilotSyms, ~, ~] = modem.modulate(txBits);

            % Build per-CR-block pilot matrix [NBlocks x N_pol]
            Nsym    = size(symbols, 1);
            NBlocks = ceil(Nsym / BlockLen);
            pilots  = zeros(NBlocks, N_pol);
            for b = 1:NBlocks
                pos       = (b - 1) * BlockLen + 1;
                posInSf   = mod(pos - 1, CPON_SF_SYMS) + 1;
                cponBlock = floor((posInSf - 1) / CPON_BLOCK_LEN) + 1;
                cponBlock = min(cponBlock, N_CPON_BLOCKS);
                pilots(b, :) = pilotSyms(cponBlock, :);
            end

            txRefBits = modem.symbolsToBits(symbols);
            rxSym     = channel.add_awgn(symbols, SNR_dB);
            rxSym     = channel.add_phase_noise(rxSym, Rs, LW);
        end

        function BER = computeBER(crSym, refBits)
            decidedSyms = modem.decideSymbols(crSym);
            rxBits      = modem.symbolsToBits(decidedSyms);
            nBits       = min(length(refBits), length(rxBits));
            BER         = sum(refBits(1:nBits) ~= rxBits(1:nBits)) / nBits;
        end

    end
end