classdef test_phase_recovery_pilotthreshold < matlab.unittest.TestCase
%TEST_PHASE_RECOVERY_PILOTTHRESHOLD
%   PilotThreshold sweep for fixed-point VV and BPS CR.
%   For each SNR, BER vs PilotThreshold is plotted with one line per
%   linewidth.  All (SNR_i, LW_j) combinations are evaluated.

    %% ================================================================
    %  Constant Parameters
    %% ================================================================
    properties (Constant)

        % Modulation / system
        M           = 4
        N_pol       = 2
        Ns          = 2^16
        Rs          = 10
        NTrials     = 20

        % SNR and linewidth grids
        SNR_dB_vec  = [15.0, 17.5, 20.0]          % dB
        LW_vec      = [600e3, 1200e3, 2400e3, 4800e3]      % Hz

        % Pilot
        PilotLen    = 4
        UsePilots   = true
        PilotThresholds = linspace(pi/8, 7*pi/8, 13)  % sweep range

        % Block
        BlockLen    = 32
        StepSize    = 1

        % VV
        VV_NTaps    = 5

        % BPS
        BPS_N       = 5
        BPS_B       = 64

        % Fixed-point
        FxpConfig   = 'fixed16'
        CordicIts   = 8

        % Plot
        Plot        = true
        PlotTrials  = false
    end

    %% ================================================================
    %  Shared Resources
    %% ================================================================
    properties
        T_vv
        T_bps
        VVFilters_fi    % cell array {1 x length(SNR_dB_vec)} - one per SNR
    end

    %% ================================================================
    %  Test Class Setup
    %% ================================================================
    methods (TestClassSetup)

        function buildMexAndSetup(testCase)

            P = testCase;

            assert(mod(P.BPS_B, 2) == 0, 'BPS_B must be even.');

            %% Fixed-point types
            testCase.T_vv  = carrier_recovery.viterbiViterbi_fxp_types(P.FxpConfig);
            testCase.T_bps = carrier_recovery.bps_fxp_types(P.FxpConfig);

            %% Build configuration (use first SNR / first LW / first threshold)
            Pbuild.N_pol          = P.N_pol;
            Pbuild.M              = P.M;
            Pbuild.PilotLen       = P.PilotLen;
            Pbuild.PilotThreshold = P.PilotThresholds(1);
            Pbuild.VV_NTaps       = P.VV_NTaps;
            Pbuild.BPS_N          = P.BPS_N;
            Pbuild.BPS_B          = P.BPS_B;
            Pbuild.BlockLen       = P.BlockLen;
            Pbuild.StepSize       = P.StepSize;
            Pbuild.FxpConfig_VV   = P.FxpConfig;
            Pbuild.FxpConfig_BPS  = P.FxpConfig;
            Pbuild.CordicIts      = P.CordicIts;

            cfg = coder.config('mex');
            cfg.GenerateReport = false;
            cfg.SaturateOnIntegerOverflow = false;

            fprintf('=== Compiling MEX binaries ===\n');
            build_carrier_recovery_viterbiViterbi_fxp_mex(Pbuild, cfg);
            build_carrier_recovery_bps_fxp_mex(Pbuild, cfg);
            fprintf('=== MEX compilation complete ===\n\n');

            %% Pre-compute VV filters for every (SNR, LW) pair
            NSNR = length(P.SNR_dB_vec);
            NLW  = length(P.LW_vec);
            testCase.VVFilters_fi = cell(NSNR, NLW);
            symEnergy = 1.0;

            for si = 1:NSNR
                for li = 1:NLW
                    VVFilter = carrier_recovery.genVVFilter(P.LW_vec(li), P.Rs, ...
                                              P.SNR_dB_vec(si), ...
                                              symEnergy, P.N_pol, P.VV_NTaps);
                    testCase.VVFilters_fi{si, li} = cast(VVFilter, ...
                                                        'like', testCase.T_vv.w);
                end
            end
        end
    end

    %% ================================================================
    %  Test
    %% ================================================================
    methods (Test)

        function test_pilotthreshold_sweep(testCase)

            P           = testCase;
            Thresholds  = P.PilotThresholds;
            SNR_dB_vec  = P.SNR_dB_vec;
            LW_vec      = P.LW_vec;
            NT          = length(Thresholds);
            NSNR        = length(SNR_dB_vec);
            NLW         = length(LW_vec);

            % Results: (trials, threshold, SNR, LW)
            BER_VV_all  = zeros(P.NTrials, NT, NSNR, NLW);
            BER_BPS_all = zeros(P.NTrials, NT, NSNR, NLW);

            for tr = 1:P.NTrials
                fprintf('--- Trial %d / %d ---\n', tr, P.NTrials);

                for si = 1:NSNR
                    for li = 1:NLW

                        SNR_dB = SNR_dB_vec(si);
                        LW     = LW_vec(li);

                        %% Channel realisation (shared across threshold sweep)
                        Nbits  = 4 * P.Ns;
                        txBits = modem.randomBits(Nbits);
                        [symbols, pilots, ~, ~] = modem.modulate(txBits);
                        txRefBits = modem.symbolsToBits(symbols);
                        rxSym = channel.add_awgn(symbols, SNR_dB);
                        rxSym = channel.add_phase_noise(rxSym, P.Rs, LW);

                        %% Cast
                        rx_vv       = cast(rxSym,  'like', testCase.T_vv.x);
                        rx_bps      = cast(rxSym,  'like', testCase.T_bps.x);
                        pilots_vv   = cast(pilots, 'like', testCase.T_vv.x);
                        pilots_bps  = cast(pilots, 'like', testCase.T_bps.x);

                        VVFilter_fi = testCase.VVFilters_fi{si, li};

                        for ti = 1:NT

                            thresh = Thresholds(ti);
                            fprintf('  SNR=%.1f dB | LW=%.0f kHz | Threshold=%.4f\n', ...
                                SNR_dB, LW/1e3, thresh);

                            %% VV
                            [cr_vv_fi, ~] = carrier_recovery.viterbiViterbi_fxp_mex( ...
                                rx_vv, P.N_pol, P.VV_NTaps, ...
                                VVFilter_fi, ...
                                pilots_vv, P.BlockLen, P.StepSize, ...
                                P.UsePilots, thresh, P.CordicIts, testCase.T_vv);

                            cr_vv = test_phase_recovery_pilotthreshold.resolvePhaseAmbiguity( ...
                                double(cr_vv_fi), txRefBits);

                            BER_VV_all(tr, ti, si, li) = ...
                                test_phase_recovery_pilotthreshold.computeBER( ...
                                cr_vv, txRefBits);

                            %% BPS
                            [cr_bps_fi, ~] = carrier_recovery.bps_fxp_mex( ...
                                rx_bps, P.BPS_N, P.N_pol, ...
                                P.M, P.BPS_B, ...
                                P.BlockLen, P.StepSize, ...
                                pilots_bps, P.UsePilots, thresh, P.CordicIts, ...
                                testCase.T_bps);

                            cr_bps = test_phase_recovery_pilotthreshold.resolvePhaseAmbiguity( ...
                                double(cr_bps_fi), txRefBits);

                            BER_BPS_all(tr, ti, si, li) = ...
                                test_phase_recovery_pilotthreshold.computeBER( ...
                                cr_bps, txRefBits);

                        end  % threshold loop
                    end  % LW loop
                end  % SNR loop
            end  % trial loop

            %% Average over trials
            BER_VV  = squeeze(mean(BER_VV_all,  1));  % (NT, NSNR, NLW)
            BER_BPS = squeeze(mean(BER_BPS_all, 1));

            %% Display numerical results
            for si = 1:NSNR
                for li = 1:NLW
                    fprintf('\nSNR=%.1f dB | LW=%.0f kHz\n', ...
                        SNR_dB_vec(si), LW_vec(li)/1e3);
                    fprintf('  BER_VV  = '); disp(BER_VV(:,si,li)');
                    fprintf('  BER_BPS = '); disp(BER_BPS(:,si,li)');
                end
            end

            %% Plot - one figure per SNR
            if P.Plot

                LW_colors = lines(NLW);

                thresh_deg = Thresholds * 180 / pi;   % x-axis in degrees

                for si = 1:NSNR

                    %% VV figure
                    fig_vv = figure('Name', ...
                        sprintf('VV BER vs PilotThreshold | SNR=%.1f dB', ...
                                SNR_dB_vec(si)), ...
                        'Position', [100 + (si-1)*60, 100 + (si-1)*60, 800, 550], ...
                        'Color', 'w');

                    ax = gca;
                    set(ax, 'FontSize', 14, 'LineWidth', 1, 'Box', 'on', 'Color', 'w');
                    hold on;

                    for li = 1:NLW
                        ber_mean = BER_VV(:, si, li);

                        if P.PlotTrials
                            for ti = 1:NT
                                semilogy(thresh_deg(ti) * ones(P.NTrials, 1), ...
                                    BER_VV_all(:, ti, si, li), 'o', ...
                                    'Color', LW_colors(li,:), ...
                                    'HandleVisibility', 'off');
                            end
                        end

                        semilogy(thresh_deg, ber_mean, 's-', ...
                            'LineWidth', 2, ...
                            'Color', LW_colors(li,:), ...
                            'DisplayName', sprintf('LW=%.0f kHz', LW_vec(li)/1e3));
                    end

                    grid on;
                    xlabel('Pilot threshold [deg]', ...
                        'Interpreter', 'latex', 'FontSize', 16);
                    ylabel('BER', ...
                        'Interpreter', 'latex', 'FontSize', 16);
                    legend('Location', 'best', ...
                        'Interpreter', 'latex', 'FontSize', 14);
                    title(sprintf(['VV | BER vs Pilot Threshold | QPSK | ' ...
                          'SNR=%.1f dB | %d trials'], ...
                          SNR_dB_vec(si), P.NTrials));

                    %% BPS figure
                    fig_bps = figure('Name', ...
                        sprintf('BPS BER vs PilotThreshold | SNR=%.1f dB', ...
                                SNR_dB_vec(si)), ...
                        'Position', [180 + (si-1)*60, 180 + (si-1)*60, 800, 550], ...
                        'Color', 'w');

                    ax = gca;
                    set(ax, 'FontSize', 14, 'LineWidth', 1, 'Box', 'on', 'Color', 'w');
                    hold on;

                    for li = 1:NLW
                        ber_mean = BER_BPS(:, si, li);

                        if P.PlotTrials
                            for ti = 1:NT
                                semilogy(thresh_deg(ti) * ones(P.NTrials, 1), ...
                                    BER_BPS_all(:, ti, si, li), 'o', ...
                                    'Color', LW_colors(li,:), ...
                                    'HandleVisibility', 'off');
                            end
                        end

                        semilogy(thresh_deg, ber_mean, 's-', ...
                            'LineWidth', 2, ...
                            'Color', LW_colors(li,:), ...
                            'DisplayName', sprintf('LW=%.0f kHz', LW_vec(li)/1e3));
                    end

                    grid on;
                    xlabel('Pilot threshold [deg]', ...
                        'Interpreter', 'latex', 'FontSize', 16);
                    ylabel('BER', ...
                        'Interpreter', 'latex', 'FontSize', 16);
                    legend('Location', 'best', ...
                        'Interpreter', 'latex', 'FontSize', 14);
                    title(sprintf(['BPS | BER vs Pilot Threshold | QPSK | ' ...
                          'SNR=%.1f dB | %d trials'], ...
                          SNR_dB_vec(si), P.NTrials));

                end  % SNR figure loop
            end  % Plot
        end
    end

    %% ================================================================
    %  Private Static Helpers
    %% ================================================================
    methods (Static, Access = private)

        function BER = computeBER(crSym, refBits)
            decidedSyms = modem.decideSymbols(crSym);
            rxBits      = modem.symbolsToBits(decidedSyms);
            BER         = sum(refBits ~= rxBits) / length(refBits);
        end

        function bestSym = resolvePhaseAmbiguity(crSym, refBits)
            bestBER = Inf;
            bestSym = crSym;

            for k = 0:3
                rotated     = crSym .* exp(-1j * k * pi/2);
                decidedSyms = modem.decideSymbols(rotated);
                rxBits      = modem.symbolsToBits(decidedSyms);
                thisBER     = sum(refBits ~= rxBits) / length(refBits);

                if thisBER < bestBER
                    bestBER = thisBER;
                    bestSym = rotated;
                end
            end
        end
    end
end