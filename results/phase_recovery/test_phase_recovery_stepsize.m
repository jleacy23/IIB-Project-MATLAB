classdef test_phase_recovery_stepsize < matlab.unittest.TestCase
%TEST_PHASE_RECOVERY_STEPSIZE
%   Step-size sweep (powers of 2) for fixed-point VV and BPS CR.
%   BER vs StepSize is plotted.

    %% ================================================================
    %  Constant Parameters
    %% ================================================================
    properties (Constant)

        % Modulation / system
        M           = 4
        N_pol       = 2
        Ns          = 2^14
        Rs          = 10
        LW          = 2400e3
        SNR_dB      = 17.5
        NTrials     = 5

        % Pilot
        PilotLen    = 4
        UsePilots   = true
        PilotThreshold = 3 * pi/4

        % Block
        BlockLen    = 128

        % Step sizes (powers of 2 up to BlockLen)
        StepSizes   = 2.^(0:log2(128))

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
        PlotTrials   = false
    end

    %% ================================================================
    %  Shared Resources
    %% ================================================================
    properties
        T_vv
        T_bps
        VVFilter_fi
    end

    %% ================================================================
    %  Test Class Setup
    %% ================================================================
    methods (TestClassSetup)

        function buildMexAndSetup(testCase)

            P = testCase;

            assert(mod(P.BPS_B,2)==0,'BPS_B must be even.');

            %% Fixed-point types
            testCase.T_vv  = carrier_recovery.viterbiViterbi_fxp_types(P.FxpConfig);
            testCase.T_bps = carrier_recovery.bps_fxp_types(P.FxpConfig);

            %% Build configuration
            Pbuild.N_pol         = P.N_pol;
            Pbuild.M             = P.M;
            Pbuild.PilotLen      = P.PilotLen;
            Pbuild.PilotThreshold = P.PilotThreshold;
            Pbuild.VV_NTaps      = P.VV_NTaps;
            Pbuild.BPS_N         = P.BPS_N;
            Pbuild.BPS_B         = P.BPS_B;
            Pbuild.BlockLen      = P.BlockLen;
            Pbuild.StepSize      = P.StepSizes(1);
            Pbuild.FxpConfig_VV  = P.FxpConfig;
            Pbuild.FxpConfig_BPS = P.FxpConfig;
            Pbuild.CordicIts     = P.CordicIts;

            cfg = coder.config('mex');
            cfg.GenerateReport = false;
            cfg.SaturateOnIntegerOverflow = false;

            fprintf('=== Compiling MEX binaries ===\n');
            build_carrier_recovery_viterbiViterbi_fxp_mex(Pbuild, cfg);
            build_carrier_recovery_bps_fxp_mex(Pbuild, cfg);
            fprintf('=== MEX compilation complete ===\n\n');

            %% VV filter
            symEnergy = 1.0;
            VVFilter  = carrier_recovery.genVVFilter(P.LW, P.Rs, P.SNR_dB, ...
                                       symEnergy, P.N_pol, P.VV_NTaps);

            testCase.VVFilter_fi = cast(VVFilter, ...
                                        'like', testCase.T_vv.w);
        end
    end

    %% ================================================================
    %  Test
    %% ================================================================
    methods (Test)

        function test_step_size_sweep(testCase)

            P = testCase;
            StepSizes = P.StepSizes;
            NS = length(StepSizes);

            BER_VV_all  = zeros(P.NTrials, NS);
            BER_BPS_all = zeros(P.NTrials, NS);

            for tr = 1:P.NTrials

                %% Channel
                Nbits = 4 * P.Ns;
                txBits = modem.randomBits(Nbits);
                [symbols, pilots, ~, ~] = modem.modulate(txBits);
                txRefBits = modem.symbolsToBits(symbols);
                rxSym = channel.add_awgn(symbols, P.SNR_dB);
                rxSym = channel.add_phase_noise(rxSym, P.Rs, P.LW);
                %% Cast
                rx_vv  = cast(rxSym,'like',testCase.T_vv.x);
                rx_bps = cast(rxSym,'like',testCase.T_bps.x);
                pilots_vv  = cast(pilots,'like',testCase.T_vv.x);
                pilots_bps = cast(pilots,'like',testCase.T_bps.x);

                for si = 1:NS

                    stepSize = StepSizes(si);
                    fprintf('StepSize = %d\n', stepSize);



                    %% VV
                    [cr_vv_fi,ThetaPU_vv] = carrier_recovery.viterbiViterbi_fxp_mex( ...
                        rx_vv, P.N_pol, P.VV_NTaps, ...
                        testCase.VVFilter_fi, ...
                        pilots_vv, P.BlockLen, stepSize, ...
                        P.UsePilots, P.PilotThreshold, P.CordicIts, testCase.T_vv);
                    

                    cr_vv = test_phase_recovery_stepsize.resolvePhaseAmbiguity( ...
                        double(cr_vv_fi), txRefBits);

                    BER_VV_all(tr,si) = ...
                        test_phase_recovery_stepsize.computeBER( ...
                        cr_vv, txRefBits);

                    %% BPS
                    [cr_bps_fi,ThetaPU_bps] = carrier_recovery.bps_fxp_mex( ...
                        rx_bps, P.BPS_N, P.N_pol, ...
                        P.M, P.BPS_B, ...
                        P.BlockLen, stepSize, ...
                        pilots_bps, P.UsePilots, P.PilotThreshold, P.CordicIts, ...
                        testCase.T_bps);

                    cr_bps = test_phase_recovery_stepsize.resolvePhaseAmbiguity( ...
                        double(cr_bps_fi), txRefBits);

                    BER_BPS_all(tr,si) = ...
                        test_phase_recovery_stepsize.computeBER( ...
                        cr_bps, txRefBits);
                end
            end

            BER_VV  = mean(BER_VV_all,1)
            BER_BPS = mean(BER_BPS_all,1)

            %% Plot
            if P.Plot
                figure('Name','BER vs StepSize', ...
                       'Position',[100 100 800 550], ...
                       'Color','w');

                set(gca,...
                    'FontSize',14,...
                    'LineWidth',1,...
                    'Box','on',...
                    'Color','w');
                
                hold on
                
                if P.PlotTrials
                    for si=1:NS
                        semilogy(StepSizes(si)*ones(P.NTrials,1), ...
                            BER_VV_all(:,si),'o','HandleVisibility','off');
                        semilogy(StepSizes(si)*ones(P.NTrials,1), ...
                            BER_BPS_all(:,si),'s','HandleVisibility','off');
                    end
                end

                semilogy(StepSizes,BER_VV,'s-','LineWidth',2,...
                    'DisplayName','VV fxp');
                semilogy(StepSizes,BER_BPS,'s-','LineWidth',2,...
                    'DisplayName','BPS fxp');

                set(gca,'XScale','log','XTick',StepSizes);
                grid on
                xlabel('Step size [symbols]',...
                    'Interpreter','latex',...
                    'FontSize',16)

                ylabel('BER',...
                    'Interpreter','latex',...
                    'FontSize',16)

                legend('Location','best',...
                    'Interpreter','latex',...
                    'FontSize',14)

                title(sprintf(['BER vs Step Size | QPSK | ' ...
                      'SNR=%.1f dB | LW=%.0f kHz | %d trials'], ...
                      P.SNR_dB, P.LW/1e3, P.NTrials));
            end
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