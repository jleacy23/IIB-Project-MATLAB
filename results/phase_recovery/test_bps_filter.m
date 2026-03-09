classdef test_bps_filter < matlab.unittest.TestCase
%TEST_BPS_FILTER
%   Sweep BPS_N filter length for fixed-point BPS CR.
%   BER vs BPS_N is plotted.

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
        NTrials     = 10

        % Pilot
        PilotLen    = 4
        UsePilots   = true
        PilotThreshold = pi/3

        % Block
        BlockLen    = 128
        StepSize    = 1

        % BPS sweep
        BPS_N_vals  = [1 3 5 7 9 11 15 21]

        % BPS params
        BPS_B       = 64

        % Fixed-point
        FxpConfig   = 'fixed16'
        CordicIts   = 8

        % Plot
        Plot        = true
        PlotTrials   = false
    end

    %% ================================================================
    % Shared Resources
    %% ================================================================
    properties
        T_bps
    end

    %% ================================================================
    % Setup
    %% ================================================================
    methods (TestClassSetup)

        function buildMexAndSetup(testCase)

            P = testCase;

            assert(mod(P.BPS_B,2)==0,'BPS_B must be even.')

            %% Types
            testCase.T_bps = carrier_recovery.bps_fxp_types(P.FxpConfig);

            %% Build config
            Pbuild.N_pol = P.N_pol;
            Pbuild.M     = P.M;

            Pbuild.PilotLen       = P.PilotLen;
            Pbuild.PilotThreshold = P.PilotThreshold;

            Pbuild.BPS_N   = P.BPS_N_vals(1);
            Pbuild.BPS_B   = P.BPS_B;

            Pbuild.BlockLen = P.BlockLen;
            Pbuild.StepSize = P.StepSize;

            Pbuild.FxpConfig_BPS = P.FxpConfig;
            Pbuild.CordicIts     = P.CordicIts;

            cfg = coder.config('mex');
            cfg.GenerateReport = false;

            fprintf('=== Compiling BPS MEX ===\n');

            build_carrier_recovery_bps_fxp_mex(Pbuild,cfg);

            fprintf('=== Compilation done ===\n\n');

        end
    end

    %% ================================================================
    % Test
    %% ================================================================
    methods (Test)

        function test_bps_filter_sweep(testCase)

            P = testCase;

            Nvals = P.BPS_N_vals;
            NN    = length(Nvals);

            BER_all = zeros(P.NTrials,NN);

            for tr = 1:P.NTrials

                %% Channel

                Nbits = 4 * P.Ns;

                txBits = modem.randomBits(Nbits);

                [symbols, pilots, ~, ~] = modem.modulate(txBits);
                txRefBits = modem.symbolsToBits(symbols);

                rxSym = channel.add_awgn(symbols,P.SNR_dB);
                rxSym = channel.add_phase_noise(rxSym,P.Rs,P.LW);

                %% Cast

                rx_bps     = cast(rxSym,'like',testCase.T_bps.x);
                pilots_bps = cast(pilots,'like',testCase.T_bps.x);

                for ni = 1:NN

                    BPS_N = Nvals(ni);

                    fprintf('BPS_N = %d\n',BPS_N);

                    %% BPS

                    [cr_bps_fi,~] = carrier_recovery.bps_fxp_mex( ...
                        rx_bps,...
                        BPS_N,...
                        P.N_pol,...
                        P.M,...
                        P.BPS_B,...
                        P.BlockLen,...
                        P.StepSize,...
                        pilots_bps,...
                        P.UsePilots,...
                        P.PilotThreshold,...
                        P.CordicIts,...
                        testCase.T_bps);

                    cr_bps = ...
                        test_bps_filter.resolvePhaseAmbiguity( ...
                        double(cr_bps_fi),...
                        txRefBits);

                    BER_all(tr,ni) = ...
                        test_bps_filter.computeBER( ...
                        cr_bps,...
                        txRefBits);

                end
            end

            BER = mean(BER_all,1);

            %% Plot

            if P.Plot

                figure( ...
                    'Name','BER vs BPS N',...
                    'Position',[100 100 800 550],...
                    'Color','w');

                hold on

                % Scatter trials
                if P.PlotTrials

                    for ni=1:NN

                        semilogy( ...
                            Nvals(ni)*ones(P.NTrials,1),...
                            BER_all(:,ni),...
                            'o',...
                            'HandleVisibility','off');

                    end

                end

                % Mean curve
                semilogy( ...
                    Nvals,...
                    BER,...
                    'o-',...
                    'LineWidth',2,...
                    'DisplayName','BPS fxp');

                grid on

                set(gca,...
                    'FontSize',14,...
                    'LineWidth',1,...
                    'Box','on',...
                    'Color','w');

                xlabel('BPS filter length $N$',...
                    'Interpreter','latex',...
                    'FontSize',16)

                ylabel('BER',...
                    'Interpreter','latex',...
                    'FontSize',16)

                legend( ...
                    'Location','best',...
                    'Interpreter','latex',...
                    'FontSize',14)

                title(sprintf([ ...
                    'BER vs BPS $N$ | QPSK | ' ...
                    'SNR=%.1f dB | LW=%.0f kHz | %d trials'],...
                    P.SNR_dB,...
                    P.LW/1e3,...
                    P.NTrials),...
                    'Interpreter','latex',...
                    'FontSize',14)

            end

        end
    end

    %% ================================================================
    % Helpers
    %% ================================================================
    methods (Static,Access=private)

        function BER = computeBER(crSym, refBits)

            decidedSyms = ...
                modem.decideSymbols(crSym);

            rxBits = ...
                modem.symbolsToBits(decidedSyms);

            BER = sum(refBits~=rxBits)/length(refBits);

        end

        function bestSym = resolvePhaseAmbiguity( ...
                crSym, refBits)

            bestBER = Inf;
            bestSym = crSym;

            for k=0:3

                rotated = crSym .* exp(-1j*k*pi/2);

                decidedSyms = ...
                    modem.decideSymbols(rotated);

                rxBits = ...
                    modem.symbolsToBits(decidedSyms);

                thisBER = ...
                    sum(refBits~=rxBits)/length(refBits);

                if thisBER < bestBER

                    bestBER = thisBER;
                    bestSym = rotated;

                end

            end

        end

    end

end