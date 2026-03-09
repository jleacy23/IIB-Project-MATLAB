classdef test_CarrierRecovery < matlab.unittest.TestCase
    % Tests for carrier recovery – floating-point and fixed-point MEX.
    %
    % Floating-point tests (testQPSK_VV, testQPSK_BPS):
    %   Apply AWGN + phase noise, run carrier recovery, plot
    %   before/after constellations, verify BER < threshold.
    %
    % Fixed-point MEX tests (testQPSK_VV_Fxp16, testQPSK_BPS_Fxp16):
    %   Cast inputs to fi, call the compiled MEX, verify BER < threshold.
    %   No MATLAB vs MEX comparison is performed.
    %
    % MEX compilation
    %   Both MEX binaries are compiled automatically in TestClassSetup
    %   (once per test run, before any test method executes).  A fresh
    %   parameter struct is constructed and its fields are overridden with
    %   the Constant properties defined here so that the compiled binary
    %   exactly matches the test configuration.
    %
    % Prerequisites
    %   - MATLAB Coder and Fixed-Point Designer toolboxes must be licensed.
%   - modem.slicer must be on the MATLAB path (used by carrier_recovery.bps_fxp).
%   - build_carrier_recovery_viterbiViterbi_fxp_mex.m and build_carrier_recovery_bps_fxp_mex.m
    %     must be on the MATLAB path.

    properties (Constant)
        % ---- Signal -------------------------------------------------
        N_pol    = 2
        Ns       = 2^13             % symbols per polarisation
        SpS      = 1                % symbol-rate processing
        BlockLen = 64
        PilotLen = 8
        M        = 4

        % ---- System -------------------------------------------------
        Rs        = 10              % symbol rate [GBd]
        SNR_dB    = 17.5            % [dB]
        Linewidth = 2400e3          % laser linewidth [Hz]
        LW        = 2400e3          % phase-noise linewidth [Hz]

        % ---- Channel (benign) ---------------------------------------
        L       = 80
        D       = 0
        CWL     = 1550
        DGDSpec = 0
        N_pmd   = 1

        % ---- Carrier recovery – shared ------------------------------
        NTaps     = 5
        UsePilots = true
        PilotThreshold = pi/3
        StepSize  = 1               % symbol-by-symbol update (full bandwidth)

        % ---- BPS-specific -------------------------------------------
        B = 64                      % number of blind test phases

        % ---- Fixed-point --------------------------------------------
        FxpConfig = 'fixed16'       % 'fixed16' | 'fixed32'

        % ---- Pass / fail --------------------------------------------
        BER_THRESHOLD = 5e-2
    end

    % ================================================================
    %  One-time setup: seed RNG and compile both MEX binaries
    % ================================================================
    methods (TestClassSetup)

        function seedRng(~)
            rng('shuffle');
        end

        function compileMex(testCase)
            % Build a parameter struct consistent with the test's own
            % Constant properties so the compiled MEX signatures match
            % exactly what the test methods will pass at runtime.
            P.N_pol         = testCase.N_pol;
            P.BlockLen      = testCase.BlockLen;
            P.PilotLen      = testCase.PilotLen;
            P.VV_NTaps      = testCase.NTaps;
            P.BPS_N         = testCase.NTaps;
            P.BPS_B         = testCase.B;
            P.M             = testCase.M;
            P.FxpConfig_VV  = testCase.FxpConfig;
            P.FxpConfig_BPS = testCase.FxpConfig;
            P.StepSize      = testCase.StepSize;
            P.PilotThreshold = testCase.PilotThreshold;

            cfg = coder.config('mex');
            cfg.GenerateReport            = false;
            cfg.SaturateOnIntegerOverflow = false;

            fprintf('  Compiling cr_viterbiViterbi_fxp_mex...\n');
            build_carrier_recovery_viterbiViterbi_fxp_mex(P, cfg);

            fprintf('  Compiling cr_bps_fxp_mex...\n');
            build_carrier_recovery_bps_fxp_mex(P, cfg);
        end

    end

    % ================================================================
    %  Floating-point tests
    % ================================================================
    methods (Test)

        function testQPSK_VV(testCase)
            [rxSym, crSym, BER, ThetaPU] = runScenario_VV(testCase);

            plotBeforeAfter(testCase, rxSym, crSym, 'QPSK | VV', BER);
            plotPhase(testCase, ThetaPU, '4-QAM VV', BER);

            testCase.verifyTrue(all(isfinite(crSym(:))), ...
                'VV carrier-recovery output contains NaN/Inf.');
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD, ...
                sprintf('QPSK VV BER %.2e exceeds threshold.', BER));
        end

        function testQPSK_BPS(testCase)
            [rxSym, crSym, BER, ThetaPU] = runScenario_BPS(testCase);

            plotBeforeAfter(testCase, rxSym, crSym, 'QPSK | BPS', BER);
            plotPhase(testCase, ThetaPU, '4-QAM BPS', BER);

            testCase.verifyTrue(all(isfinite(crSym(:))), ...
                'BPS carrier-recovery output contains NaN/Inf.');
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD, ...
                sprintf('QPSK BPS BER %.2e exceeds threshold.', BER));
        end

    end

    % ================================================================
    %  Fixed-point MEX tests
    % ================================================================
    methods (Test)

        function testQPSK_VV_Fxp16(testCase)
            [rxSym, crSym, BER, ThetaPU] = runScenarioFxp_VV(testCase, testCase.FxpConfig);

            plotBeforeAfter(testCase, rxSym, crSym, ...
                sprintf('QPSK | VV fxp (%s)', testCase.FxpConfig), BER);
            plotPhase(testCase, ThetaPU, '4-QAM VV Fxp', BER);

            testCase.verifyTrue(all(isfinite(crSym(:))), ...
                'VV fxp MEX output contains NaN/Inf.');
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD, ...
                sprintf('QPSK VV fxp16 BER %.2e exceeds threshold.', BER));
        end

        function testQPSK_BPS_Fxp16(testCase)
            [rxSym, crSym, BER, ThetaPU] = runScenarioFxp_BPS(testCase, testCase.FxpConfig);

            plotBeforeAfter(testCase, rxSym, crSym, ...
                sprintf('QPSK | BPS fxp (%s)', testCase.FxpConfig), BER);
            plotPhase(testCase, ThetaPU, '4-QAM BPS Fxp', BER);

            testCase.verifyTrue(all(isfinite(crSym(:))), ...
                'BPS fxp MEX output contains NaN/Inf.');
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD, ...
                sprintf('QPSK BPS fxp16 BER %.2e exceeds threshold.', BER));
        end

    end

    % ================================================================
    %  Private helpers – scenario runners
    % ================================================================
    methods (Access = private)

        % ---- Floating-point VV --------------------------------------
        function [rxSym, crSym, BER, ThetaPU] = runScenario_VV(testCase)
            [symbols, pilots, txRefBits, rxSym] = buildChannel(testCase);

            symEnergy = mean(abs(symbols(:)).^2);
            VVFilter  = carrier_recovery.genVVFilter(testCase.Linewidth, testCase.Rs, ...
                testCase.SNR_dB, symEnergy, testCase.N_pol, testCase.NTaps);

            [crSym, ThetaPU] = carrier_recovery.viterbiViterbi(rxSym, testCase.N_pol, ...
                VVFilter, testCase.BlockLen, testCase.StepSize, ...
                pilots, testCase.UsePilots, testCase.PilotThreshold);

            BER = computeBER(testCase, crSym, txRefBits);
            fprintf('QPSK VV BER = %.2e\n', BER);
        end

        % ---- Floating-point BPS -------------------------------------
        function [rxSym, crSym, BER, ThetaPU] = runScenario_BPS(testCase)
            [~, pilots, txRefBits, rxSym] = buildChannel(testCase);

            [crSym, ThetaPU] = carrier_recovery.bps(rxSym, testCase.NTaps, testCase.N_pol, ...
                testCase.M, testCase.B, testCase.BlockLen, testCase.StepSize, ...
                pilots, testCase.UsePilots, testCase.PilotThreshold);

            BER = computeBER(testCase, crSym, txRefBits);
            fprintf('QPSK BPS BER = %.2e\n', BER);
        end

        % ---- Fixed-point VV MEX -------------------------------------
        function [rxSym, crSym, BER, ThetaPU] = runScenarioFxp_VV(testCase, config)
            T = carrier_recovery.viterbiViterbi_fxp_types(config);

            [symbols, pilots, txRefBits, rxSym] = buildChannel(testCase);

            symEnergy = mean(abs(symbols(:)).^2);
            VVFilter  = carrier_recovery.genVVFilter(testCase.Linewidth, testCase.Rs, ...
                testCase.SNR_dB, symEnergy, testCase.N_pol, testCase.NTaps);

            rxSym_fi    = cast(rxSym,    'like', T.x);
            VVFilter_fi = cast(VVFilter, 'like', T.w);
            pilots_fi   = cast(pilots,   'like', T.x);

            [crSym_fi, ThetaPU_fi] = carrier_recovery.viterbiViterbi_fxp_mex( ...
                rxSym_fi, testCase.N_pol, testCase.NTaps, VVFilter_fi, ...
                pilots_fi, testCase.BlockLen, double(testCase.StepSize), ...
                testCase.UsePilots, testCase.PilotThreshold, T);
            ThetaPU = double(ThetaPU_fi);

            crSym = resolvePhaseAmbiguity(testCase, double(crSym_fi), txRefBits);
            BER   = computeBER(testCase, crSym, txRefBits);
            fprintf('QPSK VV fxp (%s) BER = %.2e\n', config, BER);
        end

        % ---- Fixed-point BPS MEX ------------------------------------
        function [rxSym, crSym, BER, ThetaPU] = runScenarioFxp_BPS(testCase, config)
            T = carrier_recovery.bps_fxp_types(config);

            [~, pilots, txRefBits, rxSym] = buildChannel(testCase);

            rxSym_fi  = cast(rxSym,  'like', T.x);
            pilots_fi = cast(pilots, 'like', T.x);

            [crSym_fi, ThetaPU_fi] = carrier_recovery.bps_fxp_mex( ...
                rxSym_fi, testCase.NTaps, testCase.N_pol, ...
                testCase.M, testCase.B, testCase.BlockLen, double(testCase.StepSize), ...
                pilots_fi, testCase.UsePilots, testCase.PilotThreshold, T);
            ThetaPU = double(ThetaPU_fi);

            crSym = resolvePhaseAmbiguity(testCase, double(crSym_fi), txRefBits);
            BER   = computeBER(testCase, crSym, txRefBits);
            fprintf('QPSK BPS fxp (%s) BER = %.2e\n', config, BER);
        end

        % ---- Shared channel builder ---------------------------------
        function [symbols, pilots, txRefBits, rxSym] = buildChannel(testCase)
            Nbits  = 4 * testCase.Ns;
            txBits = modem.randomBits(Nbits);
            [symbols, pilotSyms, ~, ~] = modem.modulate(txBits);
            pilots = pilotSyms(:, 1);     % single-pol pilot vector for CR
            txRefBits = modem.symbolsToBits(symbols);

            rxSym = channel.add_awgn(symbols, testCase.SNR_dB);
            rxSym = channel.add_phase_noise(rxSym, testCase.Rs, testCase.LW);
        end

        % ---- Phase ambiguity resolution -----------------------------
        function bestSym = resolvePhaseAmbiguity(testCase, crSym, txRefBits)
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

        % ---- BER computation ----------------------------------------
        function BER = computeBER(testCase, crSym, txRefBits)
            decidedSyms = modem.decideSymbols(crSym);
            rxBits      = modem.symbolsToBits(decidedSyms);
            nBits       = min(length(txRefBits), length(rxBits));
            nErrors     = sum(txRefBits(1:nBits) ~= rxBits(1:nBits));
            BER         = nErrors / nBits;
        end

        % ---- Plotting -----------------------------------------------
        function plotBeforeAfter(testCase, rxSym, crSym, titleStr, BER)
            figure('Name', titleStr, 'Position', [100 100 1200 500]);
            for p = 1:testCase.N_pol
                subplot(2, 2, (p-1)*2 + 1);
                plot(real(rxSym(:,p)), imag(rxSym(:,p)), '.', 'MarkerSize', 2);
                grid on; axis equal;
                title(sprintf('Before CR  \u2013  Pol %d', p));
                xlabel('In-Phase'); ylabel('Quadrature');

                subplot(2, 2, (p-1)*2 + 2);
                plot(real(crSym(:,p)), imag(crSym(:,p)), '.', 'MarkerSize', 2);
                grid on; axis equal;
                title(sprintf('After CR  \u2013  Pol %d', p));
                xlabel('In-Phase'); ylabel('Quadrature');
            end
            sgtitle(sprintf('QPSK: AWGN + Phase Noise  |  %s  |  BER = %.2e', ...
                titleStr, BER));
        end

        function plotPhase(testCase, ThetaPU, titleStr, BER)
            figure('Name', [titleStr ' Phase'], 'Position', [100 650 1200 400]);
            for p = 1:testCase.N_pol
                subplot(1, testCase.N_pol, p);
                plot(ThetaPU(:, p));
                grid on;
                xlabel('Symbol index');
                ylabel('\Theta_{PU} [rad]');
                title(sprintf('Phase unwrapped  –  Pol %d', p));
            end
            sgtitle(sprintf('%s Phase Estimate  |  BER = %.2e', titleStr, BER));
        end

    end
end