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
    % Prerequisites for fxp tests
    %   - cr_viterbiViterbi_fxp_mex and cr_bps_fxp_mex must be compiled.
    %     Run build_all_mex() before executing fxp tests.
    %   - qam_slicer.m must be on the MATLAB path (used by cr_bps_fxp).

    properties (Constant)
        % ---- Signal -------------------------------------------------
        N_pol    = 2
        Ns       = 2^15             % symbols per polarisation
        SpS      = 1                % symbol-rate processing
        BlockLen = 64
        PilotLen = 8

        % ---- System -------------------------------------------------
        Rs        = 32              % symbol rate [GBd]
        SNR_dB    = 20              % [dB]
        Linewidth = 1000e3          % laser linewidth [Hz]
        LW        = 1000e3          % phase-noise linewidth [Hz]

        % ---- Channel (benign) ---------------------------------------
        L       = 80
        D       = 0
        CWL     = 1550
        DGDSpec = 0
        N_pmd   = 1

        % ---- Carrier recovery – shared ------------------------------
        NTaps      = 5
        UsePilots  = true
        BlockBased = false

        % ---- BPS-specific -------------------------------------------
        B = 64                      % number of blind test phases

        % ---- Fixed-point --------------------------------------------
        FxpConfig = 'fixed16'       % 'fixed16' | 'fixed32'

        % ---- Pass / fail --------------------------------------------
        BER_THRESHOLD = 5e-2
    end

    methods (TestClassSetup)
        function seedRng(~)
            rng('shuffle');
        end
    end

    % ================================================================
    %  Floating-point tests
    % ================================================================
    methods (Test)

        function testQPSK_VV(testCase)
            M = 4;
            [rxSym, crSym, BER, ThetaPU] = runScenario_VV(testCase, M);

            plotBeforeAfter(testCase, rxSym, crSym, '4-QAM | VV', M, BER);
            plotPhase(testCase, ThetaPU, '4-QAM VV', BER);

            testCase.verifyTrue(all(isfinite(crSym(:))), ...
                'VV carrier-recovery output contains NaN/Inf.');
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD, ...
                sprintf('4-QAM VV BER %.2e exceeds threshold.', BER));
        end

        function testQPSK_BPS(testCase)
            M = 4;
            [rxSym, crSym, BER, ThetaPU] = runScenario_BPS(testCase, M);

            plotBeforeAfter(testCase, rxSym, crSym, '4-QAM | BPS', M, BER);
            plotPhase(testCase, ThetaPU, '4-QAM BPS', BER);

            testCase.verifyTrue(all(isfinite(crSym(:))), ...
                'BPS carrier-recovery output contains NaN/Inf.');
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD, ...
                sprintf('4-QAM BPS BER %.2e exceeds threshold.', BER));
        end

    end

    % ================================================================
    %  Fixed-point MEX tests
    % ================================================================
    methods (Test)

        function testQPSK_VV_Fxp16(testCase)
            M = 4;
            [rxSym, crSym, BER, ThetaPU] = runScenarioFxp_VV(testCase, M, testCase.FxpConfig);

            plotBeforeAfter(testCase, rxSym, crSym, ...
                sprintf('4-QAM | VV fxp (%s)', testCase.FxpConfig), M, BER);
            plotPhase(testCase, ThetaPU, '4-QAM VV Fxp', BER);
            testCase.verifyTrue(all(isfinite(crSym(:))), ...
                'VV fxp MEX output contains NaN/Inf.');
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD, ...
                sprintf('4-QAM VV fxp16 BER %.2e exceeds threshold.', BER));
        end

        function testQPSK_BPS_Fxp16(testCase)
            M = 4;
            [rxSym, crSym, BER, ThetaPU] = runScenarioFxp_BPS(testCase, M, testCase.FxpConfig);

            plotBeforeAfter(testCase, rxSym, crSym, ...
                sprintf('4-QAM | BPS fxp (%s)', testCase.FxpConfig), M, BER);
            plotPhase(testCase, ThetaPU, '4-QAM BPS Fxp', BER);

            testCase.verifyTrue(all(isfinite(crSym(:))), ...
                'BPS fxp MEX output contains NaN/Inf.');
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD, ...
                sprintf('4-QAM BPS fxp16 BER %.2e exceeds threshold.', BER));
        end

    end

    % ================================================================
    %  Private helpers – scenario runners
    % ================================================================
    methods (Access = private)

        % ---- Floating-point VV --------------------------------------
        function [rxSym, crSym, BER, ThetaPU] = runScenario_VV(testCase, M)
            [symbols, pilots, txBits, rxSym] = buildChannel(testCase, M);

            symEnergy = mean(abs(symbols(:)).^2);
            VVFilter  = cr_genVVFilter(testCase.Linewidth, testCase.Rs, ...
                testCase.SNR_dB, symEnergy, testCase.N_pol, testCase.NTaps);

            [crSym, ThetaPU] = cr_viterbiViterbi(rxSym, testCase.N_pol, ...
                VVFilter, testCase.BlockLen, pilots, ...
                testCase.UsePilots, testCase.BlockBased);

            BER = computeBER(testCase, crSym, txBits, M);
            fprintf('4-QAM VV BER = %.2e\n', BER);
        end

        % ---- Floating-point BPS -------------------------------------
        function [rxSym, crSym, BER, ThetaPU] = runScenario_BPS(testCase, M)
            [~, pilots, txBits, rxSym] = buildChannel(testCase, M);

            [crSym, ThetaPU] = cr_bps(rxSym, testCase.NTaps, testCase.N_pol, ...
                M, testCase.B, testCase.BlockLen, pilots, ...
                testCase.UsePilots, testCase.BlockBased);

            BER = computeBER(testCase, crSym, txBits, M);
            fprintf('4-QAM BPS BER = %.2e\n', BER);
        end

        % ---- Fixed-point VV MEX -------------------------------------
        function [rxSym, crSym, BER, ThetaPU] = runScenarioFxp_VV(testCase, M, config)
            T = cr_viterbiViterbi_fxp_types(config);

            [symbols, pilots, txBits, rxSym] = buildChannel(testCase, M);

            symEnergy = mean(abs(symbols(:)).^2);
            VVFilter  = cr_genVVFilter(testCase.Linewidth, testCase.Rs, ...
                testCase.SNR_dB, symEnergy, testCase.N_pol, testCase.NTaps);

            % Cast inputs to fi before handing to MEX
            rxSym_fi    = cast(rxSym,    'like', T.x);
            VVFilter_fi = cast(VVFilter, 'like', T.w);
            pilots_fi   = cast(pilots,   'like', T.x);

            [crSym_fi, ThetaPU_fi] = cr_viterbiViterbi_fxp_mex( ...
                rxSym_fi, testCase.N_pol, testCase.NTaps, VVFilter_fi, ...
                pilots_fi, testCase.BlockLen, ...
                testCase.UsePilots, testCase.BlockBased, T);
            ThetaPU = double(ThetaPU_fi);

            crSym = resolvePhaseAmbiguity(testCase, double(crSym_fi), txBits, M);
            BER   = computeBER(testCase, crSym, txBits, M);
            fprintf('4-QAM VV fxp (%s) BER = %.2e\n', config, BER);
        end

        % ---- Fixed-point BPS MEX ------------------------------------
        function [rxSym, crSym, BER, ThetaPU] = runScenarioFxp_BPS(testCase, M, config)
            T = cr_bps_fxp_types(config);

            [~, pilots, txBits, rxSym] = buildChannel(testCase, M);

            % Cast inputs to fi before handing to MEX
            rxSym_fi  = cast(rxSym,  'like', T.x);
            pilots_fi = cast(pilots, 'like', T.x);

            [crSym_fi, ThetaPU_fi] = cr_bps_fxp_mex( ...
                rxSym_fi, testCase.NTaps, testCase.N_pol, ...
                M, testCase.B, testCase.BlockLen, ...
                pilots_fi, testCase.UsePilots, testCase.BlockBased, T);
            ThetaPU = double(ThetaPU_fi);

            crSym = resolvePhaseAmbiguity(testCase, double(crSym_fi), txBits, M);
            BER   = computeBER(testCase, crSym, txBits, M);
            fprintf('4-QAM BPS fxp (%s) BER = %.2e\n', config, BER);
        end

        % ---- Shared channel builder ---------------------------------
        function [symbols, pilots, txBits, rxSym] = buildChannel(testCase, M)
            k      = log2(M);
            Nbits  = k * testCase.N_pol * testCase.Ns;
            txBits = qam_randomBits(Nbits, testCase.BlockLen, testCase.PilotLen, M);
            [symbols, pilots] = qam_modulate(txBits, M, testCase.N_pol, testCase.PilotLen);

            rxSym = channel_add_awgn(symbols, testCase.SNR_dB);
            rxSym = channel_add_phase_noise(rxSym, testCase.Rs, testCase.LW);
        end

        % ---- Phase ambiguity resolution -----------------------------
        function bestSym = resolvePhaseAmbiguity(testCase, crSym, txBits, M)
            % Pilots prevent cycle slips during the sequence but do not
            % resolve the absolute starting quadrant ambiguity.  Try all
            % four pi/2 rotations and keep the one with the lowest BER.
            bestBER = Inf;
            bestSym = crSym;
            for k = 0:3
                rotated     = crSym .* exp(-1j * k * pi/2);
                decidedSyms = qam_decideSymbols(rotated, M, testCase.N_pol);
                rxBits      = qam_symbolsToBits(decidedSyms, M);
                thisBER     = sum(txBits ~= rxBits) / length(txBits);
                if thisBER < bestBER
                    bestBER = thisBER;
                    bestSym = rotated;
                end
            end
        end

        % ---- BER computation ----------------------------------------
        function BER = computeBER(testCase, crSym, txBits, M)
            decidedSyms = qam_decideSymbols(crSym, M, testCase.N_pol);
            rxBits      = qam_symbolsToBits(decidedSyms, M);
            nErrors     = sum(txBits ~= rxBits);
            BER         = nErrors / length(txBits);
        end

        % ---- Plotting -----------------------------------------------
        function plotBeforeAfter(testCase, rxSym, crSym, titleStr, M, BER)
            figure('Name', titleStr, 'Position', [100 100 1200 500]);
            for p = 1:testCase.N_pol
                subplot(2, 2, (p-1)*2 + 1);
                plot(real(rxSym(:,p)), imag(rxSym(:,p)), '.', 'MarkerSize', 2);
                grid on; axis equal;
                title(sprintf('Before CR  –  Pol %d', p));
                xlabel('In-Phase'); ylabel('Quadrature');

                subplot(2, 2, (p-1)*2 + 2);
                plot(real(crSym(:,p)), imag(crSym(:,p)), '.', 'MarkerSize', 2);
                grid on; axis equal;
                title(sprintf('After CR  –  Pol %d', p));
                xlabel('In-Phase'); ylabel('Quadrature');
            end
            sgtitle(sprintf('%d-QAM: AWGN + Phase Noise  |  %s  |  BER = %.2e', ...
                M, titleStr, BER));
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
