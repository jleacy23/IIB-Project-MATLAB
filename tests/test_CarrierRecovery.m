classdef test_CarrierRecovery < matlab.unittest.TestCase
    % Tests for carrier recovery – floating-point and fixed-point MEX.
    %
    % Floating-point tests (testQPSK_VV, testQPSK_BPS, testQPSK_PilotsOnly):
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
        Ns       = 2^12             % symbols per polarisation
        SpS      = 1                % symbol-rate processing
        BlockLen = 64
        PilotLen = 8
        M        = 4

        % ---- System -------------------------------------------------
        Rs        = 30.5              % symbol rate [GBd]
        SNR_dB    = 20            % [dB]
        Linewidth = 100e3          % laser linewidth [Hz]
        LW        = 100e3          % phase-noise linewidth [Hz]
        frequency_offset = 0       %[MHz]

        % ---- Channel (benign) ---------------------------------------
        L       = 80
        D       = 0
        CWL     = 1550
        DGDSpec = 0
        N_pmd   = 1

        % ---- Carrier recovery – shared ------------------------------
        NTaps     = 5
        PilotThreshold = 5 * pi / 9
        StepSize  = 1               % symbol-by-symbol update (full bandwidth)

        % ---- BPS-specific -------------------------------------------
        B = 64                      % number of blind test phases

        % ---- Fixed-point --------------------------------------------
        CordicIts = 16              % CORDIC iterations for fxp builds
        FxpConfig = 'fixed16'       % 'fixed16' | 'fixed32'

        % ---- Pass / fail --------------------------------------------
        BER_THRESHOLD = 5e-2

        Rebuild = false
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
            P.StepSize       = testCase.StepSize;
            P.PilotThreshold = testCase.PilotThreshold;
            P.CordicIts      = testCase.CordicIts;

            cfg = coder.config('mex');
            cfg.GenerateReport            = false;
            cfg.SaturateOnIntegerOverflow = false;

            if testCase.Rebuild
                fprintf('  Compiling cr_viterbiViterbi_fxp_mex...\n');
                build_carrier_recovery_viterbiViterbi_fxp_mex(P, cfg);

                fprintf('  Compiling cr_bps_fxp_mex...\n');
                build_carrier_recovery_bps_fxp_mex(P, cfg);
            end
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

        function testQPSK_PilotsOnly(testCase)
            [rxSym, crSym, BER, ThetaPU] = runScenario_PilotsOnly(testCase);

            plotBeforeAfter(testCase, rxSym, crSym, 'QPSK | Pilots-Only', BER);
            plotPhase(testCase, ThetaPU, '4-QAM Pilots-Only', BER);

            testCase.verifyTrue(all(isfinite(crSym(:))), ...
                'Pilots-only carrier-recovery output contains NaN/Inf.');
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD, ...
                sprintf('QPSK Pilots-Only BER %.2e exceeds threshold.', BER));
        end

    end

    % ================================================================
    %  Fixed-point MEX tests
    % ================================================================
    methods (Test)

        function testQPSK_VV_Fxp16(testCase)
            [rxSym, crSym, BER, ThetaPU, crSymTheta, BERTheta] = ...
                runScenarioFxp_VV(testCase, testCase.FxpConfig);

            plotBeforeAfter(testCase, rxSym, crSym, ...
                sprintf('QPSK | VV fxp (%s)', testCase.FxpConfig), BER);
            plotBeforeAfter(testCase, rxSym, crSymTheta, ...
                sprintf('QPSK | VV float e^{-j\\theta} (%s)', testCase.FxpConfig), BERTheta);
            plotPhase(testCase, ThetaPU, '4-QAM VV Fxp', BER);

            testCase.verifyTrue(all(isfinite(crSym(:))), ...
                'VV fxp MEX output contains NaN/Inf.');
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD, ...
                sprintf('QPSK VV fxp16 BER %.2e exceeds threshold.', BER));
        end

        function testQPSK_BPS_Fxp16(testCase)
            [rxSym, crSym, BER, ThetaPU, crSymTheta, BERTheta] = ...
                runScenarioFxp_BPS(testCase, testCase.FxpConfig);

            plotBeforeAfter(testCase, rxSym, crSym, ...
                sprintf('QPSK | BPS fxp (%s)', testCase.FxpConfig), BER);
            plotBeforeAfter(testCase, rxSym, crSymTheta, ...
                sprintf('QPSK | BPS float e^{-j\\theta} (%s)', testCase.FxpConfig), BERTheta);
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

        % ---- Floating-point Pilots-Only -----------------------------
        function [rxSym, crSym, BER, ThetaPU] = runScenario_PilotsOnly(testCase)
            [~, pilots, txRefBits, rxSym] = buildChannel(testCase);

            [crSym, ThetaPU] = carrier_recovery.pilots_only(rxSym, testCase.N_pol, ...
                testCase.BlockLen, pilots);

            BER = computeBER(testCase, crSym, txRefBits);
            fprintf('QPSK Pilots-Only BER = %.2e\n', BER);
        end

        % ---- Floating-point VV --------------------------------------
        function [rxSym, crSym, BER, ThetaPU] = runScenario_VV(testCase)
            [symbols, pilots, txRefBits, rxSym] = buildChannel(testCase);

            symEnergy = mean(abs(symbols(:)).^2);
            VVFilter  = carrier_recovery.genVVFilter(testCase.Linewidth, testCase.Rs, ...
                testCase.SNR_dB, symEnergy, testCase.N_pol, testCase.NTaps);

            [crSym, ThetaPU] = carrier_recovery.viterbiViterbi(rxSym, testCase.N_pol, ...
                VVFilter, testCase.BlockLen, testCase.StepSize, ...
                pilots, testCase.PilotThreshold);

            BER = computeBER(testCase, crSym, txRefBits);
            fprintf('QPSK VV BER = %.2e\n', BER);
        end

        % ---- Floating-point BPS -------------------------------------
        function [rxSym, crSym, BER, ThetaPU] = runScenario_BPS(testCase)
            [~, pilots, txRefBits, rxSym] = buildChannel(testCase);

            [crSym, ThetaPU] = carrier_recovery.bps(rxSym, testCase.NTaps, testCase.N_pol, ...
                testCase.M, testCase.B, testCase.BlockLen, testCase.StepSize, ...
                pilots, testCase.PilotThreshold);

            BER = computeBER(testCase, crSym, txRefBits);
            fprintf('QPSK BPS BER = %.2e\n', BER);
        end

        % ---- Fixed-point VV MEX -------------------------------------
        function [rxSym, crSym, BER, ThetaPU, crSymTheta, BERTheta] = runScenarioFxp_VV(testCase, config)
            T = carrier_recovery.fxp_types(config);

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
                testCase.PilotThreshold, double(testCase.CordicIts), T);
            ThetaPU = double(ThetaPU_fi);

            crSym = resolvePhaseAmbiguity(testCase, double(crSym_fi), txRefBits);
            BER   = computeBER(testCase, crSym, txRefBits);

            crSymTheta = applyFloatPhaseCorrection(testCase, rxSym, ThetaPU);
            crSymTheta = resolvePhaseAmbiguity(testCase, crSymTheta, txRefBits);
            BERTheta   = computeBER(testCase, crSymTheta, txRefBits);

            fprintf('QPSK VV fxp (%s) BER = %.2e | float e^{-jtheta} BER = %.2e\n', ...
                config, BER, BERTheta);
        end

        % ---- Fixed-point BPS MEX ------------------------------------
        function [rxSym, crSym, BER, ThetaPU, crSymTheta, BERTheta] = runScenarioFxp_BPS(testCase, config)
            T = carrier_recovery.fxp_types(config);

            [~, pilots, txRefBits, rxSym] = buildChannel(testCase);

            rxSym_fi  = cast(rxSym,  'like', T.x);
            pilots_fi = cast(pilots, 'like', T.x);

            [crSym_fi, ThetaPU_fi] = carrier_recovery.bps_fxp_mex( ...
                rxSym_fi, testCase.NTaps, testCase.N_pol, ...
                testCase.M, testCase.B, testCase.BlockLen, double(testCase.StepSize), ...
                pilots_fi, testCase.PilotThreshold, double(testCase.CordicIts), T);
            ThetaPU = double(ThetaPU_fi);

            crSym = resolvePhaseAmbiguity(testCase, double(crSym_fi), txRefBits);
            BER   = computeBER(testCase, crSym, txRefBits);

            crSymTheta = applyFloatPhaseCorrection(testCase, rxSym, ThetaPU);
            crSymTheta = resolvePhaseAmbiguity(testCase, crSymTheta, txRefBits);
            BERTheta   = computeBER(testCase, crSymTheta, txRefBits);

            fprintf('QPSK BPS fxp (%s) BER = %.2e | float e^{-jtheta} BER = %.2e\n', ...
                config, BER, BERTheta);
        end

        % ---- Shared channel builder ---------------------------------
        function [symbols, pilots, txRefBits, rxSym] = buildChannel(testCase)
            Nbits  = 4 * testCase.Ns;
            txBits = modem.randomBits(Nbits);
            [symbols, pilotSyms, ~, ~] = modem.modulate(txBits);

            % Build per-CR-block pilot matrix: pilots(b,:) must be the known
            % TX pilot at signal position (b-1)*BlockLen+1.
            %
            % CPON places one pilot at the start of every 32-symbol block, so
            % the pilot at signal position p is pilotSyms(cponBlock,:) where
            %   cponBlock = floor((posInSubframe-1)/32) + 1
            %   posInSubframe = mod(p-1, 3712) + 1   (pattern repeats each SF)
            %
            % When BlockLen > 32 (e.g. 64), each CR block spans multiple CPON
            % blocks, so we step through CPON pilot indices by BlockLen/32.
            % BuildChannel handles any BlockLen that is a multiple of 32.
            CPON_BLOCK_LEN = 32;
            CPON_SF_SYMS   = 3712;
            Nsym    = size(symbols, 1);
            NBlocks = ceil(Nsym / testCase.BlockLen);
            pilots  = zeros(NBlocks, testCase.N_pol);
            for b = 1:NBlocks
                pos       = (b - 1) * testCase.BlockLen + 1;
                posInSf   = mod(pos - 1, CPON_SF_SYMS) + 1;
                cponBlock = floor((posInSf - 1) / CPON_BLOCK_LEN) + 1;
                pilots(b, :) = pilotSyms(cponBlock, :);
            end

            txRefBits = modem.symbolsToBits(symbols);
            rxSym = channel.add_awgn(symbols, testCase.SNR_dB);
            rxSym = channel.lo_freq_shift(rxSym, testCase.frequency_offset, testCase.Rs, testCase.SpS)
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

        % ---- Diagnostic: apply ThetaPU in floating-point -----------
        function y = applyFloatPhaseCorrection(~, rxSym, ThetaPU)
            y = rxSym;

            nRows = min(size(rxSym, 1), size(ThetaPU, 1));
            nPol  = min(size(rxSym, 2), size(ThetaPU, 2));
            if nRows == 0 || nPol == 0
                return;
            end

            rot = exp(-1j * ThetaPU(1:nRows, 1:nPol));
            y(1:nRows, 1:nPol) = rxSym(1:nRows, 1:nPol) .* rot;
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