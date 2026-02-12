classdef test_CarrierRecovery < matlab.unittest.TestCase
    % Tests for carrier recovery (Viterbi-Viterbi).
    % Applies AWGN + phase noise, runs VV carrier recovery, plots
    % before/after constellations and checks BER.

    properties (Constant)
        N_pol   = 2
        Ns      = 2^15          % symbols per polarisation
        SpS     = 1             % symbol-rate processing (no pulse shaping)

        % System
        Rs      = 32            % [GBd]
        SNR_dB  = 20            % [dB]
        Linewidth = 100e4       % laser linewidth [Hz]

        % Channel (unused impairments set to benign values)
        L       = 80            % fibre length [km]
        D       = 0             % no CD
        CWL     = 1550          % [nm]
        DGDSpec = 0             % no PMD
        N_pmd   = 1
        LW      = 100e4         % phase-noise linewidth [Hz]

        % Carrier recovery
        NTaps   = 15

        % Pass / fail
        BER_THRESHOLD = 5e-2
    end

    methods (TestMethodSetup)
        function seedRng(~)
            rng(42);
        end
    end

    % ================================================================
    methods (Test)

        % -------- 4-QAM -----------------------------
        function testQPSK(testCase)
            M = 4;
            [rxSym, crSym, BER] = runScenario(testCase, M);
            plotBeforeAfter(testCase, rxSym, crSym, ...
                '4-QAM  |  VV', M, BER);

            testCase.verifyTrue(all(isfinite(crSym(:))), ...
                'Carrier-recovery output contains NaN/Inf.');
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD, ...
                sprintf('4-QAM BER %.2e exceeds threshold.', BER));
        end

        % -------- 4-QAM fixed-point (fixed16) ------
        function testQPSK_Fxp16(testCase)
            M = 4;
            [rxSym, crSym, BER] = runScenarioFxp(testCase, M, 'fixed16');
            plotBeforeAfter(testCase, rxSym, crSym, ...
                '4-QAM  |  VV FXP fixed16', M, BER);

            testCase.verifyTrue(all(isfinite(crSym(:))), ...
                'FXP carrier-recovery output contains NaN/Inf.');
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD, ...
                sprintf('4-QAM FXP16 BER %.2e exceeds threshold.', BER));
        end

        % -------- 4-QAM fixed-point (fixed32) ------
        function testQPSK_Fxp32(testCase)
            M = 4;
            [rxSym, crSym, BER] = runScenarioFxp(testCase, M, 'fixed32');
            plotBeforeAfter(testCase, rxSym, crSym, ...
                '4-QAM  |  VV FXP fixed32', M, BER);

            testCase.verifyTrue(all(isfinite(crSym(:))), ...
                'FXP carrier-recovery output contains NaN/Inf.');
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD, ...
                sprintf('4-QAM FXP32 BER %.2e exceeds threshold.', BER));
        end
    end

    % ================================================================
    %  Helpers
    % ================================================================
    methods (Access = private)

        function [rxSym, crSym, BER] = runScenario(testCase, M)
            rng(42);

            % --- Tx ---
            k      = log2(M);
            Nbits  = k * testCase.N_pol * testCase.Ns;
            txBits = qam_randomBits(Nbits);
            symbols = qam_modulate(txBits, M, testCase.N_pol);

            % --- Channel: AWGN + phase noise ---
            rxSym = channel_add_awgn(symbols, testCase.SNR_dB);
            rxSym = channel_add_phase_noise(rxSym, testCase.Rs, testCase.LW);

            % --- Carrier Recovery (Viterbi-Viterbi) ---
            symEnergy = mean(abs(symbols(:)).^2);
            VVFilter  = cr_genVVFilter(testCase.Linewidth, testCase.Rs, ...
                testCase.SNR_dB, symEnergy, testCase.N_pol, testCase.NTaps);
            crSym = cr_viterbiViterbi(rxSym, testCase.N_pol, testCase.NTaps, ...
                VVFilter);

            % --- Demodulate & BER ---
            decidedSyms = qam_decideSymbols(crSym, M, testCase.N_pol);
            rxBits      = qam_symbolsToBits(decidedSyms, M);

            nErrors = sum(txBits ~= rxBits);
            BER     = nErrors / length(txBits);
            fprintf('%d-QAM BER = %.2e  (%d / %d)\n', ...
                     M, BER, nErrors, length(txBits));
        end

        function [rxSym, crSym, BER] = runScenarioFxp(testCase, M, config)
            rng(42);

            % --- Types ---
            T = cr_viterbiViterbi_fxp_types(config);

            % --- Tx ---
            k      = log2(M);
            Nbits  = k * testCase.N_pol * testCase.Ns;
            txBits = qam_randomBits(Nbits);
            symbols = qam_modulate(txBits, M, testCase.N_pol);

            % --- Channel: AWGN + phase noise ---
            rxSym = channel_add_awgn(symbols, testCase.SNR_dB);
            rxSym = channel_add_phase_noise(rxSym, testCase.Rs, testCase.LW);

            % --- VV filter ---
            symEnergy = mean(abs(symbols(:)).^2);
            VVFilter  = cr_genVVFilter(testCase.Linewidth, testCase.Rs, ...
                testCase.SNR_dB, symEnergy, testCase.N_pol, testCase.NTaps);

            % --- Cast inputs to fi (shared by MATLAB and MEX) ---
            rxSym_fi    = cast(rxSym,    'like', T.x);
            VVFilter_fi = cast(VVFilter, 'like', T.w);

            % % --- MATLAB fixed-point carrier recovery ---
            % crSym_ML = cr_viterbiViterbi_fxp(rxSym_fi, testCase.N_pol, ...
            %     testCase.NTaps, VVFilter_fi, T);

            % --- MEX fixed-point carrier recovery ---
            crSym_MEX = cr_viterbiViterbi_fxp_mex(rxSym_fi, testCase.N_pol, ...
                testCase.NTaps, VVFilter_fi, T);

            % % --- Verify MATLAB and MEX are bit-exact ---
            % verifyFxpMexMatchesMatlab(testCase, crSym_ML, crSym_MEX, ...
            %     sprintf('%d-QAM %s', M, config));

            crSym_dbl = double(crSym_MEX);

            % --- Phase ambiguity resolution ---
            %  VV has pi/2 phase ambiguity for QPSK.  Try all four
            %  rotations and pick the one with minimum BER.
            bestBER = Inf;
            bestSym = crSym_dbl;
            for kk = 0:3
                rotated     = crSym_dbl .* exp(-1j * kk * pi/2);
                decidedSyms = qam_decideSymbols(rotated, M, testCase.N_pol);
                rxBits      = qam_symbolsToBits(decidedSyms, M);
                nErrors     = sum(txBits ~= rxBits);
                thisBER     = nErrors / length(txBits);
                if thisBER < bestBER
                    bestBER = thisBER;
                    bestSym = rotated;
                end
            end

            BER   = bestBER;
            crSym = bestSym;
            fprintf('%d-QAM FXP (%s) BER = %.2e\n', M, config, BER);
        end

        function verifyFxpMexMatchesMatlab(testCase, eqML, eqMEX, tag)
            %VERIFYFXPMEXMATCHESMATLAB  Check MATLAB fxp and MEX are bit-exact.
            mlDbl  = double(eqML);
            mexDbl = double(eqMEX);
            maxErr = max(abs(mlDbl(:) - mexDbl(:)));
            fprintf('  [%s] max |MATLAB-MEX| = %g\n', tag, maxErr);
            testCase.verifyEqual(mexDbl, mlDbl, ...
                sprintf('[%s] MEX output differs from MATLAB fxp.', tag));
        end

        function plotBeforeAfter(testCase, rxSym, crSym, titleStr, M, BER)
            figure('Name', titleStr, 'Position', [100 100 1200 500]);
            for p = 1:testCase.N_pol
                % Before carrier recovery
                subplot(2, 2, (p-1)*2 + 1);
                plot(real(rxSym(:,p)), imag(rxSym(:,p)), '.', ...
                     'MarkerSize', 2);
                grid on; axis equal;
                title(sprintf('Before CR  –  Pol %d', p));
                xlabel('In-Phase'); ylabel('Quadrature');

                % After carrier recovery
                subplot(2, 2, (p-1)*2 + 2);
                plot(real(crSym(:,p)), imag(crSym(:,p)), '.', ...
                     'MarkerSize', 2);
                grid on; axis equal;
                title(sprintf('After CR  –  Pol %d', p));
                xlabel('In-Phase'); ylabel('Quadrature');
            end
            sgtitle(sprintf('%d-QAM: AWGN + Phase Noise  |  %s  |  BER = %.2e', ...
                             M, titleStr, BER));
        end
    end
end
