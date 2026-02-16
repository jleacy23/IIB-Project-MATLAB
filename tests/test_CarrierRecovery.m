classdef test_CarrierRecovery < matlab.unittest.TestCase
    % Tests for carrier recovery (Viterbi-Viterbi).
    % Applies AWGN + phase noise, runs VV carrier recovery, plots
    % before/after constellations and checks BER.

    properties (Constant)
        N_pol   = 2
        Ns      = 2^16          % symbols per polarisation
        SpS     = 1             % symbol-rate processing (no pulse shaping)

        % System
        Rs      = 2            % [GBd]
        SNR_dB  = 15            % [dB]
        Linewidth = 200e4       % laser linewidth [Hz]

        % Channel (unused impairments set to benign values)
        L       = 80            % fibre length [km]
        D       = 0             % no CD
        CWL     = 1550          % [nm]
        DGDSpec = 0             % no PMD
        N_pmd   = 1
        LW      = 200e4         % phase-noise linewidth [Hz]

        % Carrier recovery
        NTaps   = 15
        CR_BlockLen    = 256    % block length L for pilot-aided CS correction
        CR_NPilots     = 8      % pilot symbols P per block
        CR_CSThreshold = pi/2   % cycle-slip detection threshold [rad]

        % Pass / fail
        BER_THRESHOLD = 5e-2
    end

    methods (TestClassSetup)
        function seedRng(~)
            rng('shuffle');
        end
    end

    % ================================================================
    methods (Test)

        % -------- 4-QAM (no pilots) -----------------
        function testQPSK(testCase)
            M = 4;
            [rxSym, crSym, BER] = runScenario(testCase, M, false);
            % plotBeforeAfter(testCase, rxSym, crSym, ...
            %     '4-QAM  |  VV (no pilots)', M, BER);

            testCase.verifyTrue(all(isfinite(crSym(:))), ...
                'Carrier-recovery output contains NaN/Inf.');
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD, ...
                sprintf('4-QAM BER %.2e exceeds threshold.', BER));
        end

        % -------- 4-QAM (with pilots) ----------------
        function testQPSK_Pilots(testCase)
            M = 4;
            [rxSym, crSym, BER] = runScenario(testCase, M, true);
            % plotBeforeAfter(testCase, rxSym, crSym, ...
            %     '4-QAM  |  VV (pilots)', M, BER);

            testCase.verifyTrue(all(isfinite(crSym(:))), ...
                'Carrier-recovery output contains NaN/Inf.');
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD, ...
                sprintf('4-QAM Pilots BER %.2e exceeds threshold.', BER));
        end

        % -------- 4-QAM fixed-point (fixed16, no pilots) ------
        function testQPSK_Fxp16(testCase)
            M = 4;
            [rxSym, crSym, BER] = runScenarioFxp(testCase, M, 'fixed16', false);
            % plotBeforeAfter(testCase, rxSym, crSym, ...
            %     '4-QAM  |  VV FXP fixed16 (no pilots)', M, BER);

            testCase.verifyTrue(all(isfinite(crSym(:))), ...
                'FXP carrier-recovery output contains NaN/Inf.');
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD, ...
                sprintf('4-QAM FXP16 BER %.2e exceeds threshold.', BER));
        end

        % -------- 4-QAM fixed-point (fixed16, with pilots) ------
        function testQPSK_Fxp16_Pilots(testCase)
            M = 4;
            [rxSym, crSym, BER] = runScenarioFxp(testCase, M, 'fixed16', true);
            % plotBeforeAfter(testCase, rxSym, crSym, ...
            %     '4-QAM  |  VV FXP fixed16 (pilots)', M, BER);

            testCase.verifyTrue(all(isfinite(crSym(:))), ...
                'FXP carrier-recovery output contains NaN/Inf.');
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD, ...
                sprintf('4-QAM FXP16 Pilots BER %.2e exceeds threshold.', BER));
        end

        % -------- 4-QAM frequency recovery --------------------------
        function testQPSK_FreqRecovery(testCase)
            M      = 4;
            DeltaF = 200;       % 20 MHz frequency offset

            % --- Tx ---
            k      = log2(M);
            Nbits  = k * testCase.N_pol * testCase.Ns;
            txBits = qam_randomBits(Nbits);
            symbols = qam_modulate(txBits, M, testCase.N_pol);

            % --- Channel: AWGN + frequency shift ---
            rxSym = channel_add_awgn(symbols, testCase.SNR_dB);
            rxSym = channel_lo_freq_shift(rxSym, DeltaF, ...
                testCase.Rs, testCase.SpS);

            % --- Frequency recovery ---
            [crSym, estFreqOffset] = cr_freq_recovery(rxSym, testCase.Rs);

            % --- Demodulate & BER (resolve pi/2 ambiguity) ---
            bestBER = Inf;
            bestSym = crSym;
            for kk = 0:3
                rotated     = crSym .* exp(-1j * kk * pi/2);
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
            fprintf('4-QAM Freq Recovery BER = %.2e  (%d / %d)\n', ...
                     BER, round(BER*length(txBits)), length(txBits));

            % --- Plot ---
            plotBeforeAfter(testCase, rxSym, crSym, ...
                '4-QAM  |  Freq Recovery', M, BER);

            % --- Verify ---
            testCase.verifyTrue(all(isfinite(crSym(:))), ...
                'Freq-recovery output contains NaN/Inf.');
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD, ...
                sprintf('4-QAM Freq Recovery BER %.2e exceeds threshold.', BER));

            % --- Verify estimated frequency offset ---
            freqErr = abs(estFreqOffset - DeltaF);
            fprintf('  Estimated offset = %.6e MHz, true = %.6e MHz, error = %.6e MHz\n', ...
                     estFreqOffset, DeltaF, freqErr);
            freqTol = testCase.Rs * 1e9 / testCase.Ns;  % spectral resolution
            testCase.verifyLessThan(freqErr, freqTol, ...
                sprintf('Freq offset estimate error %.2e MHz exceeds resolution %.2e Hz.', ...
                         freqErr, freqTol));
        end
    end

    % ================================================================
    %  Helpers
    % ================================================================
    methods (Access = private)

        function [rxSym, crSym, BER] = runScenario(testCase, M, usePilots)

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

            % --- Form pilots: first P symbols of every block of L ---
            Pval     = testCase.CR_NPilots;
            Lval     = testCase.CR_BlockLen;
            NBlocks  = ceil(testCase.Ns / Lval);
            Pilots   = zeros(NBlocks * Pval, testCase.N_pol);
            for b = 1:NBlocks
                srcIdx = (b-1)*Lval + (1:Pval);
                dstIdx = (b-1)*Pval + (1:Pval);
                Pilots(dstIdx, :) = symbols(srcIdx, :);
            end

            crSym = cr_viterbiViterbi(rxSym, testCase.N_pol, testCase.NTaps, ...
                VVFilter, Pilots, Pval, Lval, testCase.CR_CSThreshold, usePilots);

            % --- Demodulate & BER ---
            decidedSyms = qam_decideSymbols(crSym, M, testCase.N_pol);
            rxBits      = qam_symbolsToBits(decidedSyms, M);

            nErrors = sum(txBits ~= rxBits);
            BER     = nErrors / length(txBits);
            fprintf('%d-QAM BER = %.2e  (%d / %d)\n', ...
                     M, BER, nErrors, length(txBits));
        end

        function [rxSym, crSym, BER] = runScenarioFxp(testCase, M, config, usePilots)

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

            % --- Form pilots: first P symbols of every block of L ---
            Pval     = testCase.CR_NPilots;
            Lval     = testCase.CR_BlockLen;
            NBlocks  = ceil(testCase.Ns / Lval);
            Pilots   = zeros(NBlocks * Pval, testCase.N_pol);
            for b = 1:NBlocks
                srcIdx = (b-1)*Lval + (1:Pval);
                dstIdx = (b-1)*Pval + (1:Pval);
                Pilots(dstIdx, :) = symbols(srcIdx, :);
            end

            % --- Cast inputs to fi (shared by MATLAB and MEX) ---
            rxSym_fi    = cast(rxSym,    'like', T.x);
            VVFilter_fi = cast(VVFilter, 'like', T.w);
            Pilots_fi   = cast(Pilots,   'like', T.x);

            % --- MEX fixed-point carrier recovery ---
            crSym_MEX = cr_viterbiViterbi_fxp_mex(rxSym_fi, testCase.N_pol, ...
                testCase.NTaps, VVFilter_fi, Pilots_fi, Pval, Lval, ...
                testCase.CR_CSThreshold, usePilots, T);

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
