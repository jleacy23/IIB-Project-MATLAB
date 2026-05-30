classdef test_CarrierRecovery < matlab.unittest.TestCase
%TEST_CARRIERRECOVERY  Fixed-point carrier-recovery tests (VV and pilots-only).
%
%   Exercises the CORDIC-based fixed-point carrier-recovery algorithms
%       carrier_recovery.viterbiViterbi_fxp
%       carrier_recovery.pilots_only_fxp
%   end-to-end through the pipeline
%
%       modem.modulate  ->  AWGN  ->  1 MHz phase noise  ->  modem.normalise
%                       ->  fixed-point CR (CORDIC)  ->  BER + constellation
%
%   The received signal is normalised into the unit box [-1,1] (modem.normalise)
%   before being cast to the fixed-point input type, so the wordlength is spent
%   on the signal of interest rather than on noise outliers.
%
%   Fixed-point configuration
%     Uniform type WL = IntBits + FL, fraction length FL.  Following the
%     swept-precision convention used in results/.../bit_width_full, the
%     number of CORDIC iterations equals the fractional precision FL (one
%     iteration per fraction bit) — angular resolution ~atan(2^-FL).
%
%   Execution mode (UseMex)
%     UseMex = true  -> builds and calls the compiled *_fxp_mex binaries
%                       (requires MATLAB Coder + Fixed-Point Designer).
%     UseMex = false -> calls the *_fxp functions directly (interpreted fi;
%                       identical arithmetic, no toolbox/build needed).
%
%   Each test verifies the post-correction BER is below threshold and plots
%   the recovered constellation (received vs recovered, per polarisation).

    properties (Constant)
        % ---- System -------------------------------------------------
        Rs             = 30.5            % symbol rate [GBd]
        N_pol          = 2
        NTaps          = 10             % VV one-sided window half-length
        BlockLen       = 32             % CPON block = CR block
        StepSize       = 32             % one phase update per block
        PilotThreshold = 5 * pi / 9
        M              = 4              % QPSK

        % ---- Channel ------------------------------------------------
        SNR_dB    = 20
        LW_Hz     = 1e6                % 1 MHz phase-noise linewidth
        NormPct   = 99                 % modem.normalise percentile

        % ---- Signal length ------------------------------------------
        NumSubframes = 4               % 1 subframe = 3712 symbols / pol

        % ---- Fixed-point: WL = IntBits + FL, CORDIC its = FL ---------
        IntBits = 16
        FL      = 12

        % ---- Pass / fail --------------------------------------------
        BER_THRESHOLD = 1e-2

        % ---- Execution / build --------------------------------------
        UseMex  = true
        Rebuild = true
    end

    properties
        Tcr        % CR fixed-point type table
        CordicIts  % CORDIC iterations (= FL)
    end

    %% ================================================================
    %  One-time setup
    %% ================================================================
    methods (TestClassSetup)

        function setupPath(~)
            here = fileparts(mfilename('fullpath'));
            addpath(genpath(fullfile(here, '..', 'src')));
            addpath(fullfile(here, '..', 'build'));
        end

        function seedRng(~)
            rng(20260530);
        end

        function setupTypes(testCase)
            fxp = struct('WL', testCase.IntBits + testCase.FL, 'FL', testCase.FL);
            testCase.Tcr       = carrier_recovery.fxp_types(fxp);
            testCase.CordicIts = testCase.FL;
        end

        function buildMex(testCase)
            % Compile both MEX binaries so their baked-in fixed-point type and
            % CORDIC-iteration constant match this test's configuration exactly.
            if ~testCase.UseMex || ~testCase.Rebuild
                return;
            end
            fxp = struct('WL', testCase.IntBits + testCase.FL, 'FL', testCase.FL);

            P.N_pol          = testCase.N_pol;
            P.BlockLen       = testCase.BlockLen;
            P.PilotLen       = 1;
            P.VV_NTaps       = testCase.NTaps;
            P.M              = testCase.M;
            P.FxpConfig_VV   = fxp;
            P.FxpConfig_PO   = fxp;
            P.StepSize       = testCase.StepSize;
            P.PilotThreshold = testCase.PilotThreshold;
            P.CordicIts      = testCase.CordicIts;

            cfg = coder.config('mex');
            cfg.GenerateReport            = false;
            cfg.SaturateOnIntegerOverflow = false;

            fprintf('  Building viterbiViterbi_fxp_mex (WL=%d, FL=%d, CordicIts=%d)...\n', ...
                testCase.IntBits + testCase.FL, testCase.FL, testCase.CordicIts);
            build_carrier_recovery_viterbiViterbi_fxp_mex(P, cfg);

            fprintf('  Building pilots_only_fxp_mex...\n');
            build_carrier_recovery_pilots_only_fxp_mex(P, cfg);
        end

    end

    %% ================================================================
    %  Tests
    %% ================================================================
    methods (Test)

        function testViterbiViterbi_Fxp(testCase)
            [rxN, pilots, txRefBits, vvFilter] = buildScenario(testCase, true);

            crSym = runVV(testCase, rxN, pilots, vvFilter);
            crSym = resolveAmbiguity(testCase, crSym, txRefBits);
            BER   = computeBER(testCase, crSym, txRefBits);
            fprintf('Viterbi-Viterbi fxp BER = %.3e\n', BER);

            plotConstellation(testCase, rxN, crSym, ...
                sprintf('Viterbi-Viterbi fxp (WL=%d, FL=%d)', ...
                    testCase.IntBits + testCase.FL, testCase.FL), BER);

            testCase.verifyTrue(all(isfinite(crSym(:))), ...
                'VV fxp output contains NaN/Inf.');
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD, ...
                sprintf('VV fxp BER %.2e exceeds threshold %.2e.', ...
                    BER, testCase.BER_THRESHOLD));
        end

        function testPilotsOnly_Fxp(testCase)
            [rxN, pilots, txRefBits, ~] = buildScenario(testCase, false);

            crSym = runPilots(testCase, rxN, pilots);
            crSym = resolveAmbiguity(testCase, crSym, txRefBits);
            BER   = computeBER(testCase, crSym, txRefBits);
            fprintf('Pilots-only fxp BER = %.3e\n', BER);

            plotConstellation(testCase, rxN, crSym, ...
                sprintf('Pilots-only fxp (WL=%d, FL=%d)', ...
                    testCase.IntBits + testCase.FL, testCase.FL), BER);

            testCase.verifyTrue(all(isfinite(crSym(:))), ...
                'Pilots-only fxp output contains NaN/Inf.');
            testCase.verifyLessThan(BER, testCase.BER_THRESHOLD, ...
                sprintf('Pilots-only fxp BER %.2e exceeds threshold %.2e.', ...
                    BER, testCase.BER_THRESHOLD));
        end

    end

    %% ================================================================
    %  Private helpers
    %% ================================================================
    methods (Access = private)

        % ---- Channel + normalise ------------------------------------
        function [rxN, pilots, txRefBits, vvFilter] = buildScenario(testCase, makeFilter)
            % modulate -> AWGN -> 1 MHz phase noise -> normalise to [-1,1].
            BITS_PER_SF = 3586 * 2 * 2;
            txBits = modem.randomBits(testCase.NumSubframes * BITS_PER_SF);
            [symbols, pilotSyms, ~, ~] = modem.modulate(txBits);

            rx = channel.add_awgn(symbols, testCase.SNR_dB);
            rx = channel.add_phase_noise(rx, testCase.Rs, testCase.LW_Hz);

            [rxN, scale] = modem.normalise(rx, testCase.NormPct);

            txRefBits = modem.symbolsToBits(symbols);
            pilots    = buildPilotMatrix(testCase, symbols, pilotSyms);

            if makeFilter
                % Wiener VV filter, matched to the normalised symbol energy.
                symbolsNorm = symbols ./ scale;
                symEnergy   = mean(abs(symbolsNorm(:)).^2);
                vvFilter    = carrier_recovery.genVVFilter(testCase.LW_Hz, testCase.Rs, ...
                    testCase.SNR_dB, symEnergy, testCase.N_pol, testCase.NTaps);
            else
                vvFilter = [];
            end
        end

        function pilots = buildPilotMatrix(testCase, symbols, pilotSyms)
            % One known pilot per CR block, taken from the CPON pilot at the
            % block's first symbol position.
            CPON_BLOCK_LEN = 32;
            CPON_SF_SYMS   = 3712;
            N_CPON_BLOCKS  = 116;
            Nsym    = size(symbols, 1);
            NBlocks = ceil(Nsym / testCase.BlockLen);
            pilots  = zeros(NBlocks, testCase.N_pol);
            for b = 1:NBlocks
                pos       = (b - 1) * testCase.BlockLen + 1;
                posInSf   = mod(pos - 1, CPON_SF_SYMS) + 1;
                cponBlock = min(floor((posInSf - 1) / CPON_BLOCK_LEN) + 1, N_CPON_BLOCKS);
                pilots(b, :) = pilotSyms(cponBlock, :);
            end
        end

        % ---- Fixed-point runners (MEX or interpreted) ---------------
        function crSym = runVV(testCase, rxN, pilots, vvFilter)
            T = testCase.Tcr;
            rx_fi  = cast(rxN,      'like', T.x);
            w_fi   = cast(vvFilter, 'like', T.w);
            pil_fi = cast(pilots,   'like', T.x);

            if testCase.UseMex
                crSym_fi = carrier_recovery.viterbiViterbi_fxp_mex( ...
                    rx_fi, testCase.N_pol, testCase.NTaps, w_fi, pil_fi, ...
                    testCase.BlockLen, double(testCase.StepSize), testCase.PilotThreshold, ...
                    double(testCase.CordicIts), T);
            else
                crSym_fi = carrier_recovery.viterbiViterbi_fxp( ...
                    rx_fi, testCase.N_pol, testCase.NTaps, w_fi, pil_fi, ...
                    testCase.BlockLen, double(testCase.StepSize), testCase.PilotThreshold, ...
                    double(testCase.CordicIts), T);
            end
            crSym = double(crSym_fi);
        end

        function crSym = runPilots(testCase, rxN, pilots)
            T = testCase.Tcr;
            rx_fi  = cast(rxN,    'like', T.x);
            pil_fi = cast(pilots, 'like', T.x);

            if testCase.UseMex
                crSym_fi = carrier_recovery.pilots_only_fxp_mex( ...
                    rx_fi, testCase.N_pol, testCase.BlockLen, pil_fi, ...
                    double(testCase.CordicIts), T);
            else
                crSym_fi = carrier_recovery.pilots_only_fxp( ...
                    rx_fi, testCase.N_pol, testCase.BlockLen, pil_fi, ...
                    double(testCase.CordicIts), T);
            end
            crSym = double(crSym_fi);
        end

        % ---- BER + QPSK phase-ambiguity resolution ------------------
        function bestSym = resolveAmbiguity(~, crSym, txRefBits)
            bestBER = Inf;
            bestSym = crSym;
            for k = 0:3
                rotated = crSym .* exp(-1j * k * pi/2);
                bits    = modem.symbolsToBits(modem.decideSymbols(rotated));
                n       = min(length(txRefBits), length(bits));
                thisBER = sum(txRefBits(1:n) ~= bits(1:n)) / n;
                if thisBER < bestBER
                    bestBER = thisBER;
                    bestSym = rotated;
                end
            end
        end

        function BER = computeBER(~, crSym, txRefBits)
            bits = modem.symbolsToBits(modem.decideSymbols(crSym));
            n    = min(length(txRefBits), length(bits));
            BER  = sum(txRefBits(1:n) ~= bits(1:n)) / n;
        end

        % ---- Recovered-constellation plot ---------------------------
        function plotConstellation(testCase, rxN, crSym, titleStr, BER)
            figure('Name', titleStr, 'Position', [100 100 1100 520]);
            for p = 1:testCase.N_pol
                subplot(2, testCase.N_pol, p);
                plot(real(rxN(:, p)), imag(rxN(:, p)), '.', 'MarkerSize', 2);
                grid on; axis equal; axis([-1.2 1.2 -1.2 1.2]);
                title(sprintf('Normalised RX  \x2013  Pol %d', p));
                xlabel('In-Phase'); ylabel('Quadrature');

                subplot(2, testCase.N_pol, testCase.N_pol + p);
                plot(real(crSym(:, p)), imag(crSym(:, p)), '.', 'MarkerSize', 2);
                grid on; axis equal;
                title(sprintf('Recovered  \x2013  Pol %d', p));
                xlabel('In-Phase'); ylabel('Quadrature');
            end
            sgtitle(sprintf('%s  |  BER = %.2e', titleStr, BER));
        end

    end
end
