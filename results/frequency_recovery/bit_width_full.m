classdef bit_width_full < matlab.unittest.TestCase
%BIT_WIDTH_FULL  Fixed-point precision sweep for differential-Kay + pilots-only.
%
%   The carrier-recovery chain is fixed to the two winning algorithms:
%       frequency recovery : differential phase and Kay (data-aided)
%       phase recovery     : pilots-only
%
%   Two precision sweeps are run, each holding one subsystem at high precision
%   (FL = FL_fixed) and sweeping the fractional bit width of the other:
%
%     Sweep 1 (phase precision) : pilots-only fractional bits swept,
%                                 differential Kay held at high precision.
%     Sweep 2 (frequency precision) : differential Kay fractional bits swept,
%                                 pilots-only held at high precision.
%
%   The differential-Kay CFO estimate is averaged over NumAvg subframes before
%   the decode subframe is corrected (each estimate from its 11 training
%   symbols), so the residual CFO entering phase recovery reflects the
%   subframe-averaged operating point rather than a single noisy estimate.
%
%   For each precision point the post-recovery BER is found by Monte-Carlo over
%   NTrials at each SNR in SNR_dB_vec; the RSNR --- the SNR at which the BER
%   crosses FEC_BER --- is interpolated per trial and reported as a mean with
%   standard deviation.  Two RSNR-vs-precision graphs (one per sweep) are
%   produced and the data saved to bit_width_full_precision_sweep.mat.
%
%   Word length is WL = IntBits + FL.  The CR (pilots-only) CORDIC iteration
%   count tracks the swept fractional bits; FR (differential Kay) keeps its
%   CORDIC iterations at FR_CordicIts because its phase/accumulator types are
%   pinned wide by freq_recovery.fxp_types.
%
%   Run with:
%       runtests('bit_width_full')

    %% ================================================================
    %  Parameters
    %% ================================================================
    properties (Constant)

        % System
        Rs          = 30.5
        N_pol       = 2
        TrainingLen = 11

        % Monte-Carlo
        NTrials     = 10

        % SNR sweep
        SNR_dB_vec  = 0 : 1 : 30

        % Bit width sweep — WL = IntBits + FL; FL_fixed is the "high precision"
        % held on the subsystem that is not being swept.
        IntBits     = 16
        FL_vec      = [2, 4, 6, 8, 10, 12, 14, 16]
        FL_fixed    = 16

        % Channel conditions
        DeltaF_Hz   = 2e9
        LW_Hz       = 1000e3
        MaxFreq     = 0.1

        % CORDIC iterations for the FR (differential Kay) build.  FR phase and
        % accumulator types are pinned wide, so this stays fixed while the FR
        % signal precision is swept.
        FR_CordicIts = 16

        % FR CFO-estimate averaging — the differential-Kay estimate is averaged
        % over this many subframes (each from its 11 training symbols) before
        % the decode subframe is corrected, matching the operating point of the
        % subframe-averaging study.
        NumAvg = 12

        % Phase recovery (pilots-only)
        BlockLen     = 32

        % FEC threshold
        FEC_BER = 2e-2

    end

    %% ================================================================
    %  Test class setup
    %% ================================================================
    methods (TestClassSetup)

        function setupPath(~)
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src'));
        end

        function setupBuildPath(~)
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'build'));
        end

        function seedRng(~)
            rng(42);
        end

    end

    %% ================================================================
    %  Test
    %% ================================================================
    methods (Test)

        function test_precision_sweep(testCase)
            P    = testCase;
            Prms = bit_width_full.extractParams(testCase);
            NFL  = length(P.FL_vec);

            % ---- Sweep 1: phase precision (pilots-only swept, diff Kay high)
            % 'Uniform' makes the FR phase/accumulator types track T.x precision.
            fxp_fr_hi = struct('WL', P.IntBits + P.FL_fixed, 'FL', P.FL_fixed, 'Uniform', true);
            bit_width_full.buildFRMexInSrc(P, fxp_fr_hi);
            T_fr_hi   = freq_recovery.fxp_types(fxp_fr_hi);

            rsnr_cr_mean = nan(NFL, 1);
            rsnr_cr_std  = nan(NFL, 1);
            for i = 1:NFL
                fl = P.FL_vec(i);
                fprintf('=== Sweep 1 (phase) FL_cr = %2d  (%d / %d) ===\n', fl, i, NFL);
                fxp_cr = struct('WL', P.IntBits + fl, 'FL', fl);
                bit_width_full.buildCRMexInSrc(P, fxp_cr);
                T_cr = carrier_recovery.fxp_types(fxp_cr);

                ber = bit_width_full.runSnrSweep(Prms, T_fr_hi, T_cr);
                cr  = bit_width_full.fecCrossingsAll(P.SNR_dB_vec, ber, P.FEC_BER);
                rsnr_cr_mean(i) = mean(cr, 'omitnan');
                rsnr_cr_std(i)  = std(cr,  'omitnan');
            end

            % ---- Sweep 2: frequency precision (diff Kay swept, pilots-only high)
            fxp_cr_hi = struct('WL', P.IntBits + P.FL_fixed, 'FL', P.FL_fixed);
            bit_width_full.buildCRMexInSrc(P, fxp_cr_hi);
            T_cr_hi   = carrier_recovery.fxp_types(fxp_cr_hi);

            rsnr_fr_mean = nan(NFL, 1);
            rsnr_fr_std  = nan(NFL, 1);
            for i = 1:NFL
                fl = P.FL_vec(i);
                fprintf('=== Sweep 2 (frequency) FL_fr = %2d  (%d / %d) ===\n', fl, i, NFL);
                fxp_fr = struct('WL', P.IntBits + fl, 'FL', fl, 'Uniform', true);
                bit_width_full.buildFRMexInSrc(P, fxp_fr);
                T_fr = freq_recovery.fxp_types(fxp_fr);

                ber = bit_width_full.runSnrSweep(Prms, T_fr, T_cr_hi);
                cr  = bit_width_full.fecCrossingsAll(P.SNR_dB_vec, ber, P.FEC_BER);
                rsnr_fr_mean(i) = mean(cr, 'omitnan');
                rsnr_fr_std(i)  = std(cr,  'omitnan');
            end

            % ---- Save + plot ------------------------------------------
            outDir = fileparts(mfilename('fullpath'));
            results.FL_vec       = P.FL_vec;
            results.rsnr_cr_mean = rsnr_cr_mean;
            results.rsnr_cr_std  = rsnr_cr_std;
            results.rsnr_fr_mean = rsnr_fr_mean;
            results.rsnr_fr_std  = rsnr_fr_std;
            save(fullfile(outDir, 'bit_width_full_precision_sweep.mat'), '-struct', 'results');

            bit_width_full.plotSweep(P.FL_vec, rsnr_cr_mean, rsnr_cr_std, ...
                'Pilots-only precision (differential Kay at high precision)', ...
                fullfile(outDir, 'pilots_precision_sweep.png'));
            bit_width_full.plotSweep(P.FL_vec, rsnr_fr_mean, rsnr_fr_std, ...
                'Differential Kay precision (pilots-only at high precision)', ...
                fullfile(outDir, 'diffkay_precision_sweep.png'));

            testCase.verifyTrue(any(isfinite(rsnr_cr_mean)) && any(isfinite(rsnr_fr_mean)), ...
                'No finite RSNR was recorded in either sweep.');
        end

    end

    %% ================================================================
    %  Static helpers
    %% ================================================================
    methods (Static, Access = private)

        % ---------------- Parameter bundle ----------------------------
        function Params = extractParams(testCase)
            Params.Rs           = testCase.Rs;
            Params.N_pol        = testCase.N_pol;
            Params.NTrials      = testCase.NTrials;
            Params.SNR_dB_vec   = testCase.SNR_dB_vec;
            Params.DeltaF_Hz    = testCase.DeltaF_Hz;
            Params.LW_Hz        = testCase.LW_Hz;
            Params.MaxFreq      = testCase.MaxFreq;
            Params.FR_CordicIts = testCase.FR_CordicIts;
            Params.BlockLen     = testCase.BlockLen;
            Params.TrainingLen  = testCase.TrainingLen;
            Params.NumAvg       = testCase.NumAvg;
        end

        % ---------------- SNR sweep -----------------------------------
        function ber = runSnrSweep(Params, T_fr, T_cr)
            % Returns per-trial BER ber(tr, si) for the differential-Kay FR +
            % pilots-only CR chain at the given fixed-point types.
            NSNR = length(Params.SNR_dB_vec);
            ber  = zeros(Params.NTrials, NSNR);

            % CR CORDIC iterations track the swept CR precision; must equal the
            % coder.Constant baked by buildCRMexInSrc (B.CordicIts = fxp_cr.FL).
            cordicIts_cr = double(T_cr.theta.FractionLength);

            for tr = 1:Params.NTrials
                for si = 1:NSNR
                    [fr_out, pilots, txRefBits] = bit_width_full.buildChannel( ...
                        Params, Params.SNR_dB_vec(si), T_fr);

                    fr_fi     = cast(fr_out, 'like', T_cr.x);
                    pilots_fi = cast(pilots, 'like', T_cr.x);

                    [cr_po, ~] = carrier_recovery.pilots_only_fxp_mex( ...
                        fr_fi, Params.N_pol, Params.BlockLen, pilots_fi, ...
                        cordicIts_cr, T_cr);
                    cr_po = bit_width_full.resolveAmbiguity(double(cr_po), txRefBits);
                    ber(tr, si) = bit_width_full.computeBER(cr_po, txRefBits);
                end
            end
        end

        % ---------------- Channel + FR --------------------------------
        function [fr_out, pilots, txRefBits] = buildChannel(P, SNR_dB, T_fr)
            CPON_BLOCK_LEN = 32;
            CPON_SF_SYMS   = 3712;
            N_CPON_BLOCKS  = 116;
            BITS_PER_SF    = 3586 * 2 * 2;
            L              = P.TrainingLen;

            % --- Decode subframe (carries the data to be recovered) ---
            txBits = modem.randomBits(BITS_PER_SF);
            [symbols, pilotSyms, training, ~] = modem.modulate(txBits);

            rx = channel.lo_freq_shift(symbols, P.DeltaF_Hz / 1e6, P.Rs, 1);
            rx = channel.add_awgn(rx, SNR_dB);
            rx = channel.add_phase_noise(rx, P.Rs, P.LW_Hz);

            tr_fi = cast(training, 'like', T_fr.x);

            % --- Average the differential-Kay CFO estimate over NumAvg
            %     subframes.  Data-aided diffkay needs only the 11 training
            %     symbols, so each estimate is formed from a short training-only
            %     observation (same CFO, independent noise).  The decode
            %     subframe's own training is the first observation.
            f_est = zeros(P.NumAvg, 1);
            for n = 1:P.NumAvg
                if n == 1
                    obs = rx(1:L, :);
                else
                    obs = channel.lo_freq_shift(training, P.DeltaF_Hz / 1e6, P.Rs, 1);
                    obs = channel.add_awgn(obs, SNR_dB);
                    obs = channel.add_phase_noise(obs, P.Rs, P.LW_Hz);
                end
                obs_fi = cast(obs, 'like', T_fr.x);
                [~, f_est(n)] = freq_recovery.differential_kay_fxp_mex( ...
                    obs_fi, tr_fi, P.Rs, P.FR_CordicIts, T_fr, true, 0, P.MaxFreq);
            end
            f_avg = mean(f_est);

            % Apply the averaged CFO correction to the decode subframe.  The
            % swept FR precision enters through the (fixed-point) per-subframe
            % estimates above; the de-rotation by the averaged estimate is done
            % in double.
            fr_out = double(channel.lo_freq_shift(rx, -f_avg / 1e6, P.Rs, 1));

            Nsym    = size(fr_out, 1);
            NBlocks = ceil(Nsym / P.BlockLen);
            pilots  = zeros(NBlocks, P.N_pol);
            for b = 1:NBlocks
                pos     = (b - 1) * P.BlockLen + 1;
                posInSf = mod(pos - 1, CPON_SF_SYMS) + 1;
                blk     = min(floor((posInSf - 1) / CPON_BLOCK_LEN) + 1, N_CPON_BLOCKS);
                pilots(b, :) = pilotSyms(blk, :);
            end

            txRefBits = modem.symbolsToBits(symbols);
        end

        function BER = computeBER(crSym, refBits)
            dec  = modem.decideSymbols(crSym);
            bits = modem.symbolsToBits(dec);
            n    = min(length(refBits), length(bits));
            BER  = sum(refBits(1:n) ~= bits(1:n)) / n;
        end

        function best = resolveAmbiguity(crSym, txRefBits)
            bestBER = Inf;
            best    = crSym;
            for k = 0:3
                rot  = crSym .* exp(-1j * k * pi/2);
                dec  = modem.decideSymbols(rot);
                bits = modem.symbolsToBits(dec);
                n    = min(length(txRefBits), length(bits));
                ber  = sum(txRefBits(1:n) ~= bits(1:n)) / n;
                if ber < bestBER
                    bestBER = ber;
                    best    = rot;
                end
            end
        end

        % ---------------- FEC-SNR crossing ----------------------------
        function snrs = fecCrossingsAll(snrDb, berPerTrial, fecBer)
            % Per-trial FEC SNR crossings. berPerTrial is [NTrials x NSNR].
            NTrials = size(berPerTrial, 1);
            snrs    = nan(NTrials, 1);
            for t = 1:NTrials
                snrs(t) = bit_width_full.fecCrossing(snrDb, berPerTrial(t, :), fecBer);
            end
        end

        function xCross = fecCrossing(x, y, yLimit)
            xCross = NaN;
            n      = length(x);
            if n < 2, return; end

            exact = find(y == yLimit, 1, 'first');
            if ~isempty(exact)
                xCross = x(exact);
                return;
            end

            for i = 1:(n - 1)
                if (y(i) - yLimit) * (y(i+1) - yLimit) < 0
                    xCross = x(i) + (yLimit - y(i)) * (x(i+1) - x(i)) / (y(i+1) - y(i));
                    return;
                end
            end
        end

        % ---------------- MEX builds (into src) -----------------------
        function buildFRMexInSrc(P, fxp_fr)
            clear mex %#ok<CLMEX>
            B.Rs           = P.Rs;
            B.N_pol        = P.N_pol;
            B.TrainingLen  = P.TrainingLen;
            B.FxpConfig_FR = fxp_fr;
            B.CordicIts    = P.FR_CordicIts;
            B.MaxFreq      = P.MaxFreq;

            cfg = coder.config('mex');
            cfg.GenerateReport = false;
            build_freq_recovery_differential_kay_fxp_mex(B, cfg);
        end

        function buildCRMexInSrc(P, fxp_cr)
            clear mex %#ok<CLMEX>
            B.N_pol        = P.N_pol;
            B.BlockLen     = P.BlockLen;
            B.FxpConfig_PO = fxp_cr;
            % pilots-only CORDIC iterations = swept fractional bits.
            B.CordicIts    = fxp_cr.FL;

            cfg = coder.config('mex');
            cfg.GenerateReport = false;
            build_carrier_recovery_pilots_only_fxp_mex(B, cfg);
        end

        % ---------------- Plotting ------------------------------------
        function plotSweep(fl_vec, rsnr_mean, rsnr_std, titleStr, outFile)
            f = figure('Name', titleStr, 'Position', [100 100 700 480]);
            errorbar(fl_vec, rsnr_mean, rsnr_std, '-o', 'LineWidth', 1.5, ...
                'MarkerFaceColor', 'auto');
            grid on;
            xlabel('Fractional bits');
            ylabel('RSNR [dB]');
            title(titleStr);
            xlim([min(fl_vec) - 1, max(fl_vec) + 1]);
            exportgraphics(f, outFile);
        end

    end
end
