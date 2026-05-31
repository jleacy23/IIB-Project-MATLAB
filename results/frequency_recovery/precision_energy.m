classdef precision_energy < matlab.unittest.TestCase
%PRECISION_ENERGY  RSNR and energy of differential-Kay + pilots-only at given
%                   fixed-point precision pairs.
%
%   For each (FR, CR) fractional-precision pair in PrecPairs = [fl_fr, fl_cr]
%   the carrier-recovery chain
%       frequency recovery : differential phase and Kay (data-aided), fl_fr bits
%       phase recovery     : pilots-only,                              fl_cr bits
%   is simulated over an AWGN + CFO + 1 MHz phase-noise channel.  The test
%   reports:
%
%     * RSNR (mean +/- std over trials) --- the SNR at which the post-recovery
%       BER crosses FEC_BER, interpolated per trial.
%     * Energy per bit, as the sum of
%         E_FR     differential-Kay estimate averaged over NumAvg subframes
%                  (op count scaled by NumAvg), amortised over one subframe,
%         E_CR     pilots-only estimate,      amortised over one block,
%         E_apply  the per-symbol CORDIC phase correction, modelled as one real
%                  multiplication at the pilots-only precision fl_cr.
%
%   The differential-Kay CFO estimate is averaged over NumAvg subframes before
%   the decode subframe is corrected, matching the subframe-averaging operating
%   point; the total estimation cost is amortised over a single subframe.
%
%   The differential-Kay FR types are uniform (phase/accumulator at the swept
%   fl_fr precision); the energy model evaluates every operation at the
%   corresponding fractional precision.
%
%   Output: precision_energy.mat (table 'tbl') + a printed summary and an
%           RSNR-vs-energy scatter (precision_energy.png).
%
%   Run:  runtests('precision_energy')

    %% ================================================================
    %  Parameters
    %% ================================================================
    properties (Constant)
        % ---- System -------------------------------------------------
        Rs          = 30.5
        N_pol       = 2
        TrainingLen = 11
        SubframeLen = 3712
        BlockLen    = 32

        % ---- Monte-Carlo --------------------------------------------
        NTrials     = 20
        SNR_dB_vec  = 0 : 1 : 30

        % ---- Channel ------------------------------------------------
        DeltaF_Hz   = 3e9
        LW_Hz       = 1e6
        MaxFreq     = 0.1

        % ---- Precision pairs [fl_fr, fl_cr] (user-specified) --------
        IntBits     = 16
        PrecPairs   = [ 10, 4; ...
                        12, 6]

        % FR CORDIC iterations (FR phase/acc precision is set by fl_fr; the
        % iteration count is held fixed).
        FR_CordicIts = 16

        % Differential-Kay CFO estimate is averaged over this many subframes
        % (each from its 11 training symbols).  The total estimation cost
        % (NumAvg estimates) is still amortised over a single subframe.
        NumAvg = 12

        % ---- FEC threshold ------------------------------------------
        FEC_BER = 2e-2

        % ---- Energy model (horowitz2014computing fit) ---------------
        EAdd_fJ      = 3.16/4
        EMult_fJ     = 3.03/4
        M            = 4
        Oversampling = 1

        % ---- Execution / build --------------------------------------
        UseMex  = true
        Rebuild = false
    end

    %% ================================================================
    %  Setup
    %% ================================================================
    methods (TestClassSetup)

        function setupPath(~)
            here = fileparts(mfilename('fullpath'));
            addpath(fullfile(here, '..', '..', 'src'));
            addpath(fullfile(here, '..', '..', 'build'));
        end

        function seedRng(~)
            rng(42);
        end

    end

    %% ================================================================
    %  Test
    %% ================================================================
    methods (Test)

        function test_precision_energy(testCase)
            P    = testCase;
            Prms = precision_energy.extractParams(testCase);
            nP   = size(P.PrecPairs, 1);

            fl_fr     = P.PrecPairs(:, 1);
            fl_cr     = P.PrecPairs(:, 2);
            rsnr_mean = nan(nP, 1);
            rsnr_std  = nan(nP, 1);
            E_fr      = nan(nP, 1);
            E_cr      = nan(nP, 1);
            E_apply   = nan(nP, 1);
            E_total   = nan(nP, 1);

            for i = 1:nP
                ffr = fl_fr(i);
                fcr = fl_cr(i);
                fprintf('=== Pair %d/%d : FR FL = %d, CR FL = %d ===\n', i, nP, ffr, fcr);

                % --- Build the two MEX at this precision pair ----------
                fxp_fr = struct('WL', P.IntBits + ffr, 'FL', ffr, 'Uniform', true);
                fxp_cr = struct('WL', P.IntBits + fcr, 'FL', fcr);
                precision_energy.buildFRMexInSrc(P, fxp_fr);
                precision_energy.buildCRMexInSrc(P, fxp_cr);
                T_fr = freq_recovery.fxp_types(fxp_fr);
                T_cr = carrier_recovery.fxp_types(fxp_cr);

                % --- RSNR over trials ----------------------------------
                ber = precision_energy.runSnrSweep(Prms, T_fr, T_cr);
                cr  = precision_energy.fecCrossingsAll(P.SNR_dB_vec, ber, P.FEC_BER);
                rsnr_mean(i) = mean(cr, 'omitnan');
                rsnr_std(i)  = std(cr,  'omitnan');

                % --- Energy per bit ------------------------------------
                E_fr(i)    = precision_energy.frEnergy(P, ffr);
                E_cr(i)    = precision_energy.crEnergy(P, fcr);
                E_apply(i) = precision_energy.applyEnergy(P, fcr);
                E_total(i) = E_fr(i) + E_cr(i) + E_apply(i);

                fprintf('  RSNR = %.2f +/- %.2f dB | E = %.2f fJ/bit (FR %.3f, CR %.3f, apply %.3f)\n', ...
                    rsnr_mean(i), rsnr_std(i), E_total(i), E_fr(i), E_cr(i), E_apply(i));
            end

            % --- Table + save ------------------------------------------
            tbl = table(fl_fr, fl_cr, rsnr_mean, rsnr_std, E_fr, E_cr, E_apply, E_total);
            outDir = fileparts(mfilename('fullpath'));
            save(fullfile(outDir, 'precision_energy.mat'), 'tbl');
            fprintf('\n');
            disp(tbl);

            % --- RSNR vs energy scatter --------------------------------
            fig = figure('Name', 'Precision vs energy', 'Position', [100 100 760 520]);
            errorbar(E_total, rsnr_mean, rsnr_std, 'o', 'LineStyle', 'none', ...
                'MarkerFaceColor', 'auto', 'LineWidth', 1.2);
            hold on; grid on;
            for i = 1:nP
                text(E_total(i), rsnr_mean(i), ...
                    sprintf('  %d/%d', fl_fr(i), fl_cr(i)), 'FontSize', 9);
            end
            xlabel('Energy per bit [fJ]');
            ylabel('RSNR [dB]');
            title('RSNR vs energy at (FR/CR) fractional-precision pairs');
            exportgraphics(fig, fullfile(outDir, 'precision_energy.png'));

            testCase.verifyTrue(any(isfinite(rsnr_mean)), ...
                'No finite RSNR recorded for any precision pair.');
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

        % ---------------- Energy --------------------------------------
        function E = frEnergy(P, fl)
            % Differential-Kay estimate.  The estimate is averaged over NumAvg
            % subframes, so its op count is scaled by NumAvg; the total is then
            % amortised over a single subframe.
            L  = P.TrainingLen;
            NM = P.NumAvg * (7 * L + 1) / P.SubframeLen;
            NA = P.NumAvg * (4 * L - 2) / P.SubframeLen;
            E  = energy.receiver(NA, NM, P.EAdd_fJ, P.EMult_fJ, P.M, P.Oversampling, fl);
        end

        function E = crEnergy(P, fl)
            % Pilots-only estimate (1 CORDIC = 1 mult + 1 add), amortised over
            % one block.
            NM = 1 / P.BlockLen;
            NA = 1 / P.BlockLen;
            E  = energy.receiver(NA, NM, P.EAdd_fJ, P.EMult_fJ, P.M, P.Oversampling, fl);
        end

        function E = applyEnergy(P, fl)
            % Per-symbol CORDIC phase correction = one real multiplication at
            % the pilots-only precision (no amortisation: applied every symbol).
            E = energy.receiver(0, 1, P.EAdd_fJ, P.EMult_fJ, P.M, P.Oversampling, fl);
        end

        % ---------------- SNR sweep -----------------------------------
        function ber = runSnrSweep(Params, T_fr, T_cr)
            NSNR = length(Params.SNR_dB_vec);
            ber  = zeros(Params.NTrials, NSNR);
            cordicIts_cr = double(T_cr.theta.FractionLength);

            for tr = 1:Params.NTrials
                for si = 1:NSNR
                    [fr_out, pilots, txRefBits] = precision_energy.buildChannel( ...
                        Params, Params.SNR_dB_vec(si), T_fr);

                    fr_fi     = cast(fr_out, 'like', T_cr.x);
                    pilots_fi = cast(pilots, 'like', T_cr.x);

                    [cr_po, ~] = carrier_recovery.pilots_only_fxp_mex( ...
                        fr_fi, Params.N_pol, Params.BlockLen, pilots_fi, ...
                        cordicIts_cr, T_cr);
                    cr_po = precision_energy.resolveAmbiguity(double(cr_po), txRefBits);
                    ber(tr, si) = precision_energy.computeBER(cr_po, txRefBits);
                end
            end
        end

        % ---------------- Channel + FR (averaged over NumAvg subframes)
        function [fr_out, pilots, txRefBits] = buildChannel(P, SNR_dB, T_fr)
            CPON_BLOCK_LEN = 32;
            CPON_SF_SYMS   = 3712;
            N_CPON_BLOCKS  = 116;
            BITS_PER_SF    = 3586 * 2 * 2;
            L              = P.TrainingLen;

            % Decode subframe (carries the data to be recovered).
            txBits = modem.randomBits(BITS_PER_SF);
            [symbols, pilotSyms, training, ~] = modem.modulate(txBits);

            rx = channel.lo_freq_shift(symbols, P.DeltaF_Hz / 1e6, P.Rs, 1);
            rx = channel.add_awgn(rx, SNR_dB);
            rx = channel.add_phase_noise(rx, P.Rs, P.LW_Hz);

            tr_fi = cast(training, 'like', T_fr.x);

            % Average the differential-Kay CFO estimate over NumAvg subframes
            % (each from its 11 training symbols, same CFO, independent noise);
            % the decode subframe's training is the first observation.
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

            % Remove the averaged CFO from the decode subframe (applied in
            % double; the swept FR precision enters via the fixed-point
            % per-subframe estimates above).
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
            NTrials = size(berPerTrial, 1);
            snrs    = nan(NTrials, 1);
            for t = 1:NTrials
                snrs(t) = precision_energy.fecCrossing(snrDb, berPerTrial(t, :), fecBer);
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
            B.CordicIts    = fxp_cr.FL;   % pilots-only CORDIC its = its precision
            cfg = coder.config('mex');
            cfg.GenerateReport = false;
            build_carrier_recovery_pilots_only_fxp_mex(B, cfg);
        end

    end
end
