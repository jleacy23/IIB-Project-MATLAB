classdef bit_width_sign_only < matlab.unittest.TestCase
%BIT_WIDTH_SIGN_ONLY  FEC-SNR vs adaptive-EQ bit width, SignOnly on/off.
%
%   For every fractional bit width fl in FL_vec the test rebuilds the
%   equalize_fxp MEX (WL = IntBits + fl, FL = fl) and runs an SNR sweep
%   on a uniform-amplitude DP-QPSK signal (no pilots, no training)
%   through AWGN + 100 kHz phase noise + PMD (80 km, 0.5 ps/sqrt(km),
%   5 sections).  The Viterbi-Viterbi carrier recovery downstream runs
%   at high precision ('fixed32') so the only fixed-point variable in
%   the experiment is the adaptive equalizer.
%
%   Pilots are not used: the CPON-spec pilot amplitude (±3 ±3j) breaks
%   CMA convergence because it sits 8× above the data-symbol modulus.
%   We therefore generate plain ±1±1j data and disable VV's pilot-aided
%   cycle-slip correction by setting PilotThreshold = 10·pi (so the
%   round(PhaseDiff/PilotThreshold) term is always zero).  The Pilots
%   matrix is passed as a dummy ones() array purely to satisfy the API.
%
%   For each (fl, SignOnly) pair the post-CR BER is computed across an
%   SNR sweep; the SNR at which the BER crosses FEC_BER = 2e-2 is
%   recovered by linear interpolation (see fecCrossing).
%
%   BER pipeline notes
%     - The adaptive EQ discards the first AEQ_NOut samples (transient).
%     - The BER is computed only over the aligned window
%         tx = symbols(NOut+1 : NOut+Nused, :)
%         rx = vv(1:Nused, :)
%     - Both polarisation swap and the pi/2 phase ambiguity are
%       resolved per-reference-pol before BER is recorded.
%
%   Output:  results/adaptive_equalizer/bit_width_sign_only.mat
%
%   Run with:
%       runtests('bit_width_sign_only')

    %% ================================================================
    %  Parameters
    %% ================================================================
    properties (Constant)

        % --- System ----------------------------------------------------
        Rs            = 30.5       % [GBd]
        N_pol         = 2
        SpS           = 1          % symbol-rate, no pulse shaping

        % --- Channel ---------------------------------------------------
        LW_Hz         = 100e3
        SNR_dB_vec    = 6 : 1 : 22

        % --- PMD -------------------------------------------------------
        L             = 80         % fibre length [km]
        DGDSpec       = 0.5        % PMD coefficient [ps/sqrt(km)]
        N_pmd         = 5          % number of PMD sections

        % --- Bit-width sweep -------------------------------------------
        IntBits       = 16
        FL_vec        = [2, 4, 6, 8, 10, 12, 14, 16]

        % --- Adaptive equalizer ----------------------------------------
        AEQ_NTaps         = 15
        AEQ_Mu            = 1e-3
        AEQ_SingleSpike   = true
        AEQ_N1            = 2000
        AEQ_NOut          = 2016        % multiple of BlockLen for clean block alignment

        % --- Viterbi-Viterbi (held at high precision) ------------------
        VV_FxpConfig    = 'fixed32'
        VV_NTaps        = 10
        BlockLen        = 32
        StepSize        = 32
        PilotThreshold  = 1000 * pi       % disables pilot-aided cycle-slip correction
        CordicIts       = 16

        % --- Test signal length & Monte-Carlo --------------------------
        NTrials         = 6
        Nsym            = 2^14          % symbols per trial per pol

        % --- FEC threshold ---------------------------------------------
        FEC_BER         = 2e-2

        % --- Diagnostics -----------------------------------------------
        %  When true, plot the AEQ + VV output constellations for the
        %  first trial at the highest SNR for each (fl, SignOnly).
        %  Also prints |aeq| and per-pol cross-correlation (a value
        %  close to 1 implies CMA mode-collapse onto a single pol).
        Debug           = true
    end

    properties
        T_vv          % VV fixed-point types (high precision)
        VVFilters     % {1 x NSNR} per-SNR Wiener filters
    end

    %% ================================================================
    %  Test class setup
    %% ================================================================
    methods (TestClassSetup)

        function setupPath(~)
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src'));
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'build'));
        end

        function seedRng(~)
            rng(42);
        end

        function buildVVMex(testCase)
            % Build the Viterbi-Viterbi MEX once at high precision.
            % AEQ MEX is rebuilt per-FL inside the main test.
            P = testCase;

            cfg = coder.config('mex');
            cfg.GenerateReport            = false;
            cfg.SaturateOnIntegerOverflow = false;

            Bp = bit_width_sign_only.buildParams(P, P.VV_FxpConfig, ...
                                                  P.VV_FxpConfig);
            fprintf('Building VV MEX  (%s) ... ', P.VV_FxpConfig);
            build_carrier_recovery_viterbiViterbi_fxp_mex(Bp, cfg);
            fprintf('OK\n');

            testCase.T_vv = carrier_recovery.fxp_types(P.VV_FxpConfig);
        end

        function calibrateVVFilters(testCase)
            % ±1±1j data has |s|² = 2.  VV is a 4th-power estimator so
            % the absolute scale only matters through the noise term.
            symEnergy = 2;
            NSNR = length(testCase.SNR_dB_vec);
            testCase.VVFilters = cell(1, NSNR);
            for si = 1:NSNR
                testCase.VVFilters{si} = carrier_recovery.genVVFilter( ...
                    testCase.LW_Hz, testCase.Rs, testCase.SNR_dB_vec(si), ...
                    symEnergy, testCase.N_pol, testCase.VV_NTaps);
            end
        end
    end

    %% ================================================================
    %  Tests
    %% ================================================================
    methods (Test)

        function test_bit_width_sign_only(testCase)
            P    = testCase;
            NFL  = length(P.FL_vec);
            NSO  = 2;          % SignOnly = {false, true}
            NSNR = length(P.SNR_dB_vec);
            signOnlyVec = [false, true];

            fecSNR  = nan(NFL, NSO);
            berCube = nan(NFL, NSO, NSNR);

            cfg = coder.config('mex');
            cfg.GenerateReport            = false;
            cfg.SaturateOnIntegerOverflow = false;

            for fi_idx = 1:NFL
                fl  = P.FL_vec(fi_idx);
                fxp = struct('WL', P.IntBits + fl, 'FL', fl);

                fprintf('============ FL = %2d  (%d/%d) ============\n', ...
                        fl, fi_idx, NFL);

                % --- Rebuild AEQ MEX for this FL ------------------------
                Bp = bit_width_sign_only.buildParams(P, fxp, P.VV_FxpConfig);
                fprintf('  Building AEQ MEX ... ');
                build_adaptive_eq_equalize_fxp_mex(Bp, cfg);
                fprintf('OK\n');
                bit_width_sign_only.clearWorkerMex();

                T_aeq = adaptive_eq.equalize_fxp_types(fxp);

                % --- SNR sweep × {SignOnly off, on} ---------------------
                for so_idx = 1:NSO
                    signOnly = signOnlyVec(so_idx);
                    fprintf('  SignOnly = %d :', signOnly);

                    ber_acc = zeros(NSNR, 1);
                    for tr = 1:P.NTrials
                        for si = 1:NSNR
                            % Plot diagnostics on the first trial at
                            % the highest SNR of each (fl, SignOnly).
                            if P.Debug && tr == 1 && si == NSNR
                                tag = sprintf( ...
                                    'fl=%d, SignOnly=%d, SNR=%.0f dB', ...
                                    fl, signOnly, P.SNR_dB_vec(si));
                            else
                                tag = '';
                            end
                            ber_acc(si) = ber_acc(si) + ...
                                bit_width_sign_only.runOnce(P, ...
                                    P.SNR_dB_vec(si), signOnly, ...
                                    T_aeq, P.T_vv, P.VVFilters{si}, tag);
                        end
                    end
                    ber_avg = ber_acc / P.NTrials;

                    nBitsTrial = (P.Nsym - P.AEQ_NOut) * P.N_pol * 2;
                    berFloor = 1 / (P.NTrials * nBitsTrial);
                    ber_avg(ber_avg == 0) = berFloor;

                    berCube(fi_idx, so_idx, :) = ber_avg;
                    fecSNR(fi_idx, so_idx) = ...
                        bit_width_sign_only.fecCrossing(P.SNR_dB_vec, ...
                                                        ber_avg, P.FEC_BER);

                    fprintf('  FEC SNR @ %.0e = %.2f dB\n', ...
                            P.FEC_BER, fecSNR(fi_idx, so_idx));
                end
            end

            % --- Save & plot ----------------------------------------
            outDir  = fileparts(mfilename('fullpath'));
            outFile = fullfile(outDir, 'bit_width_sign_only.mat');
            FL_vec_out     = P.FL_vec;     %#ok<NASGU>
            SNR_dB_vec_out = P.SNR_dB_vec; %#ok<NASGU>
            FEC_BER_out    = P.FEC_BER;    %#ok<NASGU>
            save(outFile, 'fecSNR', 'berCube', 'FL_vec_out', ...
                          'SNR_dB_vec_out', 'FEC_BER_out');
            fprintf('Saved results to %s\n', outFile);

            figure('Name', 'CMA bit-width: SignOnly on/off');
            plot(P.FL_vec, fecSNR(:,1), 'o-', 'LineWidth', 1.2, ...
                 'DisplayName', 'SignOnly = false'); hold on;
            plot(P.FL_vec, fecSNR(:,2), 's-', 'LineWidth', 1.2, ...
                 'DisplayName', 'SignOnly = true');
            grid on;
            xlabel('Fractional bit width');
            ylabel(sprintf('SNR @ BER = %.0e [dB]', P.FEC_BER));
            title(sprintf( ...
                'CMA adaptive EQ  |  LW = %.0f kHz  |  R_s = %.1f GBd', ...
                P.LW_Hz/1e3, P.Rs));
            legend('Location', 'best');

            % Verify at least one FEC crossing was found per SignOnly
            testCase.verifyTrue(any(~isnan(fecSNR(:,1))), ...
                'No FEC crossings found for SignOnly = false.');
            testCase.verifyTrue(any(~isnan(fecSNR(:,2))), ...
                'No FEC crossings found for SignOnly = true.');
        end
    end

    %% ================================================================
    %  Static helpers
    %% ================================================================
    methods (Static, Access = private)

        function symbols = randomQPSK(Nsym, N_pol)
            % Uniform ±1±1j data (|s|² = 2).  No pilots, no training.
            re = 2 * randi([0 1], Nsym, N_pol) - 1;
            im = 2 * randi([0 1], Nsym, N_pol) - 1;
            symbols = complex(re, im);
        end

        function Bp = buildParams(P, fxp_aeq, fxp_vv)
            % AEQ side
            Bp.SpS              = P.SpS;
            Bp.AEQ_NTaps        = P.AEQ_NTaps;
            Bp.AEQ_Mu           = P.AEQ_Mu;
            Bp.AEQ_SingleSpike  = P.AEQ_SingleSpike;
            Bp.AEQ_N1           = P.AEQ_N1;
            Bp.AEQ_NOut         = P.AEQ_NOut;
            Bp.AEQ_SignOnly     = false;     % build-time value only
            Bp.FxpConfig_AEQ    = fxp_aeq;

            % VV side
            Bp.N_pol            = P.N_pol;
            Bp.VV_NTaps         = P.VV_NTaps;
            Bp.BlockLen         = P.BlockLen;
            Bp.StepSize         = P.StepSize;
            Bp.PilotThreshold   = P.PilotThreshold;
            Bp.CordicIts        = P.CordicIts;
            Bp.PilotLen         = 1;
            Bp.FxpConfig_VV     = fxp_vv;
        end

        function clearWorkerMex()
            clear mex %#ok<CLMEX>
        end

        function ber = runOnce(P, SNR_dB, SignOnly, T_aeq, T_vv, VVFilter, debugTag)
            if nargin < 7, debugTag = ''; end
            % --- TX: uniform ±1±1j QPSK -----------------------------
            symbols = bit_width_sign_only.randomQPSK(P.Nsym, P.N_pol);

            % --- Channel: AWGN + phase noise + PMD ------------------
            rx = channel.add_awgn(symbols, SNR_dB);
            rx = channel.add_phase_noise(rx, P.Rs, P.LW_Hz);
            rx = channel.add_pmd(rx, P.L, P.SpS, P.Rs, ...
                                 P.DGDSpec, P.N_pmd);

            % --- Effective step size --------------------------------
            %  Nominal Mu (e.g. 1e-3) cannot be represented as a fi
            %  whose LSB exceeds it.  At FL <= 9 with LSB >= 2e-3,
            %  cast(1e-3, 'like', T.mu) floors to zero and no tap
            %  updates occur.  Floor mu to one LSB of T.mu so each
            %  configuration runs at least one quantum of step size;
            %  at FL >= 10 this reduces to the nominal value.
            lsb_mu  = 2^(-double(T_aeq.mu.FractionLength));
            mu_eff  = max(double(P.AEQ_Mu), lsb_mu);

            % --- Adaptive equalisation (fixed-point MEX) ------------
            rx_aeq_fi = cast(rx, 'like', T_aeq.x);
            aeq = adaptive_eq.equalize_fxp_mex(rx_aeq_fi, ...
                    double(P.SpS), double(P.AEQ_NTaps), ...
                    mu_eff, logical(P.AEQ_SingleSpike), ...
                    double(P.AEQ_N1), double(P.AEQ_NOut), ...
                    logical(SignOnly), T_aeq);

            aeq_d = double(aeq);
            if ~isempty(debugTag)
                fprintf('    [%s]  mu_eff = %.3e (nominal %.3e, LSB %.3e)\n', ...
                        debugTag, mu_eff, double(P.AEQ_Mu), lsb_mu);
                bit_width_sign_only.diagPrint(aeq_d, ['AEQ  ' debugTag]);
                bit_width_sign_only.diagPlot(aeq_d,  ['AEQ  ' debugTag]);
            end

            % --- Dummy pilots (correction inert via large threshold) -
            Nsym_eq    = size(aeq, 1);
            NBlocks    = ceil(Nsym_eq / P.BlockLen);
            pilots     = ones(NBlocks, P.N_pol);

            % --- Viterbi-Viterbi (fixed-point MEX) ------------------
            aeq_vv_fi = cast(double(aeq), 'like', T_vv.x);
            vvf_fi    = cast(VVFilter,     'like', T_vv.w);
            pilots_fi = cast(pilots,       'like', T_vv.x);
            [vv, ~] = carrier_recovery.viterbiViterbi_fxp_mex( ...
                        aeq_vv_fi, double(P.N_pol), ...
                        double(P.VV_NTaps), vvf_fi, pilots_fi, ...
                        double(P.BlockLen), double(P.StepSize), ...
                        double(P.PilotThreshold), ...
                        double(P.CordicIts), T_vv);

            vv_d = double(vv);
            if ~isempty(debugTag)
                bit_width_sign_only.diagPlot(vv_d, ['VV   ' debugTag]);
            end

            % --- BER over aligned symbols ---------------------------
            %  Reference window: symbols(NOut+1 : NOut+Nused, :)
            %  Recovered window: vv(1:Nused, :)
            %  Try both pol assignments and all four pi/2 rotations
            %  per reference polarisation; record the best BER.
            Nused     = size(vv_d, 1);
            symOffset = P.AEQ_NOut;
            refEnd    = min(symOffset + Nused, size(symbols, 1));
            Nuse      = refEnd - symOffset;
            refSym    = symbols(symOffset + 1 : refEnd, :);

            rotations = [1, 1j, -1, -1j];
            totalErr  = 0;
            totalBits = 0;
            for p = 1:P.N_pol
                refBitsPol = modem.symbolsToBits(refSym(:, p));
                bestPolBER = Inf;
                for q = 1:P.N_pol          % pol-swap candidate
                    for ri = 1:4           % pi/2 rotation candidate
                        vvRot   = vv_d(1:Nuse, q) * rotations(ri);
                        decRot  = modem.decideSymbols(vvRot);
                        bitsRot = modem.symbolsToBits(decRot);
                        polBER  = sum(bitsRot ~= refBitsPol) / numel(refBitsPol);
                        if polBER < bestPolBER
                            bestPolBER = polBER;
                        end
                    end
                end
                totalErr  = totalErr  + bestPolBER * numel(refBitsPol);
                totalBits = totalBits + numel(refBitsPol);
            end
            ber = totalErr / totalBits;
        end

        function diagPrint(z, tag)
            % Print magnitude stats and per-pol cross-correlation.
            % A |corr| near 1 implies CMA mode-collapse onto a single
            % pol (both columns of z carry the same input pol).
            amp     = abs(z(:));
            zc1     = z(:,1) - mean(z(:,1));
            zc2     = z(:,2) - mean(z(:,2));
            crossC  = abs(sum(zc1 .* conj(zc2))) ...
                      / sqrt(sum(abs(zc1).^2) * sum(abs(zc2).^2));
            fprintf(['    [%s]  |z|: min=%.2e mean=%.2e max=%.2e ', ...
                     ' xcorr(pol1,pol2)=%.3f\n'], ...
                    tag, min(amp), mean(amp), max(amp), crossC);
        end

        function diagPlot(z, tag)
            % Scatter plot of both pols.  Subsample to keep the figure
            % responsive; use the tail of the signal so transients are
            % out of frame.
            n      = size(z, 1);
            tail   = max(1, floor(n/4));
            idx    = (n - tail + 1) : n;
            stride = max(1, floor(tail / 4000));
            zs     = z(idx(1:stride:end), :);

            figure('Name', tag, 'Position', [60, 60, 900, 420]);
            for p = 1:size(zs, 2)
                subplot(1, 2, p);
                plot(real(zs(:,p)), imag(zs(:,p)), '.', 'MarkerSize', 3);
                grid on; axis equal;
                xlabel('I'); ylabel('Q');
                title(sprintf('Pol %d', p));
            end
            sgtitle(tag);
        end

        function xCross = fecCrossing(x, y, yLimit)
            % Linear interpolation between the two adjacent x-points
            % whose y-values straddle yLimit.  Returns NaN if no
            % crossing exists in the sweep.
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
                    xCross = x(i) + (yLimit - y(i)) ...
                             * (x(i+1) - x(i)) / (y(i+1) - y(i));
                    return;
                end
            end
        end

    end
end
