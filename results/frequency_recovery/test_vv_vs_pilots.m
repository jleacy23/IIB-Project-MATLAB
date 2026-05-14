classdef test_vv_vs_pilots < matlab.unittest.TestCase
%TEST_VV_VS_PILOTS  Compare Viterbi-Viterbi vs Pilots-only carrier
%   recovery across linewidths and SNRs.
%
%   Channel: AWGN + Wiener phase noise (no CFO, no fibre impairments),
%   so the comparison isolates the phase-recovery stage.  Uses the
%   floating-point V&V and Pilots-only implementations under
%   src/+carrier_recovery to keep the test self-contained (no MEX
%   build required).
%
%   Produces a BER vs SNR plot with V&V (solid, circles) and Pilots-only
%   (dashed, squares) for each linewidth, with a horizontal line at the
%   CPON pre-FEC limit (2e-2).
%
%   Run with:
%       runtests('test_vv_vs_pilots')

    %% ================================================================
    %  Parameters
    %% ================================================================
    properties (Constant)
        Rs           = 30.5             % symbol rate [GBd]
        N_pol        = 2

        % Monte-Carlo trials per (LW, SNR) point
        NTrials      = 50

        % Sweep grid
        LW_Hz_vec    = [100e3, 1e6, 10e6]
        SNR_dB_vec   = [0:1:30]

        % Carrier recovery
        BlockLen       = 32
        StepSize       = 32
        PilotThreshold = 5 * pi / 9
        VV_NTaps       = 10

        % CPON pre-FEC limit
        FEC_BER = 2e-2
    end

    properties
        SymEnergy
    end

    %% ================================================================
    %  Test class setup
    %% ================================================================
    methods (TestClassSetup)

        function setupPath(~)
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src'));
        end

        function seedRng(~)
            rng(42);
        end

        function calibrateSymbolEnergy(testCase)
            BITS_PER_SF = 3586 * 2 * 2;
            [tmp, ~, ~, ~] = modem.modulate(modem.randomBits(BITS_PER_SF));
            testCase.SymEnergy = mean(abs(tmp(:)).^2);
        end

    end

    %% ================================================================
    %  Tests
    %% ================================================================
    methods (Test)

        function test_compare_linewidths(testCase)
            NLW  = numel(testCase.LW_Hz_vec);
            NSNR = numel(testCase.SNR_dB_vec);

            ber_vv = nan(NLW, NSNR);
            ber_po = nan(NLW, NSNR);

            for li = 1:NLW
                lw = testCase.LW_Hz_vec(li);
                for si = 1:NSNR
                    snr = testCase.SNR_dB_vec(si);

                    VVFilter = carrier_recovery.genVVFilter(lw, testCase.Rs, snr, ...
                        testCase.SymEnergy, testCase.N_pol, testCase.VV_NTaps);

                    vv_trials = zeros(testCase.NTrials, 1);
                    po_trials = zeros(testCase.NTrials, 1);
                    for tr = 1:testCase.NTrials
                        [rxSym, pilots, refBits] = test_vv_vs_pilots.buildChannel( ...
                            testCase, snr, lw);

                        cr_vv = carrier_recovery.viterbiViterbi(rxSym, testCase.N_pol, ...
                            VVFilter, testCase.BlockLen, testCase.StepSize, ...
                            pilots, testCase.PilotThreshold);
                        cr_vv = test_vv_vs_pilots.resolveAmbiguity(cr_vv, refBits);
                        vv_trials(tr) = test_vv_vs_pilots.computeBER(cr_vv, refBits);

                        cr_po = carrier_recovery.pilots_only(rxSym, testCase.N_pol, ...
                            testCase.BlockLen, pilots);
                        cr_po = test_vv_vs_pilots.resolveAmbiguity(cr_po, refBits);
                        po_trials(tr) = test_vv_vs_pilots.computeBER(cr_po, refBits);
                    end
                    ber_vv(li, si) = mean(vv_trials);
                    ber_po(li, si) = mean(po_trials);

                    fprintf('LW = %7.0f Hz, SNR = %2d dB : V&V = %.2e, Pilots = %.2e\n', ...
                        lw, snr, ber_vv(li, si), ber_po(li, si));
                end
            end

            % Floor zero-BER samples so they remain visible on the log axis
            BITS_PER_SF = 3586 * 2 * 2;
            berFloor = 0.5 / (testCase.NTrials * BITS_PER_SF);
            ber_vv(ber_vv == 0) = berFloor;
            ber_po(ber_po == 0) = berFloor;

            test_vv_vs_pilots.plotResults(testCase, ber_vv, ber_po);
        end

    end

    %% ================================================================
    %  Static helpers
    %% ================================================================
    methods (Static, Access = private)

        function [rxSym, pilots, refBits] = buildChannel(P, SNR_dB, LW_Hz)
            CPON_BLOCK_LEN = 32;
            CPON_SF_SYMS   = 3712;
            N_CPON_BLOCKS  = 116;
            BITS_PER_SF    = 3586 * 2 * 2;

            txBits = modem.randomBits(BITS_PER_SF);
            [symbols, pilotSyms, ~, ~] = modem.modulate(txBits);

            rx = channel.add_awgn(symbols, SNR_dB);
            rx = channel.add_phase_noise(rx, P.Rs, LW_Hz);
            rxSym = rx;

            Nsym    = size(rx, 1);
            NBlocks = ceil(Nsym / P.BlockLen);
            pilots  = zeros(NBlocks, P.N_pol);
            for b = 1:NBlocks
                pos       = (b - 1) * P.BlockLen + 1;
                posInSf   = mod(pos - 1, CPON_SF_SYMS) + 1;
                cponBlock = min(floor((posInSf - 1) / CPON_BLOCK_LEN) + 1, N_CPON_BLOCKS);
                pilots(b, :) = pilotSyms(cponBlock, :);
            end

            refBits = modem.symbolsToBits(symbols);
        end

        function BER = computeBER(crSym, refBits)
            dec  = modem.decideSymbols(crSym);
            bits = modem.symbolsToBits(dec);
            n    = min(length(refBits), length(bits));
            BER  = sum(refBits(1:n) ~= bits(1:n)) / n;
        end

        function best = resolveAmbiguity(crSym, refBits)
            bestBER = Inf;
            best    = crSym;
            for k = 0:3
                rot  = crSym .* exp(-1j * k * pi/2);
                dec  = modem.decideSymbols(rot);
                bits = modem.symbolsToBits(dec);
                n    = min(length(refBits), length(bits));
                ber  = sum(refBits(1:n) ~= bits(1:n)) / n;
                if ber < bestBER
                    bestBER = ber;
                    best    = rot;
                end
            end
        end

        function plotResults(P, ber_vv, ber_po)
            NLW      = numel(P.LW_Hz_vec);
            berFloor = 1 / (P.NTrials * 3586 * 4 * 2);   % 1 error / total bits sent (×2 polarisations)

            figure('Name', 'V&V vs Pilots-only — BER across linewidths', ...
                'Color', 'w', 'Position', [100 100 900 600]);
            ax = axes;
            hold(ax, 'on'); grid(ax, 'on');
            colors = lines(NLW);

            for li = 1:NLW
                lwHz = P.LW_Hz_vec(li);
                if lwHz >= 1e6
                    lwStr = sprintf('%g MHz', lwHz / 1e6);
                else
                    lwStr = sprintf('%g kHz', lwHz / 1e3);
                end
                plot(ax, P.SNR_dB_vec, ber_vv(li, :), ...
                    'Color', colors(li, :), 'LineStyle', '-', 'Marker', 'o', ...
                    'LineWidth', 1.8, 'MarkerSize', 4, 'MarkerFaceColor', colors(li, :), ...
                    'DisplayName', sprintf('V&V, LW = %s', lwStr));
                plot(ax, P.SNR_dB_vec, ber_po(li, :), ...
                    'Color', colors(li, :), 'LineStyle', '--', 'Marker', 'o', ...
                    'LineWidth', 1.8, 'MarkerSize', 4, ...
                    'DisplayName', sprintf('Pilots, LW = %s', lwStr));
            end
            yline(ax, P.FEC_BER, 'k:', 'LineWidth', 1.5, ...
                'DisplayName', sprintf('FEC limit (%.0e)', P.FEC_BER));
            yline(ax, berFloor, 'k--', 'LineWidth', 1.0, ...
                'DisplayName', sprintf('BER floor (%.0e)', berFloor));

            set(ax, 'YScale', 'log', 'FontSize', 11, 'Box', 'on');
            xlabel(ax, 'SNR [dB]', 'FontSize', 12);
            ylabel(ax, 'BER', 'FontSize', 12);
            title(ax,  'V&V vs Pilots-only carrier recovery', ...
                'FontSize', 12, 'Interpreter', 'none');
            legend(ax, 'Location', 'best', 'FontSize', 9, 'Interpreter', 'none');

            outDir = fileparts(mfilename('fullpath'));
            saveas(gcf, fullfile(outDir, 'test_vv_vs_pilots_ber.png'));
        end

    end
end
