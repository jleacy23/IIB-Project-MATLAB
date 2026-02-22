function results = run_phase_recovery(P)
%RUN_PHASE_RECOVERY  BER vs block size sweep for block-based fixed-point
%                    carrier phase recovery (VV and BPS).
%
%   results = run_phase_recovery(P)
%
%   Input
%     P - parameter struct from phase_recovery_params().
%         If omitted, phase_recovery_params() is called automatically.
%
%   Output
%     results - struct with fields:
%       .BlockSizes   - swept block sizes [1 x NBlocks]
%       .BER_VV       - mean BER across trials, VV fxp  [1 x NBlocks]
%       .BER_BPS      - mean BER across trials, BPS fxp [1 x NBlocks]
%       .BER_VV_all   - per-trial BER, VV fxp  [NTrials x NBlocks]
%       .BER_BPS_all  - per-trial BER, BPS fxp [NTrials x NBlocks]
%
%   The test holds BlockBased = true throughout so that phase is held
%   constant for the duration of each block.  The block size therefore
%   directly controls the tracking bandwidth: smaller blocks track faster
%   at the cost of more pilot overhead; larger blocks reduce overhead but
%   increase phase error accumulation between refreshes.
%
%   MEX binaries are compiled once before the sweep begins.  The build
%   uses a temporary pipeline_params struct whose relevant fields are
%   overridden from P.  Note that the MEX binary is compiled with the
%   largest block size in P.BlockSizes so that the compiled Pilots vector
%   dimension (P.PilotLen) is consistent across all block sizes — pilots
%   are always taken from the start of each block regardless of block size.

    if nargin < 1 || isempty(P)
        P = phase_recovery_params();
    end

    %% ----------------------------------------------------------------
    %  Validate inputs
    %% ----------------------------------------------------------------
    assert(all(P.BlockSizes > P.PilotLen), ...
        'All BlockSizes must be greater than PilotLen (%d).', P.PilotLen);
    assert(mod(P.BPS_B, 2) == 0, 'BPS_B must be even.');

    %% ----------------------------------------------------------------
    %  Build fixed-point types tables
    %% ----------------------------------------------------------------
    T_vv  = cr_viterbiViterbi_fxp_types(P.FxpConfig);
    T_bps = cr_bps_fxp_types(P.FxpConfig);

    %% ----------------------------------------------------------------
    %  Compile MEX binaries
    %  Use the smallest block size for the pilot vector dimension — it is
    %  fixed at PilotLen regardless of block size, so any block size works.
    %% ----------------------------------------------------------------
    fprintf('=== Compiling MEX binaries ===\n');

    Pbuild.N_pol        = P.N_pol;
    Pbuild.M            = P.M;
    Pbuild.PilotLen     = P.PilotLen;
    Pbuild.VV_NTaps     = P.VV_NTaps;
    Pbuild.BPS_N        = P.BPS_N;
    Pbuild.BPS_B        = P.BPS_B;
    Pbuild.BlockLen     = P.BlockSizes(1);  % representative; affects only
                                            % loop-bound inference, not
                                            % the Pilots vector dimension
    Pbuild.FxpConfig_VV  = P.FxpConfig;
    Pbuild.FxpConfig_BPS = P.FxpConfig;

    cfg = coder.config('mex');
    cfg.GenerateReport           = false;
    cfg.SaturateOnIntegerOverflow = false;

    build_cr_viterbiViterbi_fxp_mex(Pbuild, cfg);
    build_cr_bps_fxp_mex(Pbuild, cfg);
    fprintf('=== MEX compilation complete ===\n\n');

    %% ----------------------------------------------------------------
    %  Pre-compute the VV filter (fixed for the sweep; block size does
    %  not affect the filter since BlockBased only controls output, not
    %  the estimation window).
    %% ----------------------------------------------------------------
    % Use a representative symbol energy of 1 (unit average power QAM)
    symEnergy = 1.0;
    VVFilter  = cr_genVVFilter(P.LW, P.Rs, P.SNR_dB, symEnergy, ...
                               P.N_pol, P.VV_NTaps);
    VVFilter_fi = cast(VVFilter, 'like', T_vv.w);

    %% ----------------------------------------------------------------
    %  Sweep
    %% ----------------------------------------------------------------
    NB = length(P.BlockSizes);
    BER_VV_all  = zeros(P.NTrials, NB);
    BER_BPS_all = zeros(P.NTrials, NB);

    for bi = 1:NB
        blockLen = P.BlockSizes(bi);
        fprintf('Block size %4d  (%d/%d)\n', blockLen, bi, NB);

        for tr = 1:P.NTrials

            %% -- Generate channel -----------------------------------
            k     = log2(P.M);
            Nbits = k * P.N_pol * P.Ns;
            txBits = qam_randomBits(Nbits, blockLen, P.PilotLen, P.M);
            [symbols, pilots] = qam_modulate(txBits, P.M, P.N_pol, P.PilotLen);

            rxSym = channel_add_awgn(symbols, P.SNR_dB);
            rxSym = channel_add_phase_noise(rxSym, P.Rs, P.LW);

            %% -- Cast to fi -----------------------------------------
            rxSym_fi_vv  = cast(rxSym, 'like', T_vv.x);
            rxSym_fi_bps = cast(rxSym, 'like', T_bps.x);
            pilots_fi_vv  = cast(pilots, 'like', T_vv.x);
            pilots_fi_bps = cast(pilots, 'like', T_bps.x);

            %% -- VV fxp MEX -----------------------------------------
            [crSym_vv_fi, ~] = cr_viterbiViterbi_fxp_mex( ...
                rxSym_fi_vv, P.N_pol, P.VV_NTaps, VVFilter_fi, ...
                pilots_fi_vv, blockLen, P.UsePilots, true, T_vv);

            crSym_vv = resolvePhaseAmbiguity(double(crSym_vv_fi), txBits, P.M, P.N_pol);
            BER_VV_all(tr, bi) = computeBER(crSym_vv, txBits, P.M, P.N_pol);

            %% -- BPS fxp MEX ----------------------------------------
            [crSym_bps_fi, ~] = cr_bps_fxp_mex( ...
                rxSym_fi_bps, P.BPS_N, P.N_pol, P.M, P.BPS_B, blockLen, ...
                pilots_fi_bps, P.UsePilots, true, T_bps);

            crSym_bps = resolvePhaseAmbiguity(double(crSym_bps_fi), txBits, P.M, P.N_pol);
            BER_BPS_all(tr, bi) = computeBER(crSym_bps, txBits, P.M, P.N_pol);

            fprintf('  Trial %d/%d  |  VV BER = %.2e  |  BPS BER = %.2e\n', ...
                tr, P.NTrials, BER_VV_all(tr, bi), BER_BPS_all(tr, bi));
        end
    end

    %% ----------------------------------------------------------------
    %  Aggregate results
    %% ----------------------------------------------------------------
    results.BlockSizes  = P.BlockSizes;
    results.BER_VV      = mean(BER_VV_all,  1);
    results.BER_BPS     = mean(BER_BPS_all, 1);
    results.BER_VV_all  = BER_VV_all;
    results.BER_BPS_all = BER_BPS_all;

    %% ----------------------------------------------------------------
    %  Plot
    %% ----------------------------------------------------------------
    if P.Plot
        plotResults(results, P);
    end

    fprintf('\n=== Sweep complete ===\n');
end

% ====================================================================
%  Local helpers
% ====================================================================

function BER = computeBER(crSym, txBits, M, NPol)
    decidedSyms = qam_decideSymbols(crSym, M, NPol);
    rxBits      = qam_symbolsToBits(decidedSyms, M);
    BER         = sum(txBits ~= rxBits) / length(txBits);
end

function bestSym = resolvePhaseAmbiguity(crSym, txBits, M, NPol)
    bestBER = Inf;
    bestSym = crSym;
    for k = 0:3
        rotated     = crSym .* exp(-1j * k * pi/2);
        decidedSyms = qam_decideSymbols(rotated, M, NPol);
        rxBits      = qam_symbolsToBits(decidedSyms, M);
        thisBER     = sum(txBits ~= rxBits) / length(txBits);
        if thisBER < bestBER
            bestBER = thisBER;
            bestSym = rotated;
        end
    end
end

function plotResults(results, P)
    figure('Name', 'BER vs Block Size', 'Position', [100 100 800 550]);

    % Per-trial scatter points
    NB = length(results.BlockSizes);
    for bi = 1:NB
        x_jitter = results.BlockSizes(bi) * ones(P.NTrials, 1);
        semilogy(x_jitter, results.BER_VV_all(:, bi), ...
            'o', 'Color', [0.3 0.6 1.0], 'MarkerSize', 4, ...
            'HandleVisibility', 'off');
        hold on;
        semilogy(x_jitter, results.BER_BPS_all(:, bi), ...
            's', 'Color', [1.0 0.5 0.2], 'MarkerSize', 4, ...
            'HandleVisibility', 'off');
    end

    % Mean BER lines
    semilogy(results.BlockSizes, results.BER_VV,  'o-', ...
        'Color', [0.1 0.4 0.9], 'LineWidth', 2, ...
        'MarkerFaceColor', [0.1 0.4 0.9], 'MarkerSize', 7, ...
        'DisplayName', sprintf('VV fxp (%s)', P.FxpConfig));
    semilogy(results.BlockSizes, results.BER_BPS, 's-', ...
        'Color', [0.9 0.35 0.05], 'LineWidth', 2, ...
        'MarkerFaceColor', [0.9 0.35 0.05], 'MarkerSize', 7, ...
        'DisplayName', sprintf('BPS fxp (%s)', P.FxpConfig));

    % Formatting
    set(gca, 'XScale', 'log', 'XTick', results.BlockSizes, ...
        'XTickLabel', arrayfun(@num2str, results.BlockSizes, ...
                               'UniformOutput', false));
    grid on;
    xlabel('Block size [symbols]');
    ylabel('BER');
    legend('Location', 'best');
    title(sprintf(['BER vs Block Size  |  Block-based CR  |  ' ...
        '%d-QAM  |  SNR = %.1f dB  |  LW = %.0f kHz  |  ' ...
        '%d trials'], ...
        P.M, P.SNR_dB, P.LW/1e3, P.NTrials));

    % Mark PilotLen overhead on top axis as context
    ax1 = gca;
    overhead = P.PilotLen ./ results.BlockSizes * 100;
    ax2 = axes('Position', ax1.Position, ...
                'XAxisLocation', 'top', ...
                'YAxisLocation', 'right', ...
                'Color', 'none', ...
                'XScale', 'log');
    ax2.XTick     = results.BlockSizes;
    ax2.XLim      = ax1.XLim;
    ax2.YTick     = [];
    ax2.XTickLabel = arrayfun(@(v) sprintf('%.0f%%', v), overhead, ...
                              'UniformOutput', false);
    xlabel(ax2, 'Pilot overhead');
    linkaxes([ax1 ax2], 'x');
    axes(ax1);  % return focus to main axes
end