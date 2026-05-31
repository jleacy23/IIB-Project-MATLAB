function process_pipeline_fxp_sweep(varargin)
%PROCESS_PIPELINE_FXP_SWEEP  RSNR + energy analysis of the full receiver sweep.
%
%   process_pipeline_fxp_sweep()
%   process_pipeline_fxp_sweep('MatFile', path, 'FECBER', 2e-2, ...
%                              'Node', '14nm', 'M', 4)
%
%   Loads pipeline_fxp_sweep.mat (the output of the pipeline_fxp_sweep test)
%   and produces three figures + a printed energy table.  The sweep runs
%   FOUR DSP combos — the two precision rows from the equalisation/recovery
%   chapter tables of report/full/full.tex,
%
%        low  precision : StaticFL/ClkFL/AdaptFL/FRFL/CRFL = 8/6/4/10/4
%        high precision :                                  = 10/8/6/12/6
%
%   each with power-of-two FFT twiddles OFF and ON (Designs rows 1..4) — under
%   an OUTER sweep over the ADC effective number of bits (ENOB_vec).
%
%   Outputs
%     1. pipeline_rsnr_vs_enob.png
%          RSNR (the SNR at which the post-receiver BER crosses the FEC
%          threshold, found per trial by log-linear interpolation) vs ENOB,
%          one errorbar line per combo (marker = mean, bar = std across
%          trials).
%     2. pipeline_rel_rsnr_energy.png
%          Scatter of each combo's actual RSNR (dB) against its absolute
%          total receiver energy (pJ/bit), at the highest ENOB.
%     3. pipeline_energy_breakdown.png
%          Stacked bar of the per-bit energy of each combo at the highest
%          ENOB, grouped into four functional blocks: front end
%          (ADC + GSOP + deskew), equalisation (static + adaptive), clock
%          recovery (Godard), and carrier recovery (frequency + phase).
%
%   Energy model
%     Per-symbol real-multiplication (RM) and real-addition (RA) counts come
%     from report/full/full.tex:
%       tab:cd_cost      static CD + matched-filter overlap-save (po2 aware),
%       tab:clk_cost     Modified-Godard metric + CORDIC + FD phase ramp,
%       tab:aeq_cost     sign-sign butterfly CMA,
%       tab:diffkay_cost differential-phase-and-Kay frequency recovery,
%       tab:pilot_cost   pilots-only carrier recovery,
%       tab:gsop_cost    Gram-Schmidt orthogonalisation (front end),
%       tab:deskew_cost  Lagrange fractional-delay deskew (front end).
%     Each stage's energy/bit is energy.receiver(NA, NM, E_A, E_M, M, 1, n)
%     (the per-symbol counts already fold in the oversampling factor eta, so
%     Oversampling = 1 is passed).  The DSP stages run at their own swept
%     word length n = IntBits + FL (the *_wl columns of the table), matching
%     process_combined_eq_clk_fxp_sweep.  The front-end stages (GSOP, deskew)
%     are assumed to run at the highest ENOB (n = max(ENOB_vec)).  The ADCs
%     contribute a flat EADC pJ per bit (default 1), independent of ENOB.
%
%   Name/Value options
%     'MatFile'     - path to the .mat (default: alongside this script)
%     'FECBER'      - FEC threshold used to score designs (default 2e-2)
%     'Node'        - '14nm' (scaled estimate, default) or '45nm' (calibrated)
%     'M'           - modulation order (default 4 for QPSK)
%     'EADC_pJ'     - ADC energy per converted bit [pJ] (default 1)
%     'DeskewOrder' - Lagrange interpolator order N for the deskew FIR
%                     (N_DS = N + 1 taps; default 4 -> 5-tap FIR)
%     'MaxENOB'     - only plot ENOB values up to this (default 8); also
%                     sets the "highest ENOB" used for figures 2 and 3
%     'SavePlots'   - true to write the three PNGs alongside the .mat
%                     (default true)

    here = fileparts(mfilename('fullpath'));
    p = inputParser;
    p.addParameter('MatFile',     fullfile(here, 'pipeline_fxp_sweep.mat'));
    p.addParameter('FECBER',      2e-2);
    p.addParameter('Node',        '14nm');
    p.addParameter('M',           4);
    p.addParameter('EADC_pJ',     1);
    p.addParameter('DeskewOrder', 4);
    p.addParameter('MaxENOB',     8);
    p.addParameter('SavePlots',   true);
    p.parse(varargin{:});
    matFile     = p.Results.MatFile;
    fecBer      = p.Results.FECBER;
    node        = p.Results.Node;
    M           = p.Results.M;
    eadcPJ      = p.Results.EADC_pJ;
    deskewOrder = p.Results.DeskewOrder;
    maxEnob     = p.Results.MaxENOB;
    savePlots   = p.Results.SavePlots;

    if ~isfile(matFile)
        error('process_pipeline_fxp_sweep:missingMat', ...
              ['Could not find %s.\n', ...
               'Run runtests(''pipeline_fxp_sweep'') first.'], matFile);
    end
    S      = load(matFile);
    tbl    = S.tbl;
    P      = S.params;
    SNR_dB = S.SNR_dB_vec(:).';
    ENOB   = sort(S.ENOB_vec(:).');
    ENOB   = ENOB(ENOB <= maxEnob);       % cap the ENOB axis (default <= 8)
    nEnob  = numel(ENOB);
    eMax   = max(ENOB);

    [EAdd, EMult] = energyCoefs(node);
    fprintf(['\n=== Pipeline fxp sweep post-processing ===\n', ...
        'FEC BER = %.0e | node = %s (E_A = %.2f*n fJ, E_M = %.2f*n^2 fJ)\n', ...
        'ADC = %.2f pJ/bit | highest ENOB = %d | M = %d\n'], ...
        fecBer, node, EAdd, EMult, eadcPJ, eMax, M);

    %% --- Per-combo (design row) metadata ---------------------------------
    designIdx = unique(tbl.design_idx(:)).';
    nDesign   = numel(designIdx);
    labels    = cell(nDesign, 1);
    shortLab  = cell(nDesign, 1);
    for d = 1:nDesign
        di  = designIdx(d);
        row = find(tbl.design_idx == di, 1);
        po2 = tbl.po2(row);
        fl  = [tbl.static_fl(row), tbl.clk_fl(row), tbl.adapt_fl(row), ...
               tbl.fr_fl(row),     tbl.cr_fl(row)];
        % eq/clk and carrier precision tiers are swept independently, so name
        % each group separately.
        eqTier = tierName(fl(1:3), [8 6 4],  [10 8 6]);
        crTier = tierName(fl(4:5), [10 4],   [12 6]);
        labels{d}   = sprintf('po2=%d, eq=%s, cr=%s [%d/%d/%d | %d/%d]', ...
            po2, eqTier, crTier, fl(1), fl(2), fl(3), fl(4), fl(5));
        shortLab{d} = sprintf('po2=%d\\newlineeq:%s cr:%s', po2, eqTier, crTier);
    end
    colors = lines(nDesign);

    %% --- RSNR for every (design, enob) cell ------------------------------
    rsnrMean = nan(nDesign, nEnob);
    rsnrStd  = nan(nDesign, nEnob);
    for d = 1:nDesign
        for e = 1:nEnob
            row = find(tbl.design_idx == designIdx(d) & ...
                       tbl.enob == ENOB(e), 1);
            if isempty(row), continue; end
            berMat  = tbl.ber{row};            % [NSNR x NTrials]
            nTrials = size(berMat, 2);
            perTr   = nan(nTrials, 1);
            for tr = 1:nTrials
                perTr(tr) = fecSnrFromBer(SNR_dB, berMat(:, tr).', fecBer);
            end
            ok = isfinite(perTr);
            if any(ok)
                rsnrMean(d, e) = mean(perTr(ok));
                rsnrStd(d, e)  = std(perTr(ok));
            end
        end
    end

    %% --- Figure 1: RSNR vs ENOB ------------------------------------------
    f1 = figure('Name', 'Pipeline: RSNR vs ENOB', ...
        'Position', [80 80 760 500]);
    ax1 = axes(f1); hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'on');
    for d = 1:nDesign
        if ~any(isfinite(rsnrMean(d, :)))
            fprintf('  skipping (never converges): %s\n', labels{d});
            continue;                       % combo never reaches FEC BER
        end
        errorbar(ax1, ENOB, rsnrMean(d, :), rsnrStd(d, :), ...
            '-o', 'Color', colors(d, :), 'LineWidth', 1.5, ...
            'MarkerSize', 6, 'MarkerFaceColor', colors(d, :), ...
            'CapSize', 8, 'DisplayName', labels{d});
    end
    xlabel(ax1, 'ADC ENOB (bits)');
    ylabel(ax1, sprintf('RSNR at BER = %.0e (dB)', fecBer));
    title(ax1, 'Required SNR vs ADC resolution');
    xticks(ax1, ENOB);
    legend(ax1, 'Location', 'northeast', 'Interpreter', 'none');
    if savePlots
        out = fullfile(here, 'pipeline_rsnr_vs_enob.png');
        exportgraphics(f1, out, 'Resolution', 200);
        fprintf('Saved %s\n', out);
    end

    %% --- Per-combo energy breakdown at the highest ENOB ------------------
    %  Common (combo-independent) blocks: ADC + GSOP + deskew, all at eMax.
    eta   = P.SpS;
    E_adc = eadcPJ;                        % flat ADC cost: pJ per bit
    [NM_g, NA_g] = gsopOps(eta);
    E_gsop  = 1e-3 * energy.receiver(NA_g, NM_g, EAdd, EMult, M, 1, eMax);
    N_DS = deskewOrder + 1;
    [NM_d, NA_d] = deskewOps(eta, N_DS);
    E_deskew = 1e-3 * energy.receiver(NA_d, NM_d, EAdd, EMult, M, 1, eMax);

    N    = P.NFFT;
    NCD  = P.NCD;
    beta = P.Rolloff;
    NTaps = P.NTapsAEQ;
    L    = P.TrainingLen;
    BlkCR = P.BlockLen_CR;
    SF   = P.SUBFRAME_SYMS;

    E_static = nan(nDesign, 1);
    E_clk    = nan(nDesign, 1);
    E_aeq    = nan(nDesign, 1);
    E_fr     = nan(nDesign, 1);
    E_cr     = nan(nDesign, 1);
    for d = 1:nDesign
        row = find(tbl.design_idx == designIdx(d) & tbl.enob == eMax, 1);
        po2 = tbl.po2(row);
        nS  = tbl.static_wl(row);
        nC  = tbl.clk_wl(row);
        nA  = tbl.adapt_wl(row);
        nF  = tbl.fr_wl(row);
        nR  = tbl.cr_wl(row);

        [NM_s, NA_s] = staticOps(N, NCD, eta, po2);
        [NM_c, NA_c] = godardOps(N, NCD, eta, beta);
        [NM_a, NA_a] = adaptOpsSignSign(NTaps);
        [NM_f, NA_f] = frOps(L, SF);
        % Carrier recovery = pilots-only estimate (per block) + the per-symbol
        % CORDIC phase rotation applied to every symbol.
        [NM_cr_e, NA_cr_e] = crEstOps(BlkCR);

        E_static(d) = 1e-3 * energy.receiver(NA_s, NM_s, EAdd, EMult, M, 1, nS);
        E_clk(d)    = 1e-3 * energy.receiver(NA_c, NM_c, EAdd, EMult, M, 1, nC);
        E_aeq(d)    = 1e-3 * energy.receiver(NA_a, NM_a, EAdd, EMult, M, 1, nA);
        E_fr(d)     = 1e-3 * energy.receiver(NA_f, NM_f, EAdd, EMult, M, 1, nF);
        E_cr(d)     = 1e-3 * ( ...
            energy.receiver(NA_cr_e, NM_cr_e, EAdd, EMult, M, 1, nR) + ...
            energy.receiver(0,       1,       EAdd, EMult, M, 1, nR));  % apply
    end
    E_front = E_gsop + E_deskew;                  % per-combo identical
    E_total = E_adc + E_front + E_static + E_clk + E_aeq + E_fr + E_cr;

    %% --- Printed energy table (pJ/bit, highest ENOB) ---------------------
    fprintf(['\n--- Energy breakdown per combo at ENOB = %d ', ...
        '(pJ/bit, %s) ---\n'], eMax, node);
    fprintf('%-26s %7s %7s %7s %7s %7s %7s %7s %8s\n', ...
        'combo', 'ADC', 'GSOP', 'deskew', 'static', 'clk', 'aeq', ...
        'fr+cr', 'TOTAL');
    for d = 1:nDesign
        fprintf('%-26s %7.2f %7.3f %7.3f %7.2f %7.2f %7.3f %7.3f %8.2f\n', ...
            labels{d}, E_adc, E_gsop, E_deskew, E_static(d), E_clk(d), ...
            E_aeq(d), E_fr(d) + E_cr(d), E_total(d));
    end

    %% --- Figure 2: relative RSNR vs relative energy (highest ENOB) -------
    eIdx     = find(ENOB == eMax, 1);
    rsnrAtMax = rsnrMean(:, eIdx);
    f2 = figure('Name', 'Pipeline: RSNR vs receiver energy', ...
        'Position', [100 100 760 520]);
    ax2 = axes(f2); hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');
    for d = 1:nDesign
        if ~isfinite(rsnrAtMax(d))
            continue;                       % non-converging combo — skip
        end
        plot(ax2, E_total(d), rsnrAtMax(d), 'o', 'MarkerSize', 9, ...
            'MarkerFaceColor', colors(d, :), 'MarkerEdgeColor', 'k', ...
            'DisplayName', labels{d});
    end
    xlabel(ax2, 'Total receiver energy (pJ/bit)');
    ylabel(ax2, sprintf('RSNR at BER = %.0e (dB)', fecBer));
    title(ax2, sprintf(['Energy-vs-performance trade-off ', ...
        '(ENOB = %d, BER = %.0e)'], eMax, fecBer));
    legend(ax2, 'Location', 'best', 'Interpreter', 'none');
    if savePlots
        out = fullfile(here, 'pipeline_rel_rsnr_energy.png');
        exportgraphics(f2, out, 'Resolution', 200);
        fprintf('Saved %s\n', out);
    end

    %% --- Figure 3: grouped stacked energy breakdown bar (highest ENOB) ---
    %  Segments grouped into four functional blocks for clarity:
    %    Front end       = ADC + GSOP + deskew  (combo-independent)
    %    Equalisation    = static eq + adaptive eq
    %    Clock recovery  = Godard
    %    Carrier recovery= frequency recovery + phase (pilots-only) recovery
    E_frontGrp   = E_adc + E_gsop + E_deskew;   % scalar
    E_equalGrp   = E_static + E_aeq;            % [nDesign x 1]
    E_clkGrp     = E_clk;
    E_carrierGrp = E_fr + E_cr;
    segs = [repmat(E_frontGrp, nDesign, 1), E_equalGrp, E_clkGrp, E_carrierGrp];
    segNames = {'Front end (ADC + GSOP + deskew)', ...
                'Equalisation (static + adaptive)', ...
                'Clock recovery', ...
                'Carrier recovery (freq + phase)'};

    keep    = isfinite(rsnrAtMax);          % only converging combos
    segsK   = segs(keep, :);
    labK    = shortLab(keep);
    nKeep   = nnz(keep);

    if nKeep == 0
        warning('process_pipeline_fxp_sweep:noConverge', ...
            'No combo converges at ENOB = %d; skipping energy bar.', eMax);
    else
        % A single group as a row vector would be drawn as separate bars, so
        % pad with an invisible NaN group to force the stacked-matrix path.
        if nKeep == 1
            plotSegs = [segsK; nan(1, size(segsK, 2))];
        else
            plotSegs = segsK;
        end
        f3 = figure('Name', 'Pipeline: energy breakdown', ...
            'Position', [120 120 820 540]);
        ax3 = axes(f3);
        hb = bar(ax3, plotSegs, 'stacked');
        for k = 1:numel(hb)
            hb(k).DisplayName = segNames{k};
        end
        grid(ax3, 'on'); box(ax3, 'on');
        xlim(ax3, [0.5, nKeep + 0.5]);
        xticks(ax3, 1:nKeep);
        xticklabels(ax3, labK);
        ylabel(ax3, 'Energy per bit (pJ)');
        title(ax3, sprintf('Receiver energy breakdown (ENOB = %d, %s)', ...
            eMax, node));
        legend(ax3, 'Location', 'eastoutside');
        if savePlots
            out = fullfile(here, 'pipeline_energy_breakdown.png');
            exportgraphics(f3, out, 'Resolution', 200);
            fprintf('Saved %s\n', out);
        end
    end
end


% =====================================================================
%  Labelling
% =====================================================================
function t = tierName(vec, lowRef, highRef)
%TIERNAME  Friendly 'low' / 'high' name for a precision sub-vector, matched
%   against the two report tiers; falls back to the raw values otherwise.
    if isequal(vec(:).', lowRef(:).')
        t = 'low';
    elseif isequal(vec(:).', highRef(:).')
        t = 'high';
    else
        t = mat2str(vec(:).');
    end
end


% =====================================================================
%  Energy / operation-count helpers (report/full/full.tex tables)
% =====================================================================
function [EAdd, EMult] = energyCoefs(node)
%ENERGYCOEFS  Linear (per-bit) and quadratic (per-bit^2) energy coefs [fJ]
%   from report/full/full.tex.  '45nm' is the Horowitz 2014 fit; '14nm'
%   applies the 75% process-scaling estimate.
    switch lower(node)
        case '45nm'
            EAdd = 3.16;  EMult = 3.03;
        case '14nm'
            EAdd = 0.79;  EMult = 0.76;
        otherwise
            error('process_pipeline_fxp_sweep:badNode', ...
                'Unknown Node ''%s'' (use ''45nm'' or ''14nm'').', node);
    end
end


function [NMult, NAdd] = staticOps(N, NCD, eta, po2)
%STATICOPS  Per-symbol RM/RA for the overlap-save CD + matched filter.
%   tab:cd_cost.  po2 replaces every FFT twiddle multiply with bit shifts.
    denom = N - NCD + 1;
    if po2
        NMult = eta * 4*N / denom;
    else
        NMult = eta * (4*N*log2(N) + 4*N) / denom;
    end
    NAdd = eta * (6*N*log2(N) + 2*N) / denom;
end


function [NMult, NAdd] = godardOps(N, NCD, eta, beta)
%GODARDOPS  Per-symbol RM/RA for the Modified-Godard metric, CORDIC arg(S),
%   and FD phase ramp.  tab:clk_cost (the FFT/IFFT is shared with the static
%   stage and counted there).
    denom = N - NCD + 1;
    NMult = 4*beta*N/denom + eta/denom + 8*N*eta/denom;
    NAdd  = 4*beta*N/denom + 0          + 4*N*eta/denom;
end


function [NMult, NAdd] = adaptOpsSignSign(NTaps)
%ADAPTOPSSIGNSIGN  Per-symbol RM/RA for the sign-sign butterfly CMA.
%   tab:aeq_cost (FIR outputs + error terms + multiplier-free update).
    NMult = 8*NTaps + 2;
    NAdd  = (4*NTaps + 2) + 1 + 8*NTaps;
end


function [NMult, NAdd] = frOps(L, subframeLen)
%FROPS  Per-symbol RM/RA for differential-phase-and-Kay frequency recovery.
%   tab:diffkay_cost gives 7L+1 RM and 4L-2 RA per subframe (L = training
%   length); a single per-subframe estimate is amortised over the subframe.
    NMult = (7*L + 1) / subframeLen;
    NAdd  = (4*L - 2) / subframeLen;
end


function [NMult, NAdd] = crEstOps(blockLen)
%CRESTOPS  Per-symbol RM/RA for the pilots-only carrier-recovery estimate.
%   tab:pilot_cost: one CORDIC (1 RM + 1 RA) per block, amortised over the
%   block.  The per-symbol phase rotation is charged separately by the caller.
    NMult = 1 / blockLen;
    NAdd  = 1 / blockLen;
end


function [NMult, NAdd] = gsopOps(eta)
%GSOPOPS  Per-symbol (per-polarisation) RM/RA for GSOP.  tab:gsop_cost:
%   3*eta RM and eta RA; block statistics are amortised to ~0.
    NMult = 3 * eta;
    NAdd  = eta;
end


function [NMult, NAdd] = deskewOps(eta, N_DS)
%DESKEWOPS  Per-symbol (per-polarisation) RM/RA for the deskew interpolator.
%   tab:deskew_cost: two real N_DS-tap FIRs (I and Q) at oversampling eta.
    NMult = 2 * eta * N_DS;
    NAdd  = 2 * eta * (N_DS - 1);
end


% =====================================================================
function snr = fecSnrFromBer(snrVec, berVec, fecBer)
%FECSNRFROMBER  Log-linear interpolation of the BER curve to find the SNR at
%   which BER == fecBer.  Returns NaN if the curve does not bracket fecBer.
    v = berVec(:).';
    v(v <= 0) = NaN;
    valid = isfinite(v);
    if nnz(valid) < 2
        snr = NaN; return;
    end
    x = snrVec(valid);
    y = log10(v(valid));
    target = log10(fecBer);
    if target < min(y) || target > max(y)
        snr = NaN; return;
    end
    for k = 1:numel(x) - 1
        if (y(k) - target) * (y(k+1) - target) <= 0
            snr = x(k) + (target - y(k)) / (y(k+1) - y(k)) * (x(k+1) - x(k));
            return;
        end
    end
    snr = NaN;
end
