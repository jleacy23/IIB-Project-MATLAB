function diag_pipeline_float(LW_override, CFO_override)
%DIAG_PIPELINE_FLOAT  Floating-point clone of pipeline_fxp_sweep.runOneTrial
%   (current Godard pipeline) with per-stage instrumentation, so the BER
%   failure can be localised to STRUCTURE (equaliser convergence /
%   alignment / FR / CR) rather than fixed-point precision.
%
%   Mirrors the current pipeline:
%       CPON Tx -> channel -> normalise -> Godard-combined eq (float) ->
%       pilot lag-align -> differential_kay FR (float) -> pilots_only CR
%       (float, per-pol) -> BER.   No ADC, no pilot rescaling.
%
%   The decisive readout is the verdict block at the end, which compares:
%     (a) ORACLE BER  - eqOut with EXACT residual-CFO removed + per-pol
%                       pilot CR (isolates equaliser+alignment quality), and
%     (b) PIPELINE BER- eqOut -> differential_kay -> pilots_only.
%   Their relationship pinpoints the broken stage.
%
%   Run:  addpath(genpath('src')); diag_pipeline_float

    here     = fileparts(mfilename('fullpath'));
    repoRoot = fileparts(fileparts(here));
    addpath(genpath(fullfile(repoRoot, 'src')));

    %% ---- Parameters (mirror pipeline_fxp_sweep) --------------------
    P.Rs=30.5; P.SpS=2; P.N_pol=2; P.D=20; P.CWL=1550; P.DGDSpec=0.1;
    P.N_pmd=1; P.PMD_seed=42; P.L_km=80; P.SFO_ppm=40; P.CFO_GHz=3;
    P.LW_Hz=1e6; P.Rolloff=0.25; P.Span=10; P.N_sub_target=8;
    P.NFFT=128; P.NCD=22;
    P.ki=3e-1; P.kp=1e-1;                  % Godard po2-off gains
    P.kiGardner=1e-7; P.kpGardner=1e-6; P.NLanesGard=32;   % Gardner po2-off
    P.NTapsAEQ=1; P.MuAEQ=1e-3; P.N1AEQ=500; P.NOutAEQ=1000;
    P.SignOnly=true; P.SingleSpike=true; P.PLanesAEQ=32;
    P.MaxFreq=0.1; P.TrainingLen=11; P.BlockLen_CR=32;
    P.NormPct=99.9; P.AlignLagMax=64; P.SUBFRAME_SYMS=3712;

    if nargin>=1 && ~isempty(LW_override),  P.LW_Hz  = LW_override;  end
    if nargin>=2 && ~isempty(CFO_override), P.CFO_GHz = CFO_override; end

    SNR_dB = 20; trialSeed = 1;
    fprintf('\n============== FLOAT GODARD PIPELINE DIAGNOSTIC ==============\n');
    fprintf('SNR=%g dB  CFO=%g GHz  LW=%g Hz  trial=%d\n', ...
        SNR_dB, P.CFO_GHz, P.LW_Hz, trialSeed);

    %% ---- 1. Channel (mirror genChannel incl. normalise) ------------
    rng(1000*trialSeed + round(SNR_dB) + 13);
    CPON_BITS_PER_SF = 3586 * P.N_pol * 2;
    nBits = P.N_sub_target * CPON_BITS_PER_SF;
    bits  = randi([0 1], nBits, 1);
    [txSymbols, pilotsRef, training, nSubTx] = modem.modulate(bits);
    NsymTx = size(txSymbols,1);
    fprintf('\n[1] Tx: %d symbols, %d subframes. pilotsRef=[%dx%d] training=[%dx%d]\n', ...
        NsymTx, nSubTx, size(pilotsRef,1), size(pilotsRef,2), ...
        size(training,1), size(training,2));
    checkTxPilotPlacement(txSymbols, pilotsRef, training, P);

    txSig = modem.rrcPulse(txSymbols, P.SpS, P.Rolloff, P.Span);
    rxSig = channel.add_chromatic_dispersion(txSig, P.L_km, P.SpS, P.Rs, P.D, P.CWL);
    s = rng; rng(P.PMD_seed);
    rxSig = channel.add_pmd(rxSig, P.L_km, P.SpS, P.Rs, P.DGDSpec, P.N_pmd);
    rng(s);
    rxSig = channel.lo_freq_shift(rxSig, P.CFO_GHz*1000, P.Rs, P.SpS);
    rxSig = channel.apply_timing_error(rxSig, P.SFO_ppm, 0, P.SpS);
    rxSig = channel.add_awgn(rxSig, SNR_dB);
    rxSig = channel.add_phase_noise(rxSig, P.Rs*P.SpS, P.LW_Hz);
    % rxSig = modem.normalise(rxSig, P.NormPct);   % no ADC

    %% ---- 2. Equaliser (FLOAT Godard) ------------------------------
    nOv = 2*ceil((P.NCD-1)/2);
    adaptOpts = struct('NTaps',P.NTapsAEQ,'Mu',P.MuAEQ,'SingleSpike',P.SingleSpike, ...
        'N1',P.N1AEQ,'NOut',P.NOutAEQ,'SignOnly',P.SignOnly,'PLanes',P.PLanesAEQ, ...
        'Mode',0,'Pilots',[],'BlockLen',P.PLanesAEQ,'SubframeBlocks',0);
    [eqOutFull, cfoBins] = eq_clk.combined_cd_fd_godard_adaptive(rxSig, P.SpS, ...
        P.NFFT, nOv, P.D, P.L_km, P.CWL, P.Rs, P.Rolloff, P.ki, P.kp, ...
        NsymTx, adaptOpts, true, false);
    binGHz = P.SpS*P.Rs/P.NFFT;
    fprintf('\n[2] Equaliser out: %d symbols. cfoBins=%.4f (bin=%.4f GHz -> %.4f GHz removed; residual target %.4f GHz)\n', ...
        size(eqOutFull,1), cfoBins, binGHz, cfoBins*binGHz, P.CFO_GHz - cfoBins*binGHz);

    %% ---- 3. Pilot lag-align (same as pipeline) --------------------
    nDropNom = P.SUBFRAME_SYMS - P.NOutAEQ;
    fprintf('\n[3] PILOT-COHERENCE vs lag (nominal drop=%d, window +/-%d):\n', ...
        nDropNom, P.AlignLagMax);
    [bestLag, Ccurve, lags] = scanCoherence(eqOutFull, nDropNom, pilotsRef, P, P.AlignLagMax);
    showLagRow(lags, Ccurve, bestLag);
    % Also test the polarisation-SWAPPED hypothesis (butterfly CMA can swap
    % pols): if swapped coherence is much higher, the pol assignment is wrong.
    pilotsSwap = pilotsRef(:, [2 1]);
    [bestLagSw, CcurveSw] = scanCoherence(eqOutFull, nDropNom, pilotsSwap, P, P.AlignLagMax);
    fprintf('     max coherence: direct=%.3f (lag %+d) | pol-swapped=%.3f (lag %+d)\n', ...
        max(Ccurve), bestLag, max(CcurveSw), bestLagSw);
    if max(CcurveSw) > max(Ccurve) + 0.05
        fprintf('     *** pol-swapped coherence higher => equaliser likely SWAPPED X/Y\n');
    end

    [eqOut, nSub, swapped] = pipeline_fxp_sweep.alignToSubframe(eqOutFull, nDropNom, pilotsRef, P);
    fprintf('     => alignToSubframe chose: swapped=%d (direct-only scan best lag was %+d), %d subframes (%d symbols)\n', ...
        swapped, bestLag, nSub, size(eqOut,1));
    if nSub < 1, fprintf('     (no usable subframes) ABORT\n'); return; end

    txOffset = P.SUBFRAME_SYMS;   % aligned stream starts at Tx subframe 2
    pilotsAll = repmat(pilotsRef, nSub, 1);

    %% ---- 4. ORACLE: exact residual-CFO removal + per-pol pilot CR --
    % Isolates EQUALISER + ALIGNMENT quality (independent of FR).  Tries the
    % direct and the X/Y-swapped pilot assignment: nothing in the current
    % pipeline resolves the butterfly-CMA polarisation-swap ambiguity, so a
    % LOW swapped-oracle with a HIGH direct-oracle pins the bug on pol swap.
    residGHz = P.CFO_GHz - cfoBins*binGHz;
    eqCfo    = removeResidCFO(eqOut, residGHz, P.Rs);
    crOracle = carrier_recovery.pilots_only(eqCfo, P.N_pol, P.BlockLen_CR, pilotsAll);
    berOracle = berAtOffset(crOracle, txSymbols, txOffset);
    pilotsAllSw = repmat(pilotsRef(:,[2 1]), nSub, 1);
    crOracleSw  = carrier_recovery.pilots_only(eqCfo, P.N_pol, P.BlockLen_CR, pilotsAllSw);
    berOracleSw = berAtOffset(crOracleSw, txSymbols, txOffset);
    fprintf('\n[4] ORACLE BER (exact CFO removed + per-pol pilot CR):\n');
    fprintf('      direct pilots = %.4e | X/Y-swapped pilots = %.4e\n', ...
        berOracle, berOracleSw);
    if berOracleSw < berOracle - 1e-2
        fprintf('      *** swapped pilots win => EQUALISER SWAPPED X/Y (pol-swap unresolved)\n');
    end
    berOracleBest = min(berOracle, berOracleSw);
    % Alignment cross-check: scan a WIDE symbol offset on the oracle path.
    fprintf('    oracle BER vs extra symbol offset around the chosen boundary:\n');
    for off = [-2 -1 0 1 2 32 -32]
        e = oracleBerAtLag(eqOutFull, nDropNom+off, cfoBins, binGHz, pilotsRef, P, txSymbols, txOffset);
        fprintf('      offset %+4d : BER=%.4e %s\n', off, e, tern(off==0,'<== chosen',''));
    end

    %% ---- 5. PIPELINE: differential_kay FR -> pilots_only CR --------
    [frOut, frHz] = freq_recovery.differential_kay(eqOut, training, P.Rs, true, 0);
    fprintf('\n[5] FR: estimate=%.4f GHz  (eq+FR removed=%.4f GHz; channel CFO=%.4f GHz)\n', ...
        frHz/1e9, cfoBins*binGHz + frHz/1e9, P.CFO_GHz);
    reportPilotPhaseSlope(eqCfo, pilotsAll, P, 'oracle (CFO-removed)');
    reportPilotPhaseSlope(frOut,  pilotsAll, P, 'after FR');

    crOut = carrier_recovery.pilots_only(frOut, P.N_pol, P.BlockLen_CR, pilotsAll);
    berPipe = berAtOffset(crOut, txSymbols, txOffset);
    fprintf('\n[6] PIPELINE BER (FR then per-pol pilot CR) = %.4e\n', berPipe);

    % Per-pol pilot phase spread after FR (large spread X-vs-Y => pol issue).
    [stdX, stdY, stdXY] = perPolPhaseSpread(frOut, pilotsAll, P);
    fprintf('    per-pol pilot phase std after FR: X=%.3f Y=%.3f  X-Y=%.3f rad\n', ...
        stdX, stdY, stdXY);

    %% ---- Constellation at each stage -----------------------------
    %  Visual check: after a clean recovery each pol should show 4 tight
    %  QPSK clusters.  A smeared RING => carrier phase not tracked; a filled
    %  BLOB => equaliser did not open the eye; clusters at the wrong angle =>
    %  alignment / pol-swap.
    stages = { eqOut,    'after EQ (aligned)'; ...
               eqCfo,    'after exact-CFO removal'; ...
               crOracle, sprintf('oracle CR (BER=%.2e)', berOracle); ...
               frOut,    sprintf('after FR (est=%.3f GHz)', frHz/1e9); ...
               crOut,    sprintf('pipeline CR (BER=%.2e)', berPipe) };
    plotConstellations(stages, P);

    %% ---- IMPAIRMENT / STAGE ISOLATION ----------------------------
    %  Re-run channel+eq+align+exact-CFO+per-pol CR with one thing toggled
    %  from the CURRENT P.  Whichever toggle drops the oracle BER is the
    %  cause of the closed eye.  'CMA frozen' (Mu=0) separates the adaptive
    %  equaliser from the static-CD+Godard-timing front end: if frozen is
    %  ALSO a blob the fault is CD/timing, if frozen is clean it is the CMA.
    fprintf('\n[ISO] oracle BER (eq+align+exact-CFO+per-pol CR), toggled from current P:\n');
    Pa = P;                                    Pa.tag = 'baseline (current P)';
    Pb = P; Pb.PmdOn = false;                  Pb.tag = 'PMD off';
    Pc = P; Pc.LW_Hz  = 0;                     Pc.tag = 'phase noise off (LW=0)';
    Pd = P; Pd.MuAEQ  = 0;                     Pd.tag = 'CMA frozen (Mu=0)';
    Pe = P; Pe.PmdOn = false; Pe.LW_Hz = 0;    Pe.tag = 'PMD off + LW=0';
    for cell_ = {Pa, Pb, Pc, Pd, Pe}
        Pi = cell_{1};
        e  = condBer(Pi, SNR_dB, trialSeed);
        fprintf('     %-26s : oracle BER = %.4e\n', Pi.tag, e);
    end

    %% ---- GARDNER cross-check (same channel realisation) ----------
    %  Swap ONLY the clock-recovery flavour (Godard -> Gardner) on the very
    %  same rxSig.  If Gardner is clean but Godard is not, the GODARD timing
    %  loop is the culprit; if both blob, the fault is shared (CMA / CD-FD /
    %  phase noise / channel), not the timing flavour.
    fprintf('\n[G] Gardner vs Godard on the SAME channel (oracle BER):\n');
    [berG, eqCfoG] = oracleFromRx(P, rxSig, NsymTx, txSymbols, pilotsRef, 'gardner');
    fprintf('     Godard oracle BER = %.4e | Gardner oracle BER = %.4e\n', berOracle, berG);
    if isfinite(berG) && berG < berOracle - 1e-2
        fprintf('     *** Gardner clean, Godard not => GODARD timing recovery is the culprit\n');
    elseif isfinite(berG) && berG > 1e-2
        fprintf('     *** both fail => NOT the clock-recovery flavour (CMA / CD-FD / phase noise)\n');
    end
    plotConstellations({ eqCfo,  'GODARD after exact-CFO'; ...
                         eqCfoG, 'GARDNER after exact-CFO' }, P);

    %% ---- 7. VERDICT ----------------------------------------------
    fprintf('\n================ VERDICT =================================\n');
    if berOracleSw < berOracle - 1e-2 && berOracleSw < 1e-2
        fprintf('  POL SWAP: direct-oracle BER=%.2e is high but X/Y-swapped\n', berOracle);
        fprintf('  oracle BER=%.2e is low => the butterfly CMA swapped the\n', berOracleSw);
        fprintf('  polarisations and NOTHING downstream resolves it.  The pilot\n');
        fprintf('  align + CR use pilotsRef in fixed X/Y order, so a swap floors\n');
        fprintf('  the BER.  FIX: make the pilot assignment swap-aware (try both\n');
        fprintf('  pilot column orders, keep the higher-coherence one for BOTH\n');
        fprintf('  alignToSubframe and the CR) -- the training sequence (distinct\n');
        fprintf('  per pol) exists precisely to resolve this.\n');
    elseif berOracleBest > 1e-2
        fprintf('  ORACLE BER is HIGH (direct=%.2e, swapped=%.2e) => the\n', berOracle, berOracleSw);
        fprintf('  EQUALISER output or the ALIGNMENT is broken (not FR/CR).\n');
        fprintf('  - Check [4] offset scan: if a non-zero offset gives low BER,\n');
        fprintf('    the lag-align picked the wrong boundary.\n');
        fprintf('  - Otherwise the single-tap CMA did not converge (NTaps=%d,\n', P.NTapsAEQ);
        fprintf('    SingleSpike) against the residual CD/PMD.\n');
    elseif berPipe > 1e-2
        fprintf('  ORACLE BER is LOW (%.2e) but PIPELINE BER is HIGH (%.2e)\n', berOracleBest, berPipe);
        fprintf('  => equaliser+alignment are FINE; the FREQUENCY RECOVERY\n');
        fprintf('  (differential_kay) is the culprit (over/under-correction or\n');
        fprintf('  a wrap near MaxFreq=%.2f).  Compare [5] FR estimate vs CFO\n', P.MaxFreq);
        fprintf('  and the pilot-phase slope before/after FR.\n');
    else
        fprintf('  BOTH oracle (%.2e) and pipeline (%.2e) BER are LOW.\n', berOracleBest, berPipe);
        fprintf('  Float structure is healthy => the fxp failure is PRECISION:\n');
        fprintf('  sweep FRFL/CRFL/AdaptFL upward and check CordicIts=FL is\n');
        fprintf('  large enough (angular resolution ~atan(2^-FL)).\n');
    end
    fprintf('=========================================================\n');
end

% =====================================================================
%  Helpers
% =====================================================================
function checkTxPilotPlacement(tx, pilotsRef, training, P)
    base = 0;
    okTrain = isequal(tx(base+(1:P.TrainingLen),:), training);
    okP = true;
    for b = 2:size(pilotsRef,1)
        if ~isequal(tx(base+(b-1)*P.BlockLen_CR+1,:), pilotsRef(b,:)), okP=false; break; end
    end
    okP1 = isequal(tx(base+1,:), pilotsRef(1,:));
    fprintf('     Tx placement: training@1..%d=%d, pilot@blockstarts=%d, TS1==pilot(1)=%d\n', ...
        P.TrainingLen, okTrain, okP, okP1);
end

function [bestLag, C, lags] = scanCoherence(eqOutFull, nDropNom, pilotsRef, P, win)
    lags = -win:win;
    C = nan(size(lags));
    N = size(eqOutFull,1);
    for i = 1:numel(lags)
        start = nDropNom + 1 + lags(i);
        if start < 1, continue; end
        nSub = floor((N - start + 1) / P.SUBFRAME_SYMS);
        if nSub < 1, continue; end
        C(i) = pipeline_fxp_sweep.pilotCoherence(eqOutFull, start, nSub, pilotsRef, P);
    end
    [~, k] = max(C);
    bestLag = lags(k);
end

function showLagRow(lags, C, bestLag)
    pick = [-2 -1 0 1 2 bestLag];
    for lg = unique(pick, 'stable')
        i = find(lags==lg, 1);
        if isempty(i) || isnan(C(i)), continue; end
        fprintf('     lag=%+3d : coherence=%.4f %s\n', lg, C(i), tern(lg==bestLag,'<== best',''));
    end
end

function y = removeResidCFO(x, residGHz, Rs)
    n = (0:size(x,1)-1).';
    y = x .* exp(-1j*2*pi*(residGHz/Rs)*n);
end

function e = oracleBerAtLag(eqOutFull, nDrop, cfoBins, binGHz, pilotsRef, P, txSymbols, txOffset)
    N = size(eqOutFull,1);
    start = nDrop + 1;
    if start < 1, e = NaN; return; end
    nSub = floor((N - start + 1) / P.SUBFRAME_SYMS);
    if nSub < 1, e = NaN; return; end
    seg = eqOutFull(start : start + nSub*P.SUBFRAME_SYMS - 1, :);
    seg = removeResidCFO(seg, P.CFO_GHz - cfoBins*binGHz, P.Rs);
    pilotsAll = repmat(pilotsRef, nSub, 1);
    cr = carrier_recovery.pilots_only(seg, P.N_pol, P.BlockLen_CR, pilotsAll);
    e = berAtOffset(cr, txSymbols, txOffset);
end

function reportPilotPhaseSlope(x, pilotsAll, P, tag)
    NB = size(pilotsAll,1);
    th = zeros(NB,1);
    for b = 1:NB
        bs = (b-1)*P.BlockLen_CR + 1;
        if bs <= size(x,1)
            th(b) = angle(sum(conj(pilotsAll(b,:)) .* x(bs,:)));
        end
    end
    dth = mod(diff(th)+pi,2*pi)-pi;
    fprintf('     [%s] pilot-phase step block-to-block: mean=%.4f std=%.4f rad/blk\n', ...
        tag, mean(dth), std(dth));
end

function [stdX, stdY, stdXY] = perPolPhaseSpread(x, pilotsAll, P)
    NB = size(pilotsAll,1);
    thX = nan(NB,1); thY = nan(NB,1);
    for b = 1:NB
        bs = (b-1)*P.BlockLen_CR + 1;
        if bs <= size(x,1)
            thX(b) = angle(conj(pilotsAll(b,1))*x(bs,1));
            thY(b) = angle(conj(pilotsAll(b,2))*x(bs,2));
        end
    end
    dXY = mod(thX-thY+pi,2*pi)-pi;
    stdX = std(thX,'omitnan'); stdY = std(thY,'omitnan'); stdXY = std(dXY,'omitnan');
end

function [berO, eqCfo] = oracleFromRx(P, rxSig, NsymTx, txSymbols, pilotsRef, blockName)
    % Run the chosen combined eq block on a GIVEN rxSig, then align +
    % exact-CFO removal + per-pol pilot CR -> oracle BER + CFO-removed const.
    nOv = 2*ceil((P.NCD-1)/2);
    ao = struct('NTaps',P.NTapsAEQ,'Mu',P.MuAEQ,'SingleSpike',P.SingleSpike, ...
        'N1',P.N1AEQ,'NOut',P.NOutAEQ,'SignOnly',P.SignOnly,'PLanes',P.PLanesAEQ, ...
        'Mode',0,'Pilots',[],'BlockLen',P.PLanesAEQ,'SubframeBlocks',0);
    switch blockName
        case 'godard'
            [eqF, cfoBins] = eq_clk.combined_cd_fd_godard_adaptive(rxSig, P.SpS, ...
                P.NFFT, nOv, P.D, P.L_km, P.CWL, P.Rs, P.Rolloff, ...
                P.ki, P.kp, NsymTx, ao, true, false);
        case 'gardner'
            [eqF, cfoBins] = eq_clk.combined_cd_fd_gardner_adaptive(rxSig, P.SpS, ...
                P.NFFT, nOv, P.D, P.L_km, P.CWL, P.Rs, P.Rolloff, ...
                P.kiGardner, P.kpGardner, NsymTx, P.NLanesGard, ao, true, false);
        otherwise
            error('oracleFromRx:block', 'unknown block %s', blockName);
    end
    nDropNom = P.SUBFRAME_SYMS - P.NOutAEQ;
    [eqOut, nSub] = pipeline_fxp_sweep.alignToSubframe(eqF, nDropNom, pilotsRef, P);
    if nSub < 1, berO = NaN; eqCfo = complex(zeros(0,P.N_pol)); return; end
    binGHz = P.SpS*P.Rs/P.NFFT;
    eqCfo  = removeResidCFO(eqOut, P.CFO_GHz - cfoBins*binGHz, P.Rs);
    cr = carrier_recovery.pilots_only(eqCfo, P.N_pol, P.BlockLen_CR, repmat(pilotsRef, nSub, 1));
    berO = berAtOffset(cr, txSymbols, P.SUBFRAME_SYMS);
end

function e = condBer(P, SNR_dB, trialSeed)
    % Compact channel+eq+align+exact-CFO+per-pol-CR oracle BER for one P.
    % Honours P.PmdOn (default true) and whatever impairment levels P carries
    % (CFO/SFO/CD/LW); a 0 level skips that channel stage.
    rng(1000*trialSeed + round(SNR_dB) + 13);
    bits = randi([0 1], P.N_sub_target*3586*P.N_pol*2, 1);
    [tx, pilotsRef, ~, ~] = modem.modulate(bits);
    NsymTx = size(tx,1);
    sig = modem.rrcPulse(tx, P.SpS, P.Rolloff, P.Span);
    sig = channel.add_chromatic_dispersion(sig, P.L_km, P.SpS, P.Rs, P.D, P.CWL);
    if ~isfield(P,'PmdOn') || P.PmdOn
        st = rng; rng(P.PMD_seed);
        sig = channel.add_pmd(sig, P.L_km, P.SpS, P.Rs, P.DGDSpec, P.N_pmd);
        rng(st);
    end
    if P.CFO_GHz ~= 0, sig = channel.lo_freq_shift(sig, P.CFO_GHz*1000, P.Rs, P.SpS); end
    if P.SFO_ppm ~= 0, sig = channel.apply_timing_error(sig, P.SFO_ppm, 0, P.SpS); end
    sig = channel.add_awgn(sig, SNR_dB);
    if P.LW_Hz ~= 0,   sig = channel.add_phase_noise(sig, P.Rs*P.SpS, P.LW_Hz); end

    nOv = 2*ceil((P.NCD-1)/2);
    ao = struct('NTaps',P.NTapsAEQ,'Mu',P.MuAEQ,'SingleSpike',P.SingleSpike, ...
        'N1',P.N1AEQ,'NOut',P.NOutAEQ,'SignOnly',P.SignOnly,'PLanes',P.PLanesAEQ, ...
        'Mode',0,'Pilots',[],'BlockLen',P.PLanesAEQ,'SubframeBlocks',0);
    [eqF, cfoBins] = eq_clk.combined_cd_fd_godard_adaptive(sig, P.SpS, P.NFFT, nOv, ...
        P.D, P.L_km, P.CWL, P.Rs, P.Rolloff, P.ki, P.kp, NsymTx, ao, true, false);

    nDropNom = P.SUBFRAME_SYMS - P.NOutAEQ;
    [eqOut, nSub] = pipeline_fxp_sweep.alignToSubframe(eqF, nDropNom, pilotsRef, P);
    if nSub < 1, e = NaN; return; end
    binGHz = P.SpS*P.Rs/P.NFFT;
    eqCfo  = removeResidCFO(eqOut, P.CFO_GHz - cfoBins*binGHz, P.Rs);
    cr = carrier_recovery.pilots_only(eqCfo, P.N_pol, P.BlockLen_CR, repmat(pilotsRef, nSub, 1));
    e = berAtOffset(cr, tx, P.SUBFRAME_SYMS);
end

function plotConstellations(stages, P)
    % One figure: rows = pipeline stages, columns = polarisations.
    nS = size(stages, 1);
    fig = figure('Name', 'Float pipeline constellations', ...
        'Position', [60 60 340*P.N_pol 220*nS]);
    maxPts = 4000;                       % subsample for a readable scatter
    for si = 1:nS
        x   = stages{si, 1};
        lab = stages{si, 2};
        for pol = 1:P.N_pol
            ax = subplot(nS, P.N_pol, (si-1)*P.N_pol + pol, 'Parent', fig);
            v = x(:, pol);
            if numel(v) > maxPts
                v = v(round(linspace(1, numel(v), maxPts)));
            end
            plot(ax, real(v), imag(v), '.', 'MarkerSize', 3);
            axis(ax, 'equal'); grid(ax, 'on'); box(ax, 'on');
            title(ax, sprintf('%s  pol %d', lab, pol), 'Interpreter', 'none');
        end
    end
    drawnow;
    outFile = fullfile(fileparts(mfilename('fullpath')), ...
        'diag_pipeline_float_const.png');
    try
        exportgraphics(fig, outFile, 'Resolution', 150);
        fprintf('\nSaved constellations to %s\n', outFile);
    catch
        fprintf('\n(could not export constellation PNG)\n');
    end
end

function b = berAtOffset(crOut, txSymbols, off)
    refEnd = min(off + size(crOut,1), size(txSymbols,1));
    refSyms = txSymbols(off+1:refEnd, :);
    m = min(size(crOut,1), size(refSyms,1));
    if m < 1, b = NaN; return; end
    b = computeBERlocal(crOut(1:m,:), refSyms(1:m,:));
end

function BER = computeBERlocal(decoded, refSyms)
    nPol = size(refSyms,2); totErr=0; totBits=0;
    for p = 1:nPol
        refBits = modem.symbolsToBits(refSyms(:,p));
        best = Inf;
        for q = 1:size(decoded,2)
            for kk = 0:31
                rot = decoded(:,q).*exp(-1j*kk*pi/16);
                bb  = modem.symbolsToBits(modem.decideSymbols(rot));
                k   = min(numel(refBits), numel(bb));
                e   = sum(refBits(1:k) ~= bb(1:k))/k;
                if e < best, best = e; end
            end
        end
        totErr = totErr + best*numel(refBits); totBits = totBits + numel(refBits);
    end
    BER = totErr/totBits;
end

function s = tern(c,a,b), if c, s=a; else, s=b; end, end
