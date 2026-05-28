function diag_pipeline_float(LW_override, CFO_override)
%DIAG_PIPELINE_FLOAT  Floating-point clone of pipeline_fxp_sweep.runOneTrial
%   with heavy per-stage instrumentation.  No MEX / no fi: every stage uses
%   the float reference so we can localise the BER failure to STRUCTURE
%   (alignment / pilot-ID / CFO / loop tuning) versus FIXED-POINT.
%
%   Run:  addpath(genpath('src')); diag_pipeline_float

    here     = fileparts(mfilename('fullpath'));
    repoRoot = fileparts(fileparts(here));
    addpath(genpath(fullfile(repoRoot, 'src')));

    %% ---- Parameters (mirror pipeline_fxp_sweep) --------------------
    P.Rs=30.5; P.SpS=2; P.N_pol=2; P.D=20; P.CWL=1550; P.DGDSpec=0.1;
    P.N_pmd=1; P.PMD_seed=12345; P.L_km=80; P.SFO_ppm=40; P.CFO_GHz=3.0;
    P.LW_Hz=1e6; P.Rolloff=0.25; P.Span=10; P.N_sub_target=8;
    P.NFFT=128; P.NCD=22; P.NLanesGard=32;
    P.ki=1e-7; P.kp=1e-6;                  % po2-off gains
    P.NTapsAEQ=1; P.MuAEQ=1e-3; P.N1AEQ=500; P.NOutAEQ=1000;
    P.SignOnly=true; P.SingleSpike=true; P.PLanesAEQ=32;
    P.MaxFreq=0.1; P.TrainingLen=11; P.BlockLen_CR=32;
    P.SUBFRAME_SYMS=3712; P.BLOCK_LEN=32; P.N_BLOCKS=116;

    if nargin>=1 && ~isempty(LW_override),  P.LW_Hz  = LW_override;  end
    if nargin>=2 && ~isempty(CFO_override), P.CFO_GHz = CFO_override; end

    SNR_dB = 40; trialSeed = 1;
    fprintf('\n================ FLOAT PIPELINE DIAGNOSTIC ===============\n');
    fprintf('SNR=%g dB  CFO=%g GHz  trial=%d\n', SNR_dB, P.CFO_GHz, trialSeed);

    %% ---- 1. Channel ------------------------------------------------
    rng(1000*trialSeed + round(SNR_dB) + 13);
    CPON_BITS_PER_SF = 3586 * P.N_pol * 2;
    nBits = P.N_sub_target * CPON_BITS_PER_SF;
    bits  = randi([0 1], nBits, 1);
    [txSymbols, pilotsRef, training, nSubTx] = modem.modulate(bits);
    NsymTx = size(txSymbols,1);
    fprintf('\n[1] Tx: %d symbols, %d subframes. pilotsRef=[%dx%d] training=[%dx%d]\n', ...
        NsymTx, nSubTx, size(pilotsRef,1), size(pilotsRef,2), size(training,1), size(training,2));

    % Sanity: verify pilot/training positions in the *Tx* stream match the
    % index formula used downstream (pilotTrainingIndices).
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

    %% ---- 2. Equaliser (FLOAT) -------------------------------------
    nOv = 2*ceil((P.NCD-1)/2);
    adaptOpts = struct('NTaps',P.NTapsAEQ,'Mu',P.MuAEQ,'SingleSpike',P.SingleSpike, ...
        'N1',P.N1AEQ,'NOut',P.NOutAEQ,'SignOnly',P.SignOnly,'PLanes',P.PLanesAEQ, ...
        'Mode',0,'Pilots',[],'BlockLen',P.PLanesAEQ,'SubframeBlocks',0);
    [eqOutFull, cfoBins] = eq_clk.combined_cd_fd_gardner_adaptive(rxSig, P.SpS, ...
        P.NFFT, nOv, P.D, P.L_km, P.CWL, P.Rs, P.Rolloff, P.ki, P.kp, ...
        NsymTx, P.NLanesGard, adaptOpts, true, false);
    fprintf('\n[2] Equaliser out: %d symbols. cfoBinsApplied=%.4f (bin=%.4f GHz -> %.4f GHz removed)\n', ...
        size(eqOutFull,1), cfoBins, P.SpS*P.Rs/P.NFFT, cfoBins*P.SpS*P.Rs/P.NFFT);

    %% ---- 3. Drop one subframe; align ------------------------------
    nDropEq = P.SUBFRAME_SYMS - P.NOutAEQ;
    eqOut = eqOutFull(nDropEq+1:end, :);
    nSubUsable = floor(size(eqOut,1)/P.SUBFRAME_SYMS);
    eqOut = eqOut(1:nSubUsable*P.SUBFRAME_SYMS, :);
    fprintf('[3] nDropEq=%d -> %d subframes usable (%d symbols)\n', ...
        nDropEq, nSubUsable, size(eqOut,1));

    %% ---- 3b. ALIGNMENT PROBE (suspect 1) --------------------------
    % Scan a lag window: for each candidate lag, compute the contrast
    % between |eqOut| at the assumed pilot positions and elsewhere.
    % A correct lag => pilots have sqrt(2) magnitude and a distinctive
    % pattern; wrong lag => indistinguishable from data.
    fprintf('\n[3b] ALIGNMENT SCAN (|pilot| contrast vs lag):\n');
    bestLag = 0; bestContrast = -inf;
    for lag = -40:40
        c = pilotContrast(eqOutFull, nDropEq, lag, nSubUsable, P);
        if ~isnan(c) && c > bestContrast, bestContrast = c; bestLag = lag; end
    end
    for lag = [-2 -1 0 1 2 bestLag]
        c = pilotContrast(eqOutFull, nDropEq, lag, nSubUsable, P);
        fprintf('     lag=%+3d : pilot-pattern correlation = %+.4f %s\n', ...
            lag, c, tern(lag==bestLag,'<== best',''));
    end
    fprintf('     => best lag = %+d (0 means drop logic is correct)\n', bestLag);

    pIdx = pipeline_fxp_sweep.pilotTrainingIndices(nSubUsable, P);
    nonP = setdiff(1:size(eqOut,1), pIdx);
    fprintf('     |eqOut(pilot/train)| = %.3f   |eqOut(data)| = %.3f   (ratio %.2f)\n', ...
        mean(abs(eqOut(pIdx,1))), mean(abs(eqOut(nonP,1))), ...
        mean(abs(eqOut(pIdx,1)))/mean(abs(eqOut(nonP,1))));

    %% ---- 4. Rescale pilots/training x3 ----------------------------
    eqOut(pIdx,:) = 3*eqOut(pIdx,:);

    %% ---- 5. Frequency recovery (FLOAT differential_kay) -----------
    [frOut, frHz] = freq_recovery.differential_kay(eqOut, 3*training, P.Rs, true, 0);
    fprintf('\n[5] FR estimate = %.4f GHz residual. (CFO removed by eq + FR = %.4f GHz; target %.4f)\n', ...
        frHz/1e9, cfoBins*P.SpS*P.Rs/P.NFFT + frHz/1e9, P.CFO_GHz);
    % Residual frequency drift after FR, inferred from pilot phase slope:
    reportPilotPhaseSlope(frOut, pilotsRef, nSubUsable, P, 'after FR');

    %% ---- 6. Carrier recovery (FLOAT pilots_only) ------------------
    pilotsAll = 3*repmat(pilotsRef, nSubUsable, 1);
    [crOut, ThetaPU] = carrier_recovery.pilots_only(frOut, P.N_pol, P.BlockLen_CR, pilotsAll);
    dth = diff(ThetaPU(1:P.BlockLen_CR:end, 1));
    dth = mod(dth+pi, 2*pi)-pi;
    fprintf('\n[6] CR per-block phase step: mean=%.4f rad  std=%.4f rad  max|.|=%.4f\n', ...
        mean(dth), std(dth), max(abs(dth)));
    fprintf('     (max|step| near pi/2=%.3f => ambiguity-jump risk)\n', pi/2);

    % Per-POL pilot phase (the CR collapses both pols into ONE estimate).
    NB = size(pilotsAll,1);
    thX = nan(NB,1); thY = nan(NB,1); thXY = nan(NB,1);
    for b = 1:NB
        bs = (b-1)*P.BlockLen_CR + 1;
        if bs <= size(frOut,1)
            thX(b)  = angle(conj(pilotsAll(b,1))*frOut(bs,1));
            thY(b)  = angle(conj(pilotsAll(b,2))*frOut(bs,2));
            thXY(b) = angle(sum(conj(pilotsAll(b,:)).*frOut(bs,:)));
        end
    end
    dXY = mod(thX-thY+pi,2*pi)-pi;
    fprintf('     per-pol phase: std(thetaX)=%.3f std(thetaY)=%.3f  std(thetaX-thetaY)=%.3f rad\n', ...
        std(thX,'omitnan'), std(thY,'omitnan'), std(dXY,'omitnan'));
    fprintf('     => if std(thetaX-thetaY) is large, the single shared CR estimate is invalid\n');
    % Per-pol BER using each pol''s OWN phase estimate (ideal pilot CR per pol):
    crIdealX = frOut(:,1).*exp(-1j*holdPhase(thX,P.BlockLen_CR,size(frOut,1)));
    crIdealY = frOut(:,2).*exp(-1j*holdPhase(thY,P.BlockLen_CR,size(frOut,1)));
    crIdeal  = [crIdealX, crIdealY];
    bIdeal = berAtOffset(crIdeal, txSymbols, P.SUBFRAME_SYMS);
    fprintf('     BER with PER-POL pilot phase (ideal per-pol CR) = %.4e\n', bIdeal);

    %% ---- 7. BER (assumed alignment) + lag scan --------------------
    txOffset = P.SUBFRAME_SYMS;
    fprintf('\n[7] BER vs small symbol lag (txOffset=%d):\n', txOffset);
    for lag = [-2 -1 0 1 2]
        b = berAtOffset(crOut, txSymbols, txOffset+lag);
        fprintf('     lag=%+d : BER=%.4e %s\n', lag, b, tern(lag==0,'<== pipeline assumption',''));
    end

    fprintf('\n================ END DIAGNOSTIC ==========================\n');
end

% ---------------------------------------------------------------------
function checkTxPilotPlacement(tx, pilotsRef, training, P)
    base = 0;  % subframe 1
    okTrain = isequal(tx(base+(1:P.TrainingLen),:), training);
    pilotPos = base + (2:P.N_BLOCKS-0)*0;  %#ok
    okP = true;
    for b = 2:P.N_BLOCKS
        if ~isequal(tx(base+(b-1)*P.BLOCK_LEN+1,:), pilotsRef(b,:)), okP=false; break; end
    end
    okP1 = isequal(tx(base+1,:), pilotsRef(1,:));
    fprintf('     Tx placement check: training@1..11=%d, pilot@blockstarts=%d, TS1==pilot(1)=%d\n', ...
        okTrain, okP, okP1);
end

function c = pilotContrast(eqOutFull, nDropEq, lag, nSubUsable, P)
    % Correlate |eqOut| at assumed pilot positions (shifted by lag) against
    % the expected on/off pilot mask, normalised.  High => pilots land right.
    start = nDropEq + 1 + lag;
    if start < 1, c = NaN; return; end
    avail = size(eqOutFull,1) - start + 1;
    nSub  = min(nSubUsable, floor(avail / P.SUBFRAME_SYMS));
    if nSub < 1, c = NaN; return; end
    seg = eqOutFull(start : start + nSub*P.SUBFRAME_SYMS - 1, :);
    pIdx = pipeline_fxp_sweep.pilotTrainingIndices(nSub, P);
    mask = false(size(seg,1),1); mask(pIdx) = true;
    a = abs(seg(:,1)) + abs(seg(:,2));
    % point-biserial correlation between magnitude and pilot mask
    c = (mean(a(mask)) - mean(a(~mask))) / (std(a) + eps);
end

function reportPilotPhaseSlope(x, pilotsRef, nSubUsable, P, tag)
    pilotsAll = repmat(pilotsRef, nSubUsable, 1);
    NB = size(pilotsAll,1);
    th = zeros(NB,1);
    for b = 1:NB
        bs = (b-1)*P.BlockLen_CR + 1;
        if bs <= size(x,1)
            th(b) = angle(sum(conj(pilotsAll(b,:)) .* x(bs,:)));
        end
    end
    dth = mod(diff(th)+pi,2*pi)-pi;
    fprintf('     [%s] pilot-phase block-to-block: mean=%.4f std=%.4f rad/blk (drift=>residual freq)\n', ...
        tag, mean(dth), std(dth));
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

function th = holdPhase(thBlk, BlockLen, Nsym)
    th = zeros(Nsym,1);
    for b = 1:numel(thBlk)
        i0 = (b-1)*BlockLen+1; i1 = min(b*BlockLen, Nsym);
        if i0 <= Nsym, th(i0:i1) = thBlk(b); end
    end
end

function s = tern(c,a,b), if c, s=a; else, s=b; end, end
