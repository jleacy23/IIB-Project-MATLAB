function diag_perpol_fix()
%DIAG_PERPOL_FIX  Confirm the shared-pol carrier-recovery bug across trials
%   by comparing carrier_recovery.pilots_only run SHARED (both pols, as the
%   pipeline does) vs PER-POL (each polarisation independently).  Pure float.

    here = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(fileparts(fileparts(here)), 'src')));

    P.Rs=30.5; P.SpS=2; P.N_pol=2; P.D=20; P.CWL=1550; P.DGDSpec=0.1;
    P.N_pmd=1; P.PMD_seed=12345; P.L_km=80; P.SFO_ppm=40; P.CFO_GHz=3.0;
    P.LW_Hz=1e6; P.Rolloff=0.25; P.Span=10; P.N_sub_target=8;
    P.NFFT=128; P.NCD=22; P.NLanesGard=32;
    P.ki=1e-7; P.kp=1e-6;
    P.NTapsAEQ=1; P.MuAEQ=1e-3; P.N1AEQ=500; P.NOutAEQ=1000;
    P.SignOnly=true; P.SingleSpike=true; P.PLanesAEQ=32;
    P.MaxFreq=0.1; P.TrainingLen=11; P.BlockLen_CR=32;
    P.SUBFRAME_SYMS=3712; P.BLOCK_LEN=32; P.N_BLOCKS=116;
    SNR_dB = 40;

    fprintf('\n  trial |  BER shared-pol |  BER per-pol  | clustered?\n');
    fprintf('  ------+-----------------+---------------+-----------\n');
    for tr = 1:4
        [crShared, crPerPol, refSyms, blkErrShared] = runTrial(P, SNR_dB, tr);
        bS = computeBERlocal(crShared, refSyms);
        bP = computeBERlocal(crPerPol, refSyms);
        fracBadBlocks = mean(blkErrShared > 0.3);
        fprintf('   %2d   |   %.4e    |  %.4e   | %.1f%% blocks >30%% err\n', ...
            tr, bS, bP, bP, 100*fracBadBlocks);
    end
    fprintf(['\n  shared-pol = carrier_recovery.pilots_only as the pipeline calls it\n', ...
             '  per-pol    = pilots_only run once per polarisation (NPol=1 each)\n']);
end

function [crShared, crPerPol, refSyms, blkErrShared] = runTrial(P, SNR_dB, trialSeed)
    rng(1000*trialSeed + round(SNR_dB) + 13);
    nBits = P.N_sub_target * (3586*P.N_pol*2);
    [txSymbols, pilotsRef, training, ~] = modem.modulate(randi([0 1],nBits,1));
    NsymTx = size(txSymbols,1);

    txSig = modem.rrcPulse(txSymbols, P.SpS, P.Rolloff, P.Span);
    rxSig = channel.add_chromatic_dispersion(txSig, P.L_km, P.SpS, P.Rs, P.D, P.CWL);
    s=rng; rng(P.PMD_seed);
    rxSig = channel.add_pmd(rxSig, P.L_km, P.SpS, P.Rs, P.DGDSpec, P.N_pmd); rng(s);
    rxSig = channel.lo_freq_shift(rxSig, P.CFO_GHz*1000, P.Rs, P.SpS);
    rxSig = channel.apply_timing_error(rxSig, P.SFO_ppm, 0, P.SpS);
    rxSig = channel.add_awgn(rxSig, SNR_dB);
    rxSig = channel.add_phase_noise(rxSig, P.Rs*P.SpS, P.LW_Hz);

    nOv = 2*ceil((P.NCD-1)/2);
    adaptOpts = struct('NTaps',P.NTapsAEQ,'Mu',P.MuAEQ,'SingleSpike',P.SingleSpike, ...
        'N1',P.N1AEQ,'NOut',P.NOutAEQ,'SignOnly',P.SignOnly,'PLanes',P.PLanesAEQ, ...
        'Mode',0,'Pilots',[],'BlockLen',P.PLanesAEQ,'SubframeBlocks',0);
    eqOutFull = eq_clk.combined_cd_fd_gardner_adaptive(rxSig, P.SpS, P.NFFT, nOv, ...
        P.D, P.L_km, P.CWL, P.Rs, P.Rolloff, P.ki, P.kp, NsymTx, P.NLanesGard, ...
        adaptOpts, true, false);

    nDropEq = P.SUBFRAME_SYMS - P.NOutAEQ;
    eqOut = eqOutFull(nDropEq+1:end, :);
    nSubUsable = floor(size(eqOut,1)/P.SUBFRAME_SYMS);
    eqOut = eqOut(1:nSubUsable*P.SUBFRAME_SYMS, :);
    pIdx = pipeline_fxp_sweep.pilotTrainingIndices(nSubUsable, P);
    eqOut(pIdx,:) = 3*eqOut(pIdx,:);

    frOut = freq_recovery.differential_kay(eqOut, 3*training, P.Rs, true, 0);
    pilotsAll = 3*repmat(pilotsRef, nSubUsable, 1);

    % SHARED (pipeline): both pols -> one phase per block
    crShared = carrier_recovery.pilots_only(frOut, P.N_pol, P.BlockLen_CR, pilotsAll);

    % PER-POL: run the same function once per polarisation (NPol=1)
    cX = carrier_recovery.pilots_only(frOut(:,1), 1, P.BlockLen_CR, pilotsAll(:,1));
    cY = carrier_recovery.pilots_only(frOut(:,2), 1, P.BlockLen_CR, pilotsAll(:,2));
    crPerPol = [cX, cY];

    off = P.SUBFRAME_SYMS;
    refEnd = min(off+size(crShared,1), size(txSymbols,1));
    refSyms = txSymbols(off+1:refEnd, :);
    m = min(size(crShared,1), size(refSyms,1));
    crShared = crShared(1:m,:); crPerPol = crPerPol(1:m,:); refSyms = refSyms(1:m,:);

    % Per-block error rate (shared, X-pol) to see if errors cluster in blocks
    nB = floor(m/P.BlockLen_CR);
    blkErrShared = zeros(nB,1);
    [~,rotS] = computeBERlocal(crShared, refSyms);  %#ok
    for b=1:nB
        idx=(b-1)*P.BlockLen_CR+(1:P.BlockLen_CR);
        rb = modem.symbolsToBits(refSyms(idx,1));
        bb = modem.symbolsToBits(modem.decideSymbols(crShared(idx,1)*exp(-1j*rotS)));
        k=min(numel(rb),numel(bb)); blkErrShared(b)=sum(rb(1:k)~=bb(1:k))/k;
    end
end

function [BER,bestRot] = computeBERlocal(decoded, refSyms)
    nPol=size(refSyms,2); totErr=0; totBits=0; bestRot=0;
    for p=1:nPol
        refBits=modem.symbolsToBits(refSyms(:,p)); best=Inf;
        for q=1:size(decoded,2)
            for kk=0:31
                rot=decoded(:,q).*exp(-1j*kk*pi/16);
                bb=modem.symbolsToBits(modem.decideSymbols(rot));
                k=min(numel(refBits),numel(bb));
                e=sum(refBits(1:k)~=bb(1:k))/k;
                if e<best, best=e; if p==1,bestRot=kk*pi/16; end, end
            end
        end
        totErr=totErr+best*numel(refBits); totBits=totBits+numel(refBits);
    end
    BER=totErr/totBits;
end
