function diag_fxp_eq()
%DIAG_FXP_EQ  Run the FIXED-POINT equaliser (interpreted, no MEX) at FL=16
%   and feed it through the validated fxp FR + CR, to confirm the residual
%   pipeline bug lives in the equaliser.  Smaller N_sub_target for speed.

    here = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(fileparts(fileparts(here)), 'src')));

    P.Rs=30.5; P.SpS=2; P.N_pol=2; P.D=20; P.CWL=1550; P.DGDSpec=0.1;
    P.N_pmd=1; P.PMD_seed=12345; P.L_km=80; P.SFO_ppm=40; P.CFO_GHz=3.0;
    P.LW_Hz=1e6; P.Rolloff=0.25; P.Span=10; P.N_sub_target=8;   % full, MEX is fast
    P.NFFT=128; P.NCD=22; P.NLanesGard=32; P.ki=1e-7; P.kp=1e-6;
    P.NTapsAEQ=1; P.MuAEQ=1e-3; P.N1AEQ=500; P.NOutAEQ=1000;
    P.SignOnly=true; P.SingleSpike=true; P.PLanesAEQ=32;
    P.MaxFreq=0.1; P.TrainingLen=11; P.CordicIts=16; P.BlockLen_CR=32;
    P.SUBFRAME_SYMS=3712; P.BLOCK_LEN=32; P.N_BLOCKS=116; P.IntBits=16;
    FL=16; SNR_dB=40; trialSeed=1;

    rng(1000*trialSeed + round(SNR_dB) + 13);
    [txSymbols, pilotsRef, training, ~] = modem.modulate(randi([0 1], P.N_sub_target*3586*P.N_pol*2,1));
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
    adaptOptsF = struct('NTaps',P.NTapsAEQ,'Mu',P.MuAEQ,'SingleSpike',P.SingleSpike, ...
        'N1',P.N1AEQ,'NOut',P.NOutAEQ,'SignOnly',P.SignOnly,'PLanes',P.PLanesAEQ, ...
        'Mode',0,'Pilots',[],'BlockLen',P.PLanesAEQ,'SubframeBlocks',0);

    %% Float eq (reference) -----------------------------------------
    eqF = eq_clk.combined_cd_fd_gardner_adaptive(rxSig, P.SpS, P.NFFT, nOv, ...
        P.D, P.L_km, P.CWL, P.Rs, P.Rolloff, P.ki, P.kp, NsymTx, P.NLanesGard, adaptOptsF, true, false);

    %% Fxp eq (interpreted) -----------------------------------------
    T_eq = eq_clk.combined_cd_fd_gardner_adaptive_fxp_types( ...
        struct('Static',struct('WL',P.IntBits+FL,'FL',FL), ...
               'Clk',   struct('WL',P.IntBits+FL,'FL',FL), ...
               'AdaptEq',struct('WL',P.IntBits+FL,'FL',FL)));
    rxSig_fi = cast(rxSig,'like',T_eq.Static.x);
    PilotsEmpty = cast(complex(zeros(0,2)),'like',T_eq.AdaptEq.y);
    adaptOptsX = struct('NTaps',double(P.NTapsAEQ),'Mu',double(P.MuAEQ), ...
        'SingleSpike',logical(P.SingleSpike),'N1',double(P.N1AEQ),'NOut',double(P.NOutAEQ), ...
        'SignOnly',logical(P.SignOnly),'UpdateStep',double(1),'PLanes',double(P.PLanesAEQ), ...
        'Mode',double(0),'Pilots',PilotsEmpty,'BlockLen',double(P.PLanesAEQ),'SubframeBlocks',double(0));
    fprintf('Running fxp equaliser MEX...\n'); t0=tic;
    eqX_fi = eq_clk.combined_cd_fd_gardner_adaptive_fxp_eq16_mex(rxSig_fi, double(P.SpS), ...
        double(P.NFFT), double(nOv), double(P.D), double(P.L_km), double(P.CWL), ...
        double(P.Rs), double(P.Rolloff), double(P.ki), double(P.kp), double(NsymTx), ...
        double(P.NLanesGard), adaptOptsX, true, false, T_eq);
    eqX = double(eqX_fi);
    fprintf('  done (%.0fs). float eq len=%d, fxp eq len=%d\n', toc(t0), size(eqF,1), size(eqX,1));

    berF = pipeFinish(eqF, txSymbols, pilotsRef, training, P, FL);
    berX = pipeFinish(eqX, txSymbols, pilotsRef, training, P, FL);
    fprintf('\n  BER (FLOAT eq -> fxp FR/CR) = %.4e\n', berF);
    fprintf('  BER (FXP   eq -> fxp FR/CR) = %.4e\n', berX);

    % If fxp-eq BER is bad, inspect the constellation magnitude/EVM.
    fprintf('\n  |eqF| mean=%.3f std=%.3f | |eqX| mean=%.3f std=%.3f\n', ...
        mean(abs(eqF(:))), std(abs(eqF(:))), mean(abs(eqX(:))), std(abs(eqX(:))));
end

function ber = pipeFinish(eqOut0, txSymbols, pilotsRef, training, P, FL)
    nDropEq = P.SUBFRAME_SYMS - P.NOutAEQ;
    eqOut = eqOut0(nDropEq+1:end, :);
    nSub = floor(size(eqOut,1)/P.SUBFRAME_SYMS);
    eqOut = eqOut(1:nSub*P.SUBFRAME_SYMS, :);
    pIdx = pipeline_fxp_sweep.pilotTrainingIndices(nSub, P);
    eqOut(pIdx,:) = 3*eqOut(pIdx,:);
    pilotsAll = 3*repmat(pilotsRef, nSub, 1);
    T_fr = freq_recovery.fxp_types(struct('WL',P.IntBits+FL,'FL',FL));
    T_cr = carrier_recovery.fxp_types(struct('WL',P.IntBits+FL,'FL',FL));
    frX = double(freq_recovery.differential_kay_fxp_fc16_mex(cast(eqOut,'like',T_fr.x), ...
        cast(3*training,'like',T_fr.x), P.Rs, P.CordicIts, T_fr, true, 0, P.MaxFreq));
    crX = double(carrier_recovery.pilots_only_fxp_fc16_mex(cast(frX,'like',T_cr.x), P.N_pol, ...
        P.BlockLen_CR, cast(pilotsAll,'like',T_cr.x), P.CordicIts, T_cr));
    off=P.SUBFRAME_SYMS; refEnd=min(off+size(crX,1),size(txSymbols,1));
    refSyms=txSymbols(off+1:refEnd,:); m=min(size(crX,1),size(refSyms,1));
    crX=crX(1:m,:); refSyms=refSyms(1:m,:);
    nPol=size(refSyms,2); te=0; tb=0;
    for p=1:nPol
        rb=modem.symbolsToBits(refSyms(:,p)); best=Inf;
        for q=1:size(crX,2)
            for kk=0:31
                bb=modem.symbolsToBits(modem.decideSymbols(crX(:,q).*exp(-1j*kk*pi/16)));
                k=min(numel(rb),numel(bb)); e=sum(rb(1:k)~=bb(1:k))/k;
                if e<best, best=e; end
            end
        end
        te=te+best*numel(rb); tb=tb+numel(rb);
    end
    ber=te/tb;
end
