function diag_fxp_stages()
%DIAG_FXP_STAGES  Localise the residual fxp pipeline bug.  Uses the
%   (validated) FLOAT equaliser output, then runs FR and CR in BOTH float
%   and fixed-point (interpreted, no MEX) at FL=16 to see which fxp stage
%   diverges.  Mirrors pipeline_fxp_sweep.runOneTrial exactly otherwise.

    here = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(fileparts(fileparts(here)), 'src')));

    P.Rs=30.5; P.SpS=2; P.N_pol=2; P.D=20; P.CWL=1550; P.DGDSpec=0.1;
    P.N_pmd=1; P.PMD_seed=12345; P.L_km=80; P.SFO_ppm=40; P.CFO_GHz=3.0;
    P.LW_Hz=1e6; P.Rolloff=0.25; P.Span=10; P.N_sub_target=8;
    P.NFFT=128; P.NCD=22; P.NLanesGard=32; P.ki=1e-7; P.kp=1e-6;
    P.NTapsAEQ=1; P.MuAEQ=1e-3; P.N1AEQ=500; P.NOutAEQ=1000;
    P.SignOnly=true; P.SingleSpike=true; P.PLanesAEQ=32;
    P.MaxFreq=0.1; P.TrainingLen=11; P.CordicIts=16; P.BlockLen_CR=32;
    P.SUBFRAME_SYMS=3712; P.BLOCK_LEN=32; P.N_BLOCKS=116;
    P.IntBits=16; FL=16;
    SNR_dB=40; trialSeed=1;

    %% Channel + FLOAT equaliser ------------------------------------
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
    adaptOpts = struct('NTaps',P.NTapsAEQ,'Mu',P.MuAEQ,'SingleSpike',P.SingleSpike, ...
        'N1',P.N1AEQ,'NOut',P.NOutAEQ,'SignOnly',P.SignOnly,'PLanes',P.PLanesAEQ, ...
        'Mode',0,'Pilots',[],'BlockLen',P.PLanesAEQ,'SubframeBlocks',0);
    eqOutFull = eq_clk.combined_cd_fd_gardner_adaptive(rxSig, P.SpS, P.NFFT, nOv, ...
        P.D, P.L_km, P.CWL, P.Rs, P.Rolloff, P.ki, P.kp, NsymTx, P.NLanesGard, adaptOpts, true, false);

    nDropEq = P.SUBFRAME_SYMS - P.NOutAEQ;
    eqOut = eqOutFull(nDropEq+1:end, :);
    nSubUsable = floor(size(eqOut,1)/P.SUBFRAME_SYMS);
    eqOut = eqOut(1:nSubUsable*P.SUBFRAME_SYMS, :);
    pIdx = pipeline_fxp_sweep.pilotTrainingIndices(nSubUsable, P);
    eqOut(pIdx,:) = 3*eqOut(pIdx,:);
    pilotsAll = 3*repmat(pilotsRef, nSubUsable, 1);

    T_fr = freq_recovery.fxp_types(struct('WL',P.IntBits+FL,'FL',FL));
    T_cr = carrier_recovery.fxp_types(struct('WL',P.IntBits+FL,'FL',FL));

    fprintf('\n================ FXP STAGE LOCALISATION (FL=%d) ===========\n', FL);
    fprintf('Float equaliser output is the common starting point.\n');

    %% ---- FR: float vs fxp -----------------------------------------
    [frF, fHzF] = freq_recovery.differential_kay(eqOut, 3*training, P.Rs, true, 0);
    eqOut_fi    = cast(eqOut, 'like', T_fr.x);
    training_fi = cast(3*training, 'like', T_fr.x);
    [frX_fi, fHzX] = freq_recovery.differential_kay_fxp(eqOut_fi, training_fi, ...
        P.Rs, P.CordicIts, T_fr, true, 0, P.MaxFreq);
    frX = double(frX_fi);
    fprintf('\n[FR] float f=%.4f GHz | fxp f=%.4f GHz | diff=%.4f GHz\n', ...
        fHzF/1e9, fHzX/1e9, (fHzF-fHzX)/1e9);

    %% ---- CR: float vs fxp (on float-FR output) --------------------
    crFF = carrier_recovery.pilots_only(frF, P.N_pol, P.BlockLen_CR, pilotsAll);
    crFX_fi = carrier_recovery.pilots_only_fxp(cast(frF,'like',T_cr.x), P.N_pol, ...
        P.BlockLen_CR, cast(pilotsAll,'like',T_cr.x), P.CordicIts, T_cr);
    crFX = double(crFX_fi);
    fprintf('[CR] on FLOAT-FR output:  BER(floatCR)=%.4e   BER(fxpCR)=%.4e\n', ...
        berAt(crFF, txSymbols, P.SUBFRAME_SYMS), berAt(crFX, txSymbols, P.SUBFRAME_SYMS));

    %% ---- CR on fxp-FR output --------------------------------------
    crXF = carrier_recovery.pilots_only(frX, P.N_pol, P.BlockLen_CR, pilotsAll);
    crXX_fi = carrier_recovery.pilots_only_fxp(cast(frX,'like',T_cr.x), P.N_pol, ...
        P.BlockLen_CR, cast(pilotsAll,'like',T_cr.x), P.CordicIts, T_cr);
    crXX = double(crXX_fi);
    fprintf('[CR] on FXP-FR   output:  BER(floatCR)=%.4e   BER(fxpCR)=%.4e\n', ...
        berAt(crXF, txSymbols, P.SUBFRAME_SYMS), berAt(crXX, txSymbols, P.SUBFRAME_SYMS));

    %% ---- Diagnose the FR phase ramp range -------------------------
    % differential_kay_fxp stores theta/max_freq in T.theta and wraps at
    % +/- pi/max_freq.  Check that range fits T.theta and that the applied
    % phase ramp actually matches the float one.
    fprintf('\n[FR detail] max_freq=%.3f -> theta wraps at +/-pi/max_freq=%.2f\n', ...
        P.MaxFreq, pi/P.MaxFreq);
    info = numerictype(T_fr.theta);
    fprintf('     T.theta: signed=%d WL=%d FL=%d -> max representable=%.2f\n', ...
        info.Signed, info.WordLength, info.FractionLength, ...
        double(fi(0,info.Signed,info.WordLength,info.FractionLength).range.max ...
                * 0 + 2^(info.WordLength-info.FractionLength-info.Signed)));
    % Compare the actual applied de-rotation float vs fxp (residual phase):
    resF = frF .* conj(eqOut) ./ max(abs(eqOut),eps);  %#ok  (informational)
    dphi = angle(frX(:,1)) - angle(frF(:,1));
    dphi = mod(dphi+pi,2*pi)-pi;
    fprintf('     angle(frX)-angle(frF) X-pol: std=%.4f rad (0 => same ramp)\n', std(dphi));

    fprintf('\n=========================================================\n');
end

function b = berAt(crOut, txSymbols, off)
    refEnd = min(off+size(crOut,1), size(txSymbols,1));
    refSyms = txSymbols(off+1:refEnd,:);
    m = min(size(crOut,1), size(refSyms,1));
    if m<1, b=NaN; return; end
    crOut=crOut(1:m,:); refSyms=refSyms(1:m,:);
    nPol=size(refSyms,2); totErr=0; totBits=0;
    for p=1:nPol
        rb=modem.symbolsToBits(refSyms(:,p)); best=Inf;
        for q=1:size(crOut,2)
            for kk=0:31
                bb=modem.symbolsToBits(modem.decideSymbols(crOut(:,q).*exp(-1j*kk*pi/16)));
                k=min(numel(rb),numel(bb)); e=sum(rb(1:k)~=bb(1:k))/k;
                if e<best, best=e; end
            end
        end
        totErr=totErr+best*numel(rb); totBits=totBits+numel(rb);
    end
    b=totErr/totBits;
end
