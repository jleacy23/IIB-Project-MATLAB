function diag_eq_localize()
%DIAG_EQ_LOCALIZE  Find which sub-stage of combined_cd_fd_gardner_adaptive_fxp
%   breaks at FL=16.  Uses the mixed-precision equaliser MEX built by
%   build_diag_mex (eqHI/eqS16/eqC16/eqA16 + eq16), each with two sub-stages
%   near-lossless (WL=48/FL=32) and one at WL=32/FL=16, isolating Static /
%   Clk / AdaptEq.  FR+CR use the FL=16 MEX (already shown clean).
%
%   Requires build_diag_mex to have been run locally first.

    here = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(fileparts(fileparts(here)), 'src')));

    P.Rs=30.5; P.SpS=2; P.N_pol=2; P.D=20; P.CWL=1550; P.DGDSpec=0.1;
    P.N_pmd=1; P.PMD_seed=12345; P.L_km=80; P.SFO_ppm=40; P.CFO_GHz=3.0;
    P.LW_Hz=1e6; P.Rolloff=0.25; P.Span=10; P.N_sub_target=8;
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
    lo  = struct('WL',P.IntBits+FL,'FL',FL);   % the FL=16 stage under test
    hi  = struct('WL',48,'FL',32);             % near-lossless (16 int bits)

    %  label, Static, Clk, AdaptEq, MEX base name (must match build_diag_mex)
    eq16name = pipeline_fxp_sweep.mexEqName(FL);
    cases = {
        'all hi (baseline)', hi, hi, hi, 'combined_cd_fd_gardner_adaptive_fxp_eqHI_mex'
        'Static = FL16',     lo, hi, hi, 'combined_cd_fd_gardner_adaptive_fxp_eqS16_mex'
        'Clk    = FL16',     hi, lo, hi, 'combined_cd_fd_gardner_adaptive_fxp_eqC16_mex'
        'AdaptEq = FL16',    hi, hi, lo, 'combined_cd_fd_gardner_adaptive_fxp_eqA16_mex'
        'all FL16',          lo, lo, lo, eq16name };

    fprintf('\n=== EQ sub-stage localisation (MEX, N_sub=%d, FL=%d) ===\n', P.N_sub_target, FL);
    for ci = 1:size(cases,1)
        T = eq_clk.combined_cd_fd_gardner_adaptive_fxp_types( ...
            struct('Static',cases{ci,2}, 'Clk',cases{ci,3}, 'AdaptEq',cases{ci,4}));
        rx_fi = cast(rxSig, 'like', T.Static.x);
        PilotsEmpty = cast(complex(zeros(0,2)), 'like', T.AdaptEq.y);
        aOpts = struct('NTaps',double(P.NTapsAEQ),'Mu',double(P.MuAEQ), ...
            'SingleSpike',logical(P.SingleSpike),'N1',double(P.N1AEQ),'NOut',double(P.NOutAEQ), ...
            'SignOnly',logical(P.SignOnly),'UpdateStep',double(1),'PLanes',double(P.PLanesAEQ), ...
            'Mode',double(0),'Pilots',PilotsEmpty,'BlockLen',double(P.PLanesAEQ),'SubframeBlocks',double(0));
        eqMex = str2func(['eq_clk.' cases{ci,5}]);
        t=tic;
        eqo = double(eqMex(rx_fi, double(P.SpS), ...
            double(P.NFFT), double(nOv), double(P.D), double(P.L_km), double(P.CWL), ...
            double(P.Rs), double(P.Rolloff), double(P.ki), double(P.kp), double(NsymTx), ...
            double(P.NLanesGard), aOpts, true, false, T));
        [ber, berPol] = pipeFinish(eqo, txSymbols, pilotsRef, training, P, FL);
        xpol = crossPol(eqo);   % ~1 => CMA singularity (both outputs same source)
        fprintf('  %-20s : BER=%.4e  perPol=[%.3f %.3f]  xpolCorr=%.3f  (%.1fs)\n', ...
            cases{ci,1}, ber, berPol(1), berPol(2), xpol, toc(t));
    end
    fprintf('================================================================\n');
    fprintf('perPol=[X Y] best BER per reference pol; one ~0 + one ~0.5 => singularity.\n');
    fprintf('xpolCorr near 1 => both equaliser outputs locked to the same polarisation.\n');
end

function c = crossPol(y)
    % Phase/lag-robust normalised cross-correlation between the two output
    % polarisations.  ~1 indicates both columns carry the same source
    % polarisation (CMA singularity); ~0 indicates clean separation.
    y1 = y(:,1) - mean(y(:,1));  y2 = y(:,2) - mean(y(:,2));
    c  = abs(sum(y1 .* conj(y2))) / (sqrt(sum(abs(y1).^2) * sum(abs(y2).^2)) + eps);
end

function [ber, berPol] = pipeFinish(eqOut0, txSymbols, pilotsRef, training, P, FL)
    nDropEq = P.SUBFRAME_SYMS - P.NOutAEQ;
    eqOut = eqOut0(nDropEq+1:end, :);
    nSub = floor(size(eqOut,1)/P.SUBFRAME_SYMS);
    if nSub < 1, ber = NaN; return; end
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
    nPol=size(refSyms,2); berPol=zeros(1,nPol);
    for p=1:nPol
        rb=modem.symbolsToBits(refSyms(:,p)); best=Inf;
        for q=1:size(crX,2)
            for kk=0:31
                bb=modem.symbolsToBits(modem.decideSymbols(crX(:,q).*exp(-1j*kk*pi/16)));
                k=min(numel(rb),numel(bb)); e=sum(rb(1:k)~=bb(1:k))/k;
                if e<best, best=e; end
            end
        end
        berPol(p)=best;
    end
    ber=mean(berPol);
end
