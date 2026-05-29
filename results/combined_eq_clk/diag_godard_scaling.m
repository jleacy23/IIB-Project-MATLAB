function diag_godard_scaling()
%DIAG_GODARD_SCALING  Measure the dynamic range of the Godard FD metric and
%   intermediate spectra (using the UNSCALED forward FFT the fxp path uses)
%   to see which fxp accumulator overflows.  Float maths; no MEX.

    here = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(fileparts(fileparts(here)), 'src')));

    % Params matching combined_eq_clk_fxp_sweep
    Rs=30.5; SpS=2; D=20; CWL=1550; DGDSpec=0.1; N_pmd=1; PMD_seed=12345;
    L_km=80; SFO_ppm=40; CFO_GHz=3.0; Rolloff=0.25; Span=10;
    NFFT=128; NCD=22; NOverlap=2*ceil((NCD-1)/2); Ns=4000; SNR_dB=20; N_pol=2;

    rng(13);
    symbols=(2*randi([0 1],Ns,N_pol)-1)+1j*(2*randi([0 1],Ns,N_pol)-1);
    tx=modem.rrcPulse(symbols,SpS,Rolloff,Span);
    rx=channel.add_chromatic_dispersion(tx,L_km,SpS,Rs,D,CWL);
    s=rng; rng(PMD_seed); rx=channel.add_pmd(rx,L_km,SpS,Rs,DGDSpec,N_pmd); rng(s);
    rx=channel.lo_freq_shift(rx,CFO_GHz*1000,Rs,SpS);
    rx=channel.apply_timing_error(rx,SFO_ppm,0,SpS);
    rx=channel.add_awgn(rx,SNR_dB);
    [rx,~]=eq_clk.coarse_cfo_fd(rx,NFFT);

    HCD=eq_clk.cd_fd_response(D,L_km,CWL,Rs,SpS,NFFT);
    HMF=eq_clk.rrc_fd_response(Rolloff,NFFT,SpS);
    Hstatic=ifftshift(HCD.*HMF);

    eta=SpS; shift=round((1-1/eta)*NFFT);
    kLo=round((1-Rolloff)/(2*eta)*NFFT)+1; kHi=round((1+Rolloff)/(2*eta)*NFFT);
    k_idx=[0:NFFT/2-1,-NFFT/2:-1].';

    stepLen=NFFT-NOverlap;
    nBlk=floor((size(rx,1)-NOverlap)/stepLen);

    maxRcorr=0; maxAbsS=0; maxAbsE=0; maxBandBin=0;
    for i=1:min(nBlk,40)
        wStart=(i-1)*stepLen+1;
        InB=rx(wStart:wStart+NFFT-1,:);
        S=0;
        for p=1:N_pol
            R=fft.fft_flp(InB(:,p),false,false);   % UNSCALED forward FFT (fxp convention)
            Rc=R.*Hstatic;                          % CD+MF (ramp~1 at start)
            maxRcorr=max(maxRcorr,max(abs(Rc)));
            maxBandBin=max(maxBandBin,max(abs(Rc(kLo:kHi))));
            S=S+sum(Rc(kLo:kHi).*conj(Rc(kLo+shift:kHi+shift)));
        end
        maxAbsS=max(maxAbsS,abs(S));
        maxAbsE=max(maxAbsE,abs(imag(S)));
    end

    fprintf('\n=== Godard FD metric dynamic range (unscaled forward FFT) ===\n');
    fprintf('  max |R_corr| (any bin)      = %.1f\n', maxRcorr);
    fprintf('  max |R_corr| (band kLo:kHi) = %.1f\n', maxBandBin);
    fprintf('  max |S| (metric)            = %.3e\n', maxAbsS);
    fprintf('  max |imag(S)| = |e|         = %.3e\n', maxAbsE);
    fprintf('\n  fxp type ranges (NIntBits=16):\n');
    fprintf('    T.Godard.metric  WL=32 FL=16 -> int range +/- %d\n', 2^15);
    fprintf('    T.Godard.lf (wide) WL=48 FL=40 -> int range +/- %d\n', 2^7);
    fprintf('\n  => metric overflow if max|S| > %d : %d\n', 2^15, maxAbsS>2^15);
    fprintf('  => loop-filter e cast overflow if |e| > %d : %d\n', 2^7, maxAbsE>2^7);
end
