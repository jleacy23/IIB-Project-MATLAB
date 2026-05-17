function smoke_mex()
    root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
    addpath(fullfile(root,'src')); addpath(fullfile(root,'build'));

    P.N_pol=2; P.Nsym=8000; P.SpS=2; P.Rs=32; P.NTaps=3; P.Mu=1e-3;
    P.SingleSpike=true; P.N1=1000; P.NOut=2000; P.L_km=80; P.DGDSpec=0.5;
    P.N_pmd=5; P.VV_NTaps=5; P.VV_BlockLen=100; P.VV_PilotThresh=1e6;
    P.VV_FxpConfig=struct('WL',64,'FL',32); P.CordicIts=16;
    P.IntBits=16; P.DpWL=48; P.DpFL=32; P.HiWL=32; P.HiFL=24;
    P.WWL=40; P.WFL=32;

    gradFL = 16;
    F = fimath('RoundingMethod','Floor','OverflowAction','Wrap', ...
        'ProductMode','SpecifyPrecision','ProductWordLength',P.DpWL, ...
        'ProductFractionLength',P.DpFL,'SumMode','SpecifyPrecision', ...
        'SumWordLength',P.DpWL,'SumFractionLength',P.DpFL);
    gWL=P.IntBits+gradFL;
    T_eq.x=fi([],1,P.HiWL,P.HiFL,F); T_eq.w=fi([],1,P.WWL,P.WFL,F);
    T_eq.y=fi([],1,P.HiWL,P.HiFL,F); T_eq.acc=fi([],1,P.HiWL,P.HiFL,F);
    T_eq.err=fi([],1,P.HiWL,P.HiFL,F); T_eq.grad=fi([],1,gWL,gradFL,F);
    T_eq.R_CMA=fi([],1,P.HiWL,P.HiFL,F);

    srcDir=fullfile(root,'src');
    clear mex
    mf=fullfile(srcDir,'+adaptive_eq',['equalize_fxp_mex.' mexext]);
    if isfile(mf), delete(mf); end
    x=fi(complex(0,0),numerictype(T_eq.x),fimath(T_eq.x));
    InType=coder.typeof(x,[Inf 2],[true false]);
    args={InType,double(P.SpS),double(P.NTaps),double(P.Mu), ...
        logical(P.SingleSpike),double(P.N1),double(P.NOut), ...
        logical(false),double(1),T_eq};
    cfg=coder.config('mex'); cfg.GenerateReport=false;
    fprintf('building equalize MEX...\n');
    codegen('-config',cfg,'adaptive_eq.equalize_fxp','-args',args, ...
        '-o',fullfile(srcDir,'+adaptive_eq','equalize_fxp_mex'));

    B.N_pol=P.N_pol; B.VV_NTaps=P.VV_NTaps; B.PilotLen=1;
    B.BlockLen=P.VV_BlockLen; B.StepSize=P.VV_BlockLen;
    B.FxpConfig_VV=P.VV_FxpConfig; B.PilotThreshold=P.VV_PilotThresh;
    B.CordicIts=P.CordicIts;
    fprintf('building VV MEX...\n');
    build_carrier_recovery_viterbiViterbi_fxp_mex(B,cfg);

    T_vv=carrier_recovery.fxp_types(P.VV_FxpConfig);

    rng(1001);
    re=2*randi([0 1],P.Nsym,P.N_pol)-1; im=2*randi([0 1],P.Nsym,P.N_pol)-1;
    syms=(re+1j*im)/sqrt(2); txSig=repelem(syms,P.SpS,1);
    rng(777); rxSig=channel.add_pmd(txSig,P.L_km,P.SpS,P.Rs,P.DGDSpec,P.N_pmd);
    rng(42);  rxSig=channel.add_awgn(rxSig,22);

    rxFi=cast(rxSig,'like',T_eq.x);
    tic;
    eqSig=adaptive_eq.equalize_fxp_mex(rxFi,double(P.SpS),double(P.NTaps), ...
        double(P.Mu),true,double(P.N1),double(P.NOut),false,double(4),T_eq);
    tEq=toc;
    eqSig=double(eqSig);

    Lflt=2*P.VV_NTaps+1;
    wVV=carrier_recovery.genVVFilter(1e3,P.Rs,22,1,P.N_pol,P.VV_NTaps);
    wVV=wVV(:); if numel(wVV)~=Lflt, wVV=ones(Lflt,1)/Lflt; end
    NB=ceil(size(eqSig,1)/P.VV_BlockLen);
    xfi=cast(eqSig,'like',T_vv.x); wfi=cast(wVV,'like',T_vv.w);
    pfi=cast(zeros(NB,P.N_pol),'like',T_vv.x);
    tic;
    [v,~]=carrier_recovery.viterbiViterbi_fxp_mex(xfi,double(P.N_pol), ...
        double(P.VV_NTaps),wfi,pfi,double(P.VV_BlockLen), ...
        double(P.VV_BlockLen),double(P.VV_PilotThresh), ...
        double(P.CordicIts),T_vv);
    tVV=toc;
    v=double(v);

    L=size(v,1); ref=syms(P.NOut+1:P.NOut+L,:);
    g=P.VV_NTaps; v=v(g+1:end-g,:); ref=ref(g+1:end-g,:);
    totErr=0; totBits=0;
    for p=1:P.N_pol
        rb=modem.symbolsToBits(ref(:,p)); best=Inf;
        for q=1:P.N_pol
            for k=0:3
                rot=v(:,q).*exp(-1j*k*pi/2);
                bb=modem.symbolsToBits(modem.decideSymbols(rot));
                n=min(numel(rb),numel(bb));
                e=sum(rb(1:n)~=bb(1:n))/n; if e<best, best=e; end
            end
        end
        totErr=totErr+best*numel(rb); totBits=totBits+numel(rb);
    end
    fprintf('MEX OK | eq=%.3fs vv=%.3fs | BER@22dB=%.3e\n', ...
        tEq,tVV,totErr/totBits);
end
