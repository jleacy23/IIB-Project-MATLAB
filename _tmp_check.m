addpath(genpath('src'));
fprintf('--- checkcode equalize.m ---\n');
checkcode('src/+adaptive_eq/equalize.m');
fprintf('--- checkcode equalize_fxp.m ---\n');
checkcode('src/+adaptive_eq/equalize_fxp.m');

rng(1);
SpS=2; NTaps=7; BlockLen=32;
Nsym = 32*20;
sym = (2*randi([0 1],Nsym,2)-1)+1j*(2*randi([0 1],Nsym,2)-1);
Pilots = zeros(Nsym/BlockLen,2);
for b=1:Nsym/BlockLen
  p = (sign(randn(1,2))*3) + 1j*(sign(randn(1,2))*3);
  sym((b-1)*BlockLen+1,:) = p; Pilots(b,:) = p;
end
x = repelem(sym,SpS,1);
x = x + 0.05*(randn(size(x))+1j*randn(size(x)));

yC0 = adaptive_eq.equalize(x,SpS,NTaps,1e-3,true,50,0,false,1);
yC  = adaptive_eq.equalize(x,SpS,NTaps,1e-3,true,50,0,false,1,0,Pilots,BlockLen);
yP  = adaptive_eq.equalize(x,SpS,NTaps,1e-3,true,50,0,false,1,1,Pilots,BlockLen);
fprintf('sizes  legacy=%dx%d  cmaSkip=%dx%d  pilot=%dx%d\n', size(yC0),size(yC),size(yP));
fprintf('finite legacy=%d cmaSkip=%d pilot=%d\n', all(isfinite(yC0(:))), all(isfinite(yC(:))), all(isfinite(yP(:))));
fprintf('CMA skip differs from legacy: %d (expect 1)\n', any(abs(yC(:)-yC0(:))>0));

yLegacy1 = adaptive_eq.equalize(x,SpS,NTaps,1e-3,true,50,0,false,1);
yNew1    = adaptive_eq.equalize(x,SpS,NTaps,1e-3,true,50,0,false,1,0,[],[]);
fprintf('backward-compat identical: %d\n', isequal(yLegacy1,yNew1));
fprintf('DONE\n');
