function diag_fft_precision()
%DIAG_FFT_PRECISION  Measure the static-stage (forward FFT -> CD/MF mask ->
%   inverse FFT) accuracy vs fractional length, to see if the 1/2-per-stage
%   forward-FFT truncation is what causes the FL<=8 cliff.  Interpreted; no MEX.

    here = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(fileparts(fileparts(here)), 'src')));

    Rs=30.5; SpS=2; D=20; CWL=1550; L_km=80; Rolloff=0.25; NFFT=128; NIntBits=16;

    HCD=eq_clk.cd_fd_response(D,L_km,CWL,Rs,SpS,NFFT);
    HMF=eq_clk.rrc_fd_response(Rolloff,NFFT,SpS);
    Hstatic=ifftshift(HCD.*HMF);

    rng(1);
    x = (randn(NFFT,1)+1j*randn(NFFT,1))/sqrt(2);   % ~unit-power block

    % Float reference: new-convention round trip (forward 1/N, inverse x1)
    Xf = fft.fft_flp(x,false,false);
    yf = fft.fft_flp(Xf.*Hstatic,true,false);

    fprintf('\n=== Static round-trip NRMSE vs FL  (CD/MF, NFFT=%d) ===\n', NFFT);
    fprintf('  FL | NEW (wide acc, 1/N output) | OLD (unscaled fwd, 1/N end)\n');
    fprintf('  ---+----------------------------+---------------------------\n');
    for FL = [16 12 10 8 6 4]
        T = cd_eq.equalize_fxp_types(struct('WL',NIntBits+FL,'FL',FL));
        Tf.x=T.x; Tf.tw=T.tw; Tf.acc=T.acc;
        Hfi = cast(Hstatic,'like',T.hcd);
        x_fi = cast(x,'like',T.x);

        % NEW convention (current fft_fxp)
        Xn = fft.fft_fxp(x_fi,false,false,Tf);
        Yn = complex(zeros(NFFT,1,'like',T.acc));
        for k=1:NFFT, Yn(k)=Xn(k)*Hfi(k); end
        yn = fft.fft_fxp(Yn,true,false,Tf);
        nrmseNew = norm(double(yn)-yf)/norm(yf);

        % OLD convention (local fi FFT: unscaled fwd, 1/N on inverse)
        Xo = fftFi(x_fi,false,Tf);
        Yo = complex(zeros(NFFT,1,'like',T.acc));
        for k=1:NFFT, Yo(k)=Xo(k)*Hfi(k); end
        yo = fftFi(Yo,true,Tf);
        nrmseOld = norm(double(yo)-yf)/norm(yf);

        fprintf('  %2d | %26.4e | %25.4e\n', FL, nrmseNew, nrmseOld);
    end
    fprintf('\n  NEW keeps the spectrum ~O(1) (no overflow, needed for Godard);\n');
    fprintf('  OLD keeps the spectrum ~O(N) (better low-FL SNR, but Godard overflows).\n');
end

function y = builtinScale(x, inverse, N)
    if inverse, y = ifft(x) * N; else, y = fft(x) / N * N; end %#ok
    if ~inverse, y = fft(x); end   % old fwd == unscaled fft
    if inverse,  y = ifft(x) * N; end % old inv == ifft*N (since the H-mult input was unscaled fft)
end

function X = fftFi(x, inverse, T)
%FFTFI  OLD-convention fi radix-2: unscaled forward, 1/N on inverse.
    N = size(x,1); numStages = round(log2(double(N)));
    X = complex(zeros(N,1,'like',T.acc));
    for i = 0:N-1
        rev = bitrevLocal(i,numStages);
        X(rev+1) = cast(x(i+1),'like',T.acc);
    end
    for s = 1:numStages
        halfLen=2^(s-1); fullLen=2^s; numGroups=N/fullLen;
        for k = 0:halfLen-1
            theta = -2*pi*double(k)/double(fullLen);
            if inverse, theta = -theta; end
            W = complex(cast(cos(theta),'like',T.tw), cast(sin(theta),'like',T.tw));
            for g = 0:numGroups-1
                it = g*fullLen+k+1; ib = it+halfLen;
                u = X(it); t = W*X(ib);
                X(it) = u + t;      % no per-stage scaling
                X(ib) = u - t;
            end
        end
    end
    if inverse
        M = cast(log2(double(N)),'int32');
        for i = 1:N, X(i) = bitshift(X(i), -M); end   % 1/N once at the end
    end
end

function rev = bitrevLocal(idx, nbits)
    rev=0; val=idx;
    for b=1:nbits, rev=rev*2+mod(val,2); val=floor(val/2); end
end
