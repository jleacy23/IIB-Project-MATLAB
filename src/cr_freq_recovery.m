function [z, Delta_f] = cr_freq_recovery(y, Rs)
    % unit conversion
    rs = Rs * 1e9;  % GBd -> Bd

    f = (-1/2+1/length(y):1/length(y):1/2)*rs;
    Ts =1/rs;
    SignalSpectrum= fftshift(abs(fft(y(:,1).^4)));
    % Restrict peak search to positive frequencies
    posIdx = f >= 0;
    specPos = SignalSpectrum;
    specPos(~posIdx) = 0;
    delta_f = (1/4)*f(specPos==max(specPos));

    k = repmat((0:length(y)-1).',1,size(y,2));

    z = y.*exp(-1i*2*pi*delta_f*Ts*k);
    Delta_f = delta_f * 1e-6;
end

