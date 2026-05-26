function HCD = cd_fd_response(D, L, CLambda, Rs, SpSIn, NFFT)
%CD_FD_RESPONSE  Chromatic-dispersion frequency response (fftshifted grid).
%
%   HCD = cd_fd_response(D, L, CLambda, Rs, SpSIn, NFFT)
%
%   Returns the (NFFT x 1) frequency response of a CD-compensating filter
%   on the centred grid n = -NFFT/2 .. NFFT/2-1 (i.e. the same ordering
%   used by cd_eq.equalize).  Apply with fftshifted spectra, or use
%   ifftshift(HCD) when multiplying the natural-order output of fft().

    c = 299792458;

    D_si       = D * 1e-6;
    L_si       = L * 1e3;
    CLambda_si = CLambda * 1e-9;
    Rs_si      = Rs * 1e9;

    n  = (-NFFT/2:NFFT/2-1)';
    fN = SpSIn * Rs_si / 2;

    HCD = exp(-1i*pi*CLambda_si^2*D_si*L_si/c * (n*2*fN/NFFT).^2);
end
