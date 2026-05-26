function HMF = rrc_fd_response(beta, NFFT, SpS)
%RRC_FD_RESPONSE  Root-raised-cosine matched-filter frequency response.
%
%   HMF = rrc_fd_response(beta, NFFT, SpS)
%
%   Analytical RRC magnitude response on the centred fftshifted grid
%   n = -NFFT/2 .. NFFT/2-1.  Normalised frequency f*T = n*SpS/NFFT.
%   Real, even, zero-group-delay — multiplying R = fftshift(fft(x)) by
%   HMF realises a zero-delay matched filter.

    n   = (-NFFT/2:NFFT/2-1)';
    fn  = n * SpS / NFFT;
    af  = abs(fn);

    HMF = zeros(NFFT, 1);
    lo  = (1 - beta) / 2;
    hi  = (1 + beta) / 2;

    HMF(af <= lo) = 1;
    mask = (af > lo) & (af <= hi);
    HMF(mask) = sqrt(0.5 * (1 + cos(pi/beta * (af(mask) - lo))));
end
