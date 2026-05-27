function [rxOut, cfoBins] = coarse_cfo_fd(rxIn, NFFT)
%COARSE_CFO_FD  One-shot coarse CFO estimation and time-domain correction.
%
%   [rxOut, cfoBins] = coarse_cfo_fd(rxIn, NFFT)
%
%   Estimate the CFO once from the magnitude-squared centroid of the
%   averaged FFT spectrum (Welch periodogram over all NFFT-length
%   non-overlapping blocks of the input) and apply the corresponding
%   continuous-valued time-domain phasor to the whole signal.
%
%   Averaging is necessary because the centroid of a single NFFT-block
%   FFT of a random-data QPSK signal has a per-block standard deviation
%   of order a bin (~0.5 GHz at NFFT = 128, Rs = 30.5 GBd, SpS = 2),
%   which dominates the CFO range of interest.  Averaging across the
%   whole input reduces this variance by ~sqrt(nBlocks) and yields a
%   sub-bin estimate that the time-domain phasor applies with no
%   rounding loss.
%
%   The matched filter is NOT applied before the centroid: it is
%   centered at DC and would asymmetrically attenuate the high-frequency
%   edge of a CFO-shifted signal, biasing the centroid toward DC.  CD
%   is all-pass in magnitude so it does not bias |R|^2.
%
%   Inputs
%     rxIn  - input signal [samples x NPol]
%     NFFT  - block size used for the periodogram
%
%   Outputs
%     rxOut   - corrected signal [samples x NPol]
%     cfoBins - estimated CFO in FFT bins (signed).  Bin spacing in Hz
%               is SpS*Rs/NFFT.

    nSamples = size(rxIn, 1);
    nBlk     = floor(nSamples / NFFT);

    k_idx = [0:NFFT/2-1, -NFFT/2:-1].';
    mag2  = zeros(NFFT, 1);
    for i = 1:nBlk
        s = (i-1) * NFFT + 1;
        e = s + NFFT - 1;
        R = fft(rxIn(s:e, :), [], 1);
        mag2 = mag2 + sum(abs(R).^2, 2);
    end

    den = sum(mag2);
    if den > 0
        cfoBins = sum(k_idx .* mag2) / den;
    else
        cfoBins = 0;
    end

    n      = (0 : nSamples - 1).';
    phasor = exp(-1j * 2*pi * cfoBins * n / NFFT);
    rxOut  = rxIn .* phasor;
end
