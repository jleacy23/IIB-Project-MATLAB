function [y, frequency_offset] = fft_search(x, training, Rs, K)
%FFT_SEARCH  FFT-based frequency offset estimator and corrector.
%
%   [y, frequency_offset] = fft_search(x, training, Rs)
%   [y, frequency_offset] = fft_search(x, training, Rs, K)
%
%   De-rotates the received training symbols against the known sequence,
%   zero-pads by factor K to interpolate the frequency axis, takes the
%   FFT, finds the peak bin, averages the offset estimate across
%   polarisations, then applies the correction to the full subframe.
%
%   Inputs
%     x        - input subframe   [Nsym x NPol]
%     training - known training symbols [L x NPol]
%     Rs       - symbol rate [GBd]
%     K        - zero-padding (pruning) factor (integer >= 1); default 8
%
%   Outputs
%     y                - frequency-corrected subframe [Nsym x NPol]
%     frequency_offset - estimated frequency offset [kHz]

    if nargin < 4
        K = 8;
    end

    [L, N_pol] = size(training);
    Nsym       = size(x, 1);
    Nfft       = K * L;

    %% De-rotate received training against known sequence
    % z[n] ≈ |training[n]|^2 * exp(j*2*pi*(f0/Rs)*n)  (noiseless)
    x_tr = x(1:L, :);                               % [L x NPol]
    z    = x_tr .* conj(training);                   % [L x NPol]

    %% Zero-pad and FFT — one spectrum per polarisation
    Z    = fft(z, Nfft, 1);                          % [Nfft x NPol]

    %% Coarse search: peak bin per polarisation
    % fftshift centres DC so bin 1 = -Rs/2, bin (Nfft/2+1) = 0
    Zs        = fftshift(Z, 1);                      % [Nfft x NPol] complex
    Zmag_s    = abs(Zs);                             % [Nfft x NPol]
    [~, idx]  = max(Zmag_s, [], 1);                  % [1 x NPol]  1-indexed

    %% Fine interpolation using the three bins around the peak
    % Use 1-based mod to wrap at the spectrum edges
    idx_m = mod(idx - 2, Nfft) + 1;   % k-1 (with wraparound)
    idx_p = mod(idx,     Nfft) + 1;   % k+1 (with wraparound)

    delta = zeros(1, N_pol);
    for p = 1:N_pol
        Xm = Zs(idx_m(p), p);   % X_{k-1}
        Xk = Zs(idx(p),   p);   % X_k
        Xp = Zs(idx_p(p), p);   % X_{k+1}

        % TODO: fill in interpolation formula using Xm, Xk, Xp
        delta(p) = -real((Xp - Xm) / (2*Xk - Xm - Xp));
    end

    % Refined fractional bin index (still 1-indexed, after fftshift)
    idx_fine = idx + delta;                          % [1 x NPol]

    % Bin k (1-indexed, after fftshift) -> frequency: (k-1-Nfft/2)/Nfft * Rs
    f_per_pol = (idx_fine - 1 - Nfft/2) ./ Nfft * (Rs * 1e9);  % [1 x NPol] Hz

    %% Average offset across polarisations then correct full subframe
    frequency_offset = mean(f_per_pol);              % scalar [Hz]

    n = (0 : Nsym-1).';                              % [Nsym x 1]
    y = x .* exp(-1j * 2*pi * frequency_offset / (Rs*1e9) .* n);

    frequency_offset = frequency_offset / 1e3;       % Hz -> kHz
end
