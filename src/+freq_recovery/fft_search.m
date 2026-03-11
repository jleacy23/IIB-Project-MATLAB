function [y, frequency_offset] = fft_search(x, training, Rs, K, data_aided, D)
%FFT_SEARCH  FFT-based frequency offset estimator and corrector.
%
%   [y, frequency_offset] = fft_search(x, training, Rs)
%   [y, frequency_offset] = fft_search(x, training, Rs, K)
%   [y, frequency_offset] = fft_search(x, training, Rs, K, data_aided, D)
%
%   Builds an observation sequence, zero-pads by factor K, takes the FFT,
%   finds the peak bin with parabolic interpolation, averages across
%   polarisations, then applies the correction to the full subframe.
%
%   When data_aided = true (default), the observation sequence is formed
%   by de-rotating the L training symbols against the known sequence.
%   When data_aided = false, the D data symbols immediately following the
%   training block are raised to the 4th power to remove modulation.
%
%   Inputs
%     x          - input subframe   [Nsym x NPol]
%     training   - known training symbols [L x NPol]
%     Rs         - symbol rate [GBd]
%     K          - zero-padding factor (integer >= 1); default 8
%     data_aided - true = training-aided (default), false = blind 4th-power
%     D          - number of data symbols in blind mode (required when
%                  data_aided = false)
%
%   Outputs
%     y                - frequency-corrected subframe [Nsym x NPol]
%     frequency_offset - estimated frequency offset [Hz]

    if nargin < 4 || isempty(K), K = 8; end
    if nargin < 5, data_aided = true; end

    [L, N_pol] = size(training);
    Nsym       = size(x, 1);

    if data_aided
        %% Training-aided: de-rotate with known sequence
        z    = x(1:L, :) .* conj(training);   % [L x NPol]
        No   = L;
    else
        %% Blind: 4th-power of D data symbols after training block
        x_data = x(L+1 : L+D, :);             % [D x NPol]
        z      = x_data .^ 4;                 % [D x NPol]
        No     = D;
    end

    Nfft = K * No;

    %% Zero-pad and FFT — one spectrum per polarisation
    Z  = fft(z, Nfft, 1);                           % [Nfft x NPol]

    %% Coarse search: peak bin per polarisation
    Zs       = fftshift(Z, 1);
    Zmag_s   = abs(Zs);
    [~, idx] = max(Zmag_s, [], 1);                  % [1 x NPol]  1-indexed

    %% Fine interpolation using the three bins around the peak
    idx_m = mod(idx - 2, Nfft) + 1;
    idx_p = mod(idx,     Nfft) + 1;

    delta = zeros(1, N_pol);
    for p = 1:N_pol
        Xm = Zs(idx_m(p), p);
        Xk = Zs(idx(p),   p);
        Xp = Zs(idx_p(p), p);
        delta(p) = -real((Xp - Xm) / (2*Xk - Xm - Xp));
    end

    idx_fine  = idx + delta;
    f_per_pol = (idx_fine - 1 - Nfft/2) ./ Nfft * (Rs * 1e9);  % [1 x NPol] Hz

    %% In blind mode the 4th power maps f0 -> 4*f0; undo the factor
    if ~data_aided
        f_per_pol = f_per_pol / 4;
    end

    %% Average offset across polarisations then correct full subframe
    frequency_offset = mean(f_per_pol);

    n = (0 : Nsym-1).';
    y = x .* exp(-1j * 2*pi * frequency_offset / (Rs*1e9) .* n);
end
