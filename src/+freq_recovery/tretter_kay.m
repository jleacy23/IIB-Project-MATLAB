function [y, frequency_offset] = tretter_kay(x, training, Rs)
%TRETTER_KAY  Tretter/Kay frequency offset estimator and corrector.
%
%   [y, frequency_offset] = tretter_kay(x, training, Rs)
%
%   Estimates the carrier frequency offset from the weighted phase
%   differences of the training-sequence correlation.  The offset is
%   estimated independently for each polarisation and then averaged
%   before correcting the full subframe.
%
%   Inputs
%     x        - input subframe  [Nsym x NPol]
%     training - training symbols [L x NPol]
%     Rs       - symbol rate [GBd]
%
%   Outputs
%     y                - frequency-corrected subframe [Nsym x NPol]
%     frequency_offset - estimated frequency offset [kHz]

    [L, N_pol] = size(training);
    Nsym       = size(x, 1);

    %% Weights w(k) = 6k(L-k) / (L(L^2-1)) — one column per polarisation
    k       = (1 : L-1).';                           % [(L-1) x 1]
    w_col   = 6 .* k .* (L - k) / (L * (L^2 - 1));  % [(L-1) x 1]
    weights = repmat(w_col, 1, N_pol);                % [(L-1) x NPol]

    %% De-rotate received training with known sequence
    x_tr = x(1:L, :);                                % [L x NPol]
    z    = x_tr .* conj(training);                   % [L x NPol]

    %% Weighted phase-difference sum — one offset per polarisation
    dz        = z(2:end, :) .* conj(z(1:end-1, :)); % [(L-1) x NPol]
    f_coeff   = (Rs * 1e9) / (2 * pi);
    f_per_pol = f_coeff * sum(weights .* angle(dz), 1); % [1 x NPol] Hz

    %% Average offset across polarisations then correct full subframe
    frequency_offset = mean(f_per_pol);              % scalar [Hz]

    n = (0 : Nsym-1).';                              % [Nsym x 1]
    y = x .* exp(-1j * 2*pi * frequency_offset / (Rs*1e9) .* n);

    frequency_offset = frequency_offset / 1e3;       % Hz -> kHz
end










