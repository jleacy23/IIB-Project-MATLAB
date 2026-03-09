function [y, frequency_offset] = fitz(x, training, Rs, N)
%FITZ  Fitz (1994) autocorrelation-based frequency offset estimator.
%
%   [y, frequency_offset] = fitz(x, training, Rs, N)
%
%   De-rotates the received training symbols against the known sequence to
%   give z(k), computes the lag-N autocorrelation R(N) = sum z(k)*conj(z(k-N)),
%   and estimates the frequency offset as
%
%       f0 = Rs / (2*pi*N) * angle( R(N) )      [Hz]
%
%   The estimate is computed independently for each polarisation and then
%   averaged before correcting the full subframe.
%
%   Inputs
%     x        - input subframe  [Nsym x NPol]
%     training - known training symbols [L x NPol]
%     Rs       - symbol rate [GBd]
%     N        - autocorrelation lag (integer, 1 <= N < L); default floor(L/2)
%
%   Outputs
%     y                - frequency-corrected subframe [Nsym x NPol]
%     frequency_offset - estimated frequency offset [kHz]

    [L, N_pol] = size(training);
    Nsym       = size(x, 1);

    if nargin < 4
        N = floor(L / 2);
    end

    assert(N >= 1 && N < L, 'Lag N must satisfy 1 <= N < L.');

    %% De-rotate received training with known sequence
    x_tr = x(1:L, :);                                % [L x NPol]
    z    = x_tr .* conj(training);                   % [L x NPol]

    %% Lag-N autocorrelation — one value per polarisation
    % R(N) = sum_{k=N+1}^{L} z(k) * conj(z(k-N))
    R = sum(z(N+1:end, :) .* conj(z(1:end-N, :)), 1);  % [1 x NPol]

    %% Frequency estimate per polarisation
    f_per_pol = (Rs * 1e9) / (2 * pi * N) * angle(R);   % [1 x NPol] Hz

    %% Average over polarisations then correct full subframe
    frequency_offset = mean(f_per_pol);               % scalar [Hz]

    n = (0 : Nsym-1).';                               % [Nsym x 1]
    y = x .* exp(-1j * 2*pi * frequency_offset / (Rs*1e9) .* n);

    frequency_offset = frequency_offset / 1e3;        % Hz -> kHz
end
