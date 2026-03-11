function [y, frequency_offset] = differential(x, training, Rs, data_aided, D)
%DIFFERENTIAL  Simple differential-phase frequency offset estimator.
%
%   [y, frequency_offset] = differential(x, training, Rs)
%   [y, frequency_offset] = differential(x, training, Rs, data_aided, D)
%
%   Forms the same observation sequence as the Tretter-Kay estimator and
%   computes successive phase differences.  Rather than applying the
%   optimal Kay weights, each difference is treated equally and the sum
%   is divided by the observation length No:
%
%       f0_hat = Rs / (2*pi) * (1/No) * sum_k angle( z(k+1)*conj(z(k)) )
%
%   The estimate is computed independently for each polarisation and then
%   averaged before correcting the full subframe.
%
%   When data_aided = true (default), the observation sequence is formed
%   by de-rotating the L training symbols against the known sequence.
%   When data_aided = false, the D data symbols immediately following the
%   training block are raised to the 4th power to remove modulation.
%
%   Inputs
%     x          - input subframe  [Nsym x NPol]
%     training   - known training symbols [L x NPol]
%     Rs         - symbol rate [GBd]
%     data_aided - true = training-aided (default), false = blind 4th-power
%     D          - number of data symbols to use in blind mode (required
%                  when data_aided = false)
%
%   Outputs
%     y                - frequency-corrected subframe [Nsym x NPol]
%     frequency_offset - estimated frequency offset [Hz]

    if nargin < 4, data_aided = true; end

    [L, N_pol] = size(training);
    Nsym       = size(x, 1);

    if data_aided
        %% Training-aided: de-rotate with known sequence
        z  = x(1:L, :) .* conj(training);   % [L x NPol]
        No = L;
    else
        %% Blind: 4th-power of the D data symbols after training block
        x_data = x(L+1 : L+D, :);           % [D x NPol]
        z      = x_data .^ 4;               % [D x NPol]
        No     = D;
    end

    %% Successive phase differences
    dz = z(2:end, :) .* conj(z(1:end-1, :));   % [No-1 x NPol]

    %% Average phase difference, scaled by observation length No (not No-1)
    f_coeff   = (Rs * 1e9) / (2 * pi);
    f_per_pol = f_coeff * sum(angle(dz), 1) / No;   % [1 x NPol] Hz

    %% In blind mode the 4th power maps f0 -> 4*f0; undo the factor
    if ~data_aided
        f_per_pol = f_per_pol / 4;
    end

    %% Average offset across polarisations then correct full subframe
    frequency_offset = mean(f_per_pol);

    n = (0 : Nsym-1).';
    y = x .* exp(-1j * 2*pi * frequency_offset / (Rs*1e9) .* n);
end