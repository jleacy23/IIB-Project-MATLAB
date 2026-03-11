function [y, frequency_offset] = differential_kay(x, training, Rs, data_aided, D)
%DIFFERENTIAL_KAY  Two-stage differential-phase + Tretter-Kay estimator.
%
%   [y, frequency_offset] = differential_kay(x, training, Rs)
%   [y, frequency_offset] = differential_kay(x, training, Rs, data_aided, D)
%
%   Stage 1 – Coarse:  runs the simple differential-phase estimator to
%   obtain a coarse frequency estimate f_coarse and applies that correction
%   to the full subframe.
%
%   Stage 2 – Fine:  runs the Tretter-Kay (Kay 1989) weighted
%   phase-difference estimator on the coarsely corrected subframe to
%   estimate the residual offset f_fine.
%
%   The combined correction is applied to the original input:
%
%       frequency_offset = f_coarse + f_fine
%
%   Both stages operate in the same mode (training-aided or blind 4th-power)
%   as selected by data_aided and D.
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

    Nsym  = size(x, 1);
    Rs_Hz = Rs * 1e9;

    %% Stage 1: coarse estimate via differential-phase
    if data_aided
        [x_coarse, f_coarse] = freq_recovery.differential(x, training, Rs);
    else
        [x_coarse, f_coarse] = freq_recovery.differential(x, training, Rs, false, D);
    end

    %% Stage 2: fine estimate via Tretter-Kay on the coarsely corrected signal
    if data_aided
        [~, f_fine] = freq_recovery.tretter_kay(x_coarse, training, Rs);
    else
        [~, f_fine] = freq_recovery.tretter_kay(x_coarse, training, Rs, false, D);
    end

    %% Combine and apply total correction to the original input
    frequency_offset = f_coarse + f_fine;

    n = (0 : Nsym-1).';
    y = x .* exp(-1j * 2*pi * frequency_offset / Rs_Hz .* n);
end
