function s_dec = slicer(s_rot, M) %#codegen
%SLICER  Nearest-neighbour decision for square M-QAM, unit average power.
%
%   s_dec = slicer(s_rot, M)
%
%   Inputs
%     s_rot - complex double rotated sample (scalar)
%     M     - QAM order (4, 16, 64, 256, ...)  Must be a perfect square.
%
%   Output
%     s_dec - complex double nearest constellation point (unit average power)
%
    % Normalisation factor and constellation half-extent
    k      = sqrt((2 * (M - 1)) / 3);   % scale: normalised -> unnormalised
    extent = sqrt(M) - 1;               % maximum unnormalised level

    % Real axis decision
    x_re  = real(s_rot) * k;
    d_re  = 2 * round((x_re - 1) / 2) + 1;
    d_re  = min(extent, max(-extent, d_re));

    % Imaginary axis decision
    x_im  = imag(s_rot) * k;
    d_im  = 2 * round((x_im - 1) / 2) + 1;
    d_im  = min(extent, max(-extent, d_im));

    % Scale back to unit average power
    s_dec = complex(d_re / k, d_im / k);
end