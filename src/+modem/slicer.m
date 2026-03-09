function s_dec = slicer(s_rot) %#codegen
%SLICER  QPSK hard decision — returns nearest ±1 ±1j constellation point.
%
%   s_dec = slicer(s_rot)
%
%   For QPSK the decision is simply the sign of the real and imaginary
%   parts.  This is far cheaper than the generic M-QAM slicer.
%
%   Input
%     s_rot - complex sample (scalar, vector or array)
%
%   Output
%     s_dec - complex nearest QPSK point (±1 ±1j), same size as s_rot

    d_re = sign(real(s_rot));
    d_im = sign(imag(s_rot));

    % Handle exact zero (map to +1 by convention)
    d_re(d_re == 0) = 1;
    d_im(d_im == 0) = 1;

    s_dec = complex(d_re, d_im);
end