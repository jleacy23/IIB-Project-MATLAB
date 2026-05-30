function [normSignal, scale] = normalise(rxSignal, pct)
%NORMALISE Scale a received signal into the range [-1, 1] on each axis.
%   normSignal = normalise(rxSignal, pct) scales rxSignal so that its real
%   and imaginary parts fall within [-1, 1]. The scaling reference is the
%   pct-th percentile of the sample magnitudes (per channel), rather than the
%   peak, so that a small number of large outliers do not squash the bulk of
%   the data towards zero. Samples beyond the percentile are clipped to the
%   unit box [-1, 1] on each axis.
%
%   rxSignal   : [N x M] complex array (M independent channels/polarizations)
%   pct        : percentile of |rxSignal| used as the normalisation reference,
%                in (0, 100]. Optional, defaults to 99.
%
%   normSignal : [N x M] complex array scaled (and clipped) to the unit box
%   scale      : [1 x M] per-channel scale factors applied (normSignal is
%                rxSignal ./ scale before clipping)

    if nargin < 2 || isempty(pct)
        pct = 99;
    end

    if pct <= 0 || pct > 100
        error('modem:normalise:badPercentile', ...
            'pct must be in the interval (0, 100].');
    end

    % Per-channel percentile of the sample magnitude as the scale reference.
    % Sort each column ascending and pick the value at the percentile rank
    % (linear interpolation between order statistics, no toolbox dependency).
    mags = sort(abs(rxSignal), 1);
    N = size(mags, 1);
    rank = max(1, pct / 100 * N);   % 1-based fractional rank
    lo = floor(rank);
    hi = min(lo + 1, N);
    frac = rank - lo;
    scale = mags(lo, :) + frac .* (mags(hi, :) - mags(lo, :));

    % Guard against a degenerate (all-zero) channel.
    scale(scale == 0) = 1;

    normSignal = rxSignal ./ scale;

    % Clip outliers beyond the percentile back into the unit box per axis.
    realPart = min(max(real(normSignal), -1), 1);
    imagPart = min(max(imag(normSignal), -1), 1);

    normSignal = realPart + 1i * imagPart;
end
