function Y = apply_timing_error(X, SFO, tau0, SpS)
%APPLY_TIMING_ERROR  Resample with sample-frequency offset and constant phase.
%
%   Y = apply_timing_error(X, SFO, tau0, SpS)
%
%   Models a receiver whose sampling clock differs from the transmitter
%   clock by a constant rate offset (SFO, "frequency") and a fixed
%   timing phase (tau0).  Output sample n is taken from the input at
%   the fractional position
%
%       t(n) = (n + tau0 * SpS) / (1 + SFO * 1e-6)        (input samples)
%
%   so the receiver's nominal sampling instants drift linearly with
%   respect to the transmitter's grid.
%
%   Inputs
%     X    - input signal [Nsamp x Npol], complex
%     SFO  - sample-frequency offset [ppm].  Positive => RX clock faster
%            than TX clock (RX samples earlier in TX time per output sample).
%     tau0 - constant timing offset [symbol periods, fractional]
%     SpS  - input oversampling factor [samples per symbol]
%
%   Output
%     Y    - resampled signal [Nsamp x Npol].  Positions that fall
%            outside the input span are zero (spline extrapolation off).

    [N, P] = size(X);
    alpha  = SFO * 1e-6;
    n      = (0:N-1).';
    tIn    = n;
    tOut   = (n + tau0 * SpS) / (1 + alpha);

    Y = zeros(N, P, 'like', X);
    for p = 1:P
        Y(:, p) = interp1(tIn, X(:, p), tOut, 'spline', 0);
    end
end
