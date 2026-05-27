function Y = add_iq_skew(X, Skew, Rs, SpS)
%ADD_IQ_SKEW  Introduce a timing skew between the I and Q components.
%
%   Y = add_iq_skew(X, Skew, Rs, SpS)
%
%   Models a receiver in which the in-phase and quadrature sampling paths
%   are not perfectly aligned in time.  The Q component is delayed by
%   +Skew/2 and the I component by -Skew/2, so the mean sampling instant
%   of each polarisation is preserved.
%
%   Inputs
%     X    - input signal [Nsamp x Npol], complex
%     Skew - I/Q timing skew [ps].  Positive => Q sampled later than I.
%     Rs   - symbol rate [GBd]
%     SpS  - oversampling factor [samples per symbol]
%
%   Output
%     Y    - skewed signal [Nsamp x Npol].  Samples that fall outside the
%            input span are set to zero.
%
%   The skew is applied independently to each polarisation via spline
%   interpolation in the time domain, matching the style used by
%   channel.apply_timing_error.
%
%   See also channel.apply_timing_error.
    [N, P]   = size(X);
    Ts       = 1 / (Rs * 1e9 * SpS);    % sample period [s]
    dSamp    = (Skew * 1e-12) / Ts;     % skew in input samples
    n        = (0:N-1).';
    tI       = n - dSamp/2;             % I sampled dSamp/2 earlier
    tQ       = n + dSamp/2;             % Q sampled dSamp/2 later

    Y = zeros(N, P, 'like', X);
    for p = 1:P
        I = interp1(n, real(X(:, p)), tI, 'spline', 0);
        Q = interp1(n, imag(X(:, p)), tQ, 'spline', 0);
        Y(:, p) = I + 1i * Q;
    end
end
