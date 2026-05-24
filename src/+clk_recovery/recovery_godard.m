function Out = recovery_godard(In, NSymb, N, beta)
%recovery_godard  Feedforward clock recovery using Modified Godard.
%
%   Out = recovery_godard(In, NSymb, N, beta)
%
%   Both the timing estimate and the timing correction are performed in
%   the frequency domain on the same FFT block of N samples.
%
%   The Modified Godard feedforward estimate (Josten et al., Appl. Sci.
%   2017, eq. 5 with the arg/(2π) feedforward normalisation) is
%
%       τ̂_b / T = arg( Σ R(k)·R*(k+(1-1/η)N) ) / (2π)
%
%   summed over the excess-bandwidth bins.  The same block FFT is then
%   corrected with the linear phase ramp from eq. 1,
%
%       Y_corr(k) = Y(k) · exp(-j 2π k τ̂_samp / N),
%
%   and inverse-transformed.  No interpolator, loop filter, or NCO is
%   required.
%
%   Inputs
%     In    - input signal at 2 Sa/symbol (column vector, one polarisation)
%     NSymb - number of transmitted symbols (output limited to NSymb*2)
%     N     - FFT block size
%     beta  - pulse-shaping roll-off factor (0 < beta <= 1)
%
%   Output
%     Out   - clock-recovered signal (column vector, 2 Sa/symbol)

    eta   = 2;                            % input oversampling (Sa/symbol)
    shift = round((1 - 1/eta) * N);       % MG bin shift (N/2 for eta = 2)

    % MG summation bounds (eq. 5), converted to MATLAB 1-based indexing
    kLo = round((1 - beta) / (2*eta) * N) + 1;
    kHi = round((1 + beta) / (2*eta) * N);

    % Signed FFT bin index for the phase ramp: 0..N/2-1 then -N/2..-1
    k_idx = [0:N/2-1, -N/2:-1].';

    In      = In(:);
    LIn     = length(In);
    nBlocks = floor(LIn / N);

    Out = zeros(LIn, 1);

    for b = 1:nBlocks
        idx = (b-1)*N + (1:N);
        R   = fft(In(idx));

        % Feedforward MG estimate (units of T) → input-sample delay
        prod    = R(kLo:kHi) .* conj(R(kLo+shift:kHi+shift));
        tauT    = angle(sum(prod)) / (2*pi);
        tauSamp = tauT * eta;

        % Frequency-domain phase-ramp correction (eq. 1)
        R_corr   = R .* exp(-1j * 2*pi * k_idx * tauSamp / N);
        Out(idx) = ifft(R_corr);
    end

    % Pass-through any tail samples that don't fill a complete block
    tail = nBlocks*N + 1 : LIn;
    if ~isempty(tail)
        Out(tail) = In(tail);
    end

    % Limit output length to NSymb*2
    if NSymb*2 < length(Out)
        Out = Out(1:NSymb*2);
    end
end
