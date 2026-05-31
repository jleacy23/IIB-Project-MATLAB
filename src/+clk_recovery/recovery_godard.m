function Out = recovery_godard(In, NSymb, N, beta, ki, kp, G)
%recovery_godard  Feedback clock recovery using the Modified Godard
%   frequency-domain timing estimator with a PI loop filter.
%
%   Out = recovery_godard(In, NSymb, N, beta)
%   Out = recovery_godard(In, NSymb, N, beta, ki, kp)
%   Out = recovery_godard(In, NSymb, N, beta, ki, kp, G)
%
%   The Modified Godard frequency-domain timing metric (Josten et al.,
%   Appl. Sci. 2017, eq. 5) is computed over the excess-bandwidth bins
%   of an FFT block of N input samples,
%
%       S = Σ R(k) · R*(k + (1 - 1/η) N),   k ∈ [kLo, kHi].
%
%   The accumulated timing estimate τ̂_samp is split NCO-style into an
%   integer and a fractional part,
%
%       d_int = round(τ̂_samp),    μ = τ̂_samp - d_int  ∈ [-0.5, 0.5),
%
%   so that the correction is a frequency-domain phase ramp carrying only
%   the bounded fractional part,
%
%       R_corr(k) = R(k) · exp(-j 2π k μ / N),
%
%   while the integer part d_int is folded into the block read pointer
%   (the window is read d_int samples earlier).  This keeps the phase-ramp
%   slope bounded to ±½ sample no matter how far the clock drifts, so the
%   ramp never wraps samples circularly inside the FFT block — the
%   frequency-domain analogue of an NCO whose fractional interval feeds a
%   Farrow interpolator (cf. clk_recovery.recovery).
%
%   The per-block timing estimate τ̂_samp is produced by a PI loop filter
%   driven by the small-angle approximation of arg(S) on the already-
%   corrected block,
%
%       e_b = imag(S_corr),
%
%   where S_corr is the MG metric evaluated on R_corr.  Because the loop
%   closes on the residual error after correction, imag ≈ arg is self-
%   consistent once converged.  The integrator is initialised at zero
%   and acquires the steady-state offset via the integral path.  The PI
%   update is
%
%       LF_I ← LF_I + ki·e_b,    τ̂_samp ← kp·e_b + LF_I,
%
%   followed by IFFT and overlap-save with guard G.
%
%   Inputs
%     In     - input signal at 2 Sa/symbol (column vector)
%     NSymb  - number of transmitted symbols (output limited to NSymb*2)
%     N      - FFT block size
%     beta   - pulse-shaping roll-off factor (0 < beta <= 1)
%     ki, kp - integral and proportional gains of the PI loop filter.
%              τ̂_samp is in input-sample units, e_b is unnormalised
%              imag(S_corr), so gains absorb the |S| scale.
%              Defaults: ki = 1e-5, kp = 1e-4.
%     G      - overlap-save guard length per block edge (default N/4).
%              Must be < N/2.
%
%   Output
%     Out    - clock-recovered signal (column vector, 2 Sa/symbol)

    if nargin < 5 || isempty(ki)
        ki = 1e-5;
    end
    if nargin < 6 || isempty(kp)
        kp = 1e-4;
    end
    if nargin < 7 || isempty(G)
        G = round(N / 4);
    end

    eta   = 2;                            % input oversampling (Sa/symbol)
    shift = round((1 - 1/eta) * N);       % MG bin shift (N/2 for eta = 2)
    kLo   = round((1 - beta) / (2*eta) * N) + 1;
    kHi   = round((1 + beta) / (2*eta) * N);

    In  = In(:);
    LIn = length(In);

    k_idx = [0:N/2-1, -N/2:-1].';       % signed FFT bin index
    M     = N - 2*G;                    % valid output samples per block

    InPad   = [zeros(G,1); In; zeros(G,1)];
    LPad    = length(InPad);
    nBlocks = floor((LPad - N) / M) + 1;

    Out = zeros(nBlocks * M, 1);

    LF_I    = 0;
    tauSamp = 0;

    for b = 1:nBlocks
        % --- NCO-style integer/fractional split of the current estimate ---
        %  Integer part is consumed by the read pointer; only the bounded
        %  fractional residual mu drives the phase ramp.
        dInt = round(tauSamp);
        mu   = tauSamp - dInt;          % residual in [-0.5, 0.5)

        % --- Read window, integer timing folded into the read pointer -----
        %  rdStart drifts away from the nominal stride as the clock slips;
        %  indices outside the padded input are zero-filled.
        rdStart = (b - 1) * M + 1 - dInt;
        idx     = (rdStart : rdStart + N - 1).';
        valid   = idx >= 1 & idx <= LPad;
        block   = zeros(N, 1);
        block(valid) = InPad(idx(valid));

        R = fft(block);

        % Apply fractional-only phase-ramp correction (bounded slope)
        R_corr = R .* exp(-1j * 2*pi * k_idx * mu / N);

        % Residual error on the corrected block (small-angle proxy for arg)
        S = sum(R_corr(kLo:kHi) .* conj(R_corr(kLo+shift:kHi+shift)));
        e = imag(S);

        % PI loop filter
        LF_I    = LF_I + ki * e;
        tauSamp = kp * e + LF_I;

        outBlk = ifft(R_corr);

        outStart = (b - 1) * M + 1;
        Out(outStart : outStart + M - 1) = outBlk(G + 1 : G + M);
    end

    Out = Out(1 : min(end, LIn));

    if NSymb*2 < length(Out)
        Out = Out(1:NSymb*2);
    end
end
