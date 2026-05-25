function Out = recovery_godard(In, NSymb, N, beta, mode, varargin)
%recovery_godard  Clock recovery using the Modified Godard estimator.
%
%   Out = recovery_godard(In, NSymb, N, beta)
%   Out = recovery_godard(In, NSymb, N, beta, 'feedforward')
%   Out = recovery_godard(In, NSymb, N, beta, 'feedforward', G)
%   Out = recovery_godard(In, NSymb, N, beta, 'feedback', ki, kp)
%   Out = recovery_godard(In, NSymb, N, beta, 'feedback', ki, kp, G)
%
%   Both modes share the Modified Godard frequency-domain timing metric
%   (Josten et al., Appl. Sci. 2017, eq. 5) computed over the
%   excess-bandwidth bins of an FFT block of N input samples,
%
%       S = Σ R(k) · R*(k + (1 - 1/η) N),   k ∈ [kLo, kHi].
%
%   --- Feedforward mode (default) ---
%   The per-block estimate
%
%       δ_b = arg(S) / (2π) · η                 (input samples, |δ_b| ≤ η/2)
%
%   is split into an integer slip and a fractional residual.  Integer
%   slips are absorbed by advancing the input read pointer
%
%       inStart_b = (b-1)·M + 1 + nSlip,    nSlip ← nSlip + round(δ_b),
%
%   so successive blocks see only the residual.  The full δ_b
%   (magnitude ≤ η/2 = 1 sample for η=2) is applied as a linear phase
%   ramp Y_corr(k) = Y(k)·exp(-j 2π k δ_b / N) followed by an IFFT.
%   Because the residual is bounded sub-sample regardless of cumulative
%   drift, the overlap-save guard G only needs to cover ~1 sample of
%   circular-shift wrap (cf. the old "G ≥ max cumulative |τ̂_samp|"
%   constraint, which scaled with record length).
%
%   --- Feedback mode ---
%   The correction is still a frequency-domain phase ramp, but the per-
%   block timing estimate is produced by a PI loop filter driven by the
%   small-angle approximation of arg(S),
%
%       e_b = imag(S_corr),
%
%   where S_corr is the metric evaluated on the already-corrected block
%   R_corr(k) = R(k)·exp(-j 2π k τ̂_samp / N).  Because the loop closes
%   on the residual error after correction, imag ≈ arg is self-consistent
%   once converged.  The integrator is seeded with the feedforward
%   arg(S)/(2π) estimate from the first block so that the small-angle
%   assumption is valid from block 1.  The PI update is
%
%       LF_I ← LF_I + ki·e_b,    τ̂_samp ← kp·e_b + LF_I,
%
%   followed by IFFT and overlap-save (same guard G as feedforward).
%
%   Inputs
%     In     - input signal at 2 Sa/symbol (column vector)
%     NSymb  - number of transmitted symbols (output limited to NSymb*2)
%     N      - FFT block size
%     beta   - pulse-shaping roll-off factor (0 < beta <= 1)
%     mode   - 'feedforward' (default) or 'feedback'
%   Feedforward extra args:
%     G      - overlap-save guard length per block edge (default N/4).
%              Must be < N/2 and ≥ ceil(η/2) samples (the maximum
%              fractional residual after integer-slip absorption);
%              independent of total record drift.
%   Feedback extra args:
%     ki, kp - integral and proportional gains of the PI loop filter.
%              τ̂_samp is in input-sample units, e_b is unnormalised
%              imag(S_corr), so gains absorb the |S| scale.
%     G      - overlap-save guard length per block edge (default N/4).
%              Same role as feedforward.
%
%   Output
%     Out    - clock-recovered signal (column vector, 2 Sa/symbol)

    if nargin < 5 || isempty(mode)
        mode = 'feedforward';
    end

    eta   = 2;                            % input oversampling (Sa/symbol)
    shift = round((1 - 1/eta) * N);       % MG bin shift (N/2 for eta = 2)
    kLo   = round((1 - beta) / (2*eta) * N) + 1;
    kHi   = round((1 + beta) / (2*eta) * N);

    In  = In(:);
    LIn = length(In);

    switch lower(mode)
        case 'feedforward'
            if isempty(varargin) || isempty(varargin{1})
                G = round(N / 4);
            else
                G = varargin{1};
            end
            Out = godard_ff(In, LIn, N, eta, shift, kLo, kHi, G);

        case 'feedback'
            if numel(varargin) < 2
                error('recovery_godard:MissingGains', ...
                      'Feedback mode requires ki and kp.');
            end
            ki = varargin{1};
            kp = varargin{2};
            if numel(varargin) < 3 || isempty(varargin{3})
                G = round(N / 4);
            else
                G = varargin{3};
            end
            Out = godard_fb(In, LIn, N, eta, shift, kLo, kHi, ki, kp, G);

        otherwise
            error('recovery_godard:UnknownMode', ...
                  'mode must be ''feedforward'' or ''feedback''.');
    end

    if NSymb*2 < length(Out)
        Out = Out(1:NSymb*2);
    end
end


function Out = godard_ff(In, LIn, N, eta, shift, kLo, kHi, G)
%godard_ff  Feedforward Modified Godard: integer-sample read-pointer slip
%           plus fractional-residual phase-ramp correction.

    k_idx = [0:N/2-1, -N/2:-1].';       % signed FFT bin index
    M     = N - 2*G;                    % valid output samples per block

    InPad      = [zeros(G,1); In; zeros(G,1)];
    LPad       = length(InPad);
    nBlocksMax = floor((LPad - N) / M) + 1;

    Out   = zeros(nBlocksMax * M, 1);
    nSlip = 0;                          % cumulative integer-sample slip
    nOut  = 0;

    for b = 1:nBlocksMax
        inStart = (b - 1) * M + 1 + nSlip;

        % Stop when the read window (or slip) leaves the padded input.
        if inStart < 1 || inStart + N - 1 > LPad
            break;
        end

        block = InPad(inStart : inStart + N - 1);
        R     = fft(block);

        prod  = R(kLo:kHi) .* conj(R(kLo+shift:kHi+shift));
        delta = (angle(sum(prod)) / (2*pi)) * eta;   % residual [samples]

        % Apply the full residual as a fractional phase ramp.  |delta|
        % ≤ η/2, so the circular-shift wrap stays inside G.
        R_corr = R .* exp(-1j * 2*pi * k_idx * delta / N);
        outBlk = ifft(R_corr);

        Out(nOut + 1 : nOut + M) = outBlk(G + 1 : G + M);
        nOut = nOut + M;

        % Absorb the integer part of this block's residual into the
        % read pointer for the next block.
        nSlip = nSlip + round(delta);
    end

    Out = Out(1 : nOut);
    Out = Out(1 : min(end, LIn));
end


function Out = godard_fb(In, LIn, N, eta, shift, kLo, kHi, ki, kp, G)
%godard_fb  Feedback Modified Godard: PI loop filter on imag(S_corr),
%           frequency-domain phase-ramp correction with overlap-save.

    k_idx = [0:N/2-1, -N/2:-1].';
    M     = N - 2*G;

    InPad   = [zeros(G,1); In; zeros(G,1)];
    LPad    = length(InPad);
    nBlocks = floor((LPad - N) / M) + 1;

    Out = zeros(nBlocks * M, 1);

    % Seed the integrator with the feedforward estimate from block 1 so
    % imag(S) ≈ arg(S) holds from the first PI update.
    block1 = InPad(1:N);
    R1     = fft(block1);
    S1     = sum(R1(kLo:kHi) .* conj(R1(kLo+shift:kHi+shift)));
    tauSamp = (angle(S1) / (2*pi)) * eta;
    LF_I    = tauSamp;

    for b = 1:nBlocks
        inStart = (b - 1) * M + 1;
        block   = InPad(inStart : inStart + N - 1);
        R       = fft(block);

        % Apply current cumulative phase-ramp correction
        R_corr = R .* exp(-1j * 2*pi * k_idx * tauSamp / N);

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
end
