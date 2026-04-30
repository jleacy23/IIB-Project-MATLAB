function E = receiver(NAdd, NMult, EAdd, EMult, M, Oversampling, n)
    % Calculates receiver energy per bit for one stage of the receiver.
    % NAdd, NMult: real adds/multiplies per symbol (per sample if oversampled)
    % EAdd, EMult: per-bit energy coefficients (E_A = EAdd*n, E_M = EMult*n^2)
    % M:           modulation order
    % Oversampling: f/Rs (use 1 for symbol-rate stages)
    % n:           word length in bits

    Es = 2 * Oversampling * (NAdd * EAdd * n + NMult * EMult * n^2);
    E  = Es / log2(M);
end
