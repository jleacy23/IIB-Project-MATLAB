function bits = qam_randomBits(Nbits, L, P, M)
%QAM_RANDOMBITS  Generate random binary bits.
%   Nbits  - number of bits to generate
%   L      - block length for pilot insertion
%   P      - number of pilot symbols at the start of every block
%   M      - QAM constellation order (must be a power of 2)

    % convery symbol legnths to bits
    L = L * log2(M);
    P = P * log2(M);
    Nbits = Nbits + mod(Nbits, L); % pad to multiple of block size
    bits = randi([0 1], Nbits, 1);
    pilot_bits = bits(1:P);

    % fill every block with pilot bits
    NBlocks = ceil(Nbits / L);
    for b = 1:NBlocks
        bits((b-1)*L+1:(b-1)*L+P) = pilot_bits;
    end
end
