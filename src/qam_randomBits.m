function bits = qam_randomBits(Nbits)
%QAM_RANDOMBITS  Generate random binary bits.
    bits = randi([0 1], Nbits, 1);
end
