function bits = randomBits(Nbits)
%RANDOMBITS  Generate random binary data bits for DP-QPSK.
%
%   bits = randomBits(Nbits)
%
%   Input
%     Nbits - desired number of data bits (will be rounded up so that the
%             resulting subframe count is integer after modulate() pads)
%
%   Output
%     bits  - [Nbits x 1] column vector of random 0/1 values

    bits = randi([0 1], Nbits, 1);
end
