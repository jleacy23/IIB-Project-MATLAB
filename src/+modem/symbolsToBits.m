function bits = symbolsToBits(symbols, M)
%SYMBOLSTOBITS  Demodulate QAM symbols to bits.
%
%   bits = symbolsToBits(symbols, M)

    k = log2(M);

    idx = qamdemod(symbols(:), M, ...
        'UnitAveragePower', true, ...
        'OutputType', 'integer');

    bitsMat = de2bi(idx, k, 'left-msb');
    bits = reshape(bitsMat.', [], 1);
end
