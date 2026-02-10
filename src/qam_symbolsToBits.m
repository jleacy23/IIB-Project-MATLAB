function bits = qam_symbolsToBits(symbols, M)
%QAM_SYMBOLSTOBITS  Demodulate QAM symbols to bits.
%
%   bits = qam_symbolsToBits(symbols, M)

    k = log2(M);

    idx = qamdemod(symbols(:), M, ...
        'UnitAveragePower', true, ...
        'OutputType', 'integer');

    bitsMat = de2bi(idx, k, 'left-msb');
    bits = reshape(bitsMat.', [], 1);
end
