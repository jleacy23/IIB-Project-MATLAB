function decidedSymbols = decideSymbols(rxSymbols, M, N_pol)
%DECIDESYMBOLS  Hard-decision demodulation to nearest constellation point.
%
%   decidedSymbols = decideSymbols(rxSymbols, M, N_pol)

    [Ns, Np] = size(rxSymbols);

    if Np ~= N_pol
        error('Polarization count mismatch.');
    end

    idx = qamdemod(rxSymbols(:), M, ...
        'UnitAveragePower', false, ...
        'OutputType', 'integer');

    decidedSymbols = qammod(idx, M, ...
        'UnitAveragePower', false, ...
        'InputType', 'integer');

    decidedSymbols = reshape(decidedSymbols, Ns, Np);
end
