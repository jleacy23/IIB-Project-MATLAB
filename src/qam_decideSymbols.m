function decidedSymbols = qam_decideSymbols(rxSymbols, M, N_pol)
%QAM_DECIDESYMBOLS  Hard-decision demodulation to nearest constellation point.
%
%   decidedSymbols = qam_decideSymbols(rxSymbols, M, N_pol)

    [Ns, Np] = size(rxSymbols);

    if Np ~= N_pol
        error('Polarization count mismatch.');
    end

    idx = qamdemod(rxSymbols(:), M, ...
        'UnitAveragePower', true, ...
        'OutputType', 'integer');

    decidedSymbols = qammod(idx, M, ...
        'UnitAveragePower', true, ...
        'InputType', 'integer');

    decidedSymbols = reshape(decidedSymbols, Ns, Np);
end
