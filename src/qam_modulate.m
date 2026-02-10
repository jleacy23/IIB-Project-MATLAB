function symbols = qam_modulate(bits, M, N_pol)
%QAM_MODULATE  Map bits to QAM symbols.
%
%   symbols = qam_modulate(bits, M, N_pol)
%
%   Inputs
%     bits  - column vector of binary bits
%     M     - QAM constellation order (must be a power of 2)
%     N_pol - number of polarisations
%
%   Output
%     symbols - [Ns_pol x N_pol] complex QAM symbols (unit average power)

    if mod(log2(M),1) ~= 0
        error('M must be power of 2.');
    end

    k = log2(M);
    bits = bits(:);

    Ns_total = floor(length(bits) / k);
    Ns_pol   = floor(Ns_total / N_pol);

    if Ns_pol == 0
        error('Not enough bits.');
    end

    bits = bits(1:Ns_pol * N_pol * k);
    bits = reshape(bits, k, []).';
    symIdx = bi2de(bits, 'left-msb');

    syms = qammod(symIdx, M, ...
        'UnitAveragePower', true, ...
        'InputType', 'integer');

    symbols = reshape(syms, Ns_pol, N_pol);
end
