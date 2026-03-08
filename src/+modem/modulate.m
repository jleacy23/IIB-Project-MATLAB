function [symbols, pilots] = modulate(bits, M, N_pol, P)
%MODULATE  Map bits to QAM symbols.
%
%   symbols = modulate(bits, M, N_pol)
%
%   Inputs
%     bits  - column vector of binary bits
%     M     - QAM constellation order (must be a power of 2)
%     N_pol - number of polarisations
%     L     - block length for pilot insertion
%     P     - number of pilot symbols at the start of every block
%
%   Output
%     symbols - [Ns_pol x N_pol] complex QAM symbols (unit average power)
%     pilots  - [P x 1] pilot symbols

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
        'UnitAveragePower', false, ...
        'InputType', 'integer');
    pilots = syms(1:P); % this sequence is repeated every block over both polarisations.
    symbols = reshape(syms, Ns_pol, N_pol);
end
