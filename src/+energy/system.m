function E = system(E_tx, E_rx, K)
    % System energy per bit [J/bit] for a 1:K downstream PON, summing the
    % wall-plug transmitter energy and the receiver energy across all K
    % ONUs.  Implements eq:sys_energy of report/energy/energy.tex.
    %
    % E_tx : wall-plug transmitter energy per bit (e.g. from
    %        energy.transmitter_shot or energy.transmitter_amplified)
    % E_rx : per-ONU receiver energy per bit (sum over stages of
    %        energy.receiver)
    % K    : number of ONUs (= passive splitter ratio)

    E = E_tx / K + E_rx;
end
