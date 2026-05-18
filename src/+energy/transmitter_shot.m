function E = transmitter_shot(SNR, B, lambda_m, loss_dB, nsr0, K, M, eta)
    % Wall-plug transmitter energy per bit [J/bit] for a shot-noise-
    % limited downstream PON link (no in-line amplification).
    %
    % Implements the flat-loss shot-noise NSR relation of eq:snr_shot,
    % including the explicit 1:K split, the initial NSR contributed at
    % the OLT, and the per-bit energy of eq:tx_bit_energy:
    %
    %     NSR  = 2*K*Gamma*h*nu*B / P_tx + NSR_0
    %     P_tx = 2*K*Gamma*h*nu*B / (NSR - NSR_0)
    %     E_tx = P_tx / (K * B * log2(M))
    %          = 2*Gamma*h*nu / ((NSR - NSR_0) * log2(M))   [B = Rs]
    %
    % where NSR = 1/SNR, Gamma is the flat OLT-to-ONU power loss, NSR_0
    % is the noise-to-signal ratio already present at the OLT, and the
    % OLT power is shared across the K ONUs, as defined by the network
    % model of report/full/full.tex.  The K factor cancels in E_tx but
    % is kept explicit so the model matches the report term by term.
    %
    % SNR      : linear required SNR  (vectorised allowed)
    % B        : signal bandwidth [Hz] = symbol rate Rs under Nyquist
    % lambda_m : optical carrier wavelength [m]
    % loss_dB  : flat OLT-to-ONU power loss Gamma [dB]
    % nsr0     : initial (linear) NSR at the OLT
    % K        : passive splitter ratio (number of ONUs)
    % M        : modulation order
    % eta      : laser wall-plug efficiency in (0, 1]

    h  = 6.62607015e-34;
    c0 = 299792458;
    nu = c0 / lambda_m;

    Gamma = 10.^(loss_dB / 10);

    NSR = 1 ./ SNR;
    if any(NSR <= nsr0)
        error('transmitter_shot:infeasible', ...
            'Required NSR (%g) is below the OLT NSR floor (%g).', ...
            min(NSR), nsr0);
    end

    P_tx = 2 * K * Gamma * h * nu * B ./ (NSR - nsr0);
    E    = P_tx ./ (K * B * log2(M) * eta);
end
