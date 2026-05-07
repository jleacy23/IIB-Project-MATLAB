function E = transmitter_shot(SNR, B, lambda_m, K, alpha_dbkm, L_km, M, eta)
    % Wall-plug transmitter energy per bit [J/bit] for a shot-noise-
    % limited downstream PON link (no in-line amplification).
    %
    % Implements the SNR <-> P_tx relation of eq. for shot-noise SNR,
    % the optical-energy form eq:tx_energy, and the wall-plug form
    % eq:tx_energy_true of report/energy/energy.tex.
    %
    % SNR        : linear required SNR  (vectorised allowed)
    % B          : signal bandwidth [Hz] = symbol rate Rs under Nyquist
    % lambda_m   : optical carrier wavelength [m]
    % K          : passive splitter ratio (number of ONUs)
    % alpha_dbkm : fibre attenuation [dB/km]
    % L_km       : fibre length [km]
    % M          : modulation order
    % eta        : laser wall-plug efficiency in (0, 1]

    h  = 6.62607015e-34;
    c0 = 299792458;
    nu = c0 / lambda_m;

    P_tx = 2 * K * h * nu * B * 10.^(alpha_dbkm * L_km / 10) .* SNR;
    E    = P_tx ./ (B * log2(M) * eta);
end
