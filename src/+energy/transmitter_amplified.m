function E = transmitter_amplified(SNR, B, lambda_m, n_spans, NF_db, G_db, C_NLI, M, eta)
    % Wall-plug transmitter energy per bit [J/bit] for an amplified link
    % limited by ASE + non-linear interference.
    %
    % Solves eq:snr_nl of report/energy/energy.tex for the transmit power
    % P_tx (a depressed cubic) and combines with eq:tx_energy and
    % eq:tx_energy_true.  Returns the lower-power branch — the higher-
    % power root lies above the NLI-limited optimum.
    %
    % SNR     : linear required SNR (scalar)
    % B       : signal bandwidth [Hz] = symbol rate Rs under Nyquist
    % lambda_m: optical carrier wavelength [m]
    % n_spans : number of identical fibre spans (each with one EDFA)
    % NF_db   : amplifier noise figure [dB]
    % G_db    : amplifier gain [dB] (assumed to compensate one span loss)
    % C_NLI   : NLI coefficient [W^-2 Hz^2] from the GN model
    % M       : modulation order
    % eta     : laser wall-plug efficiency in (0, 1]

    h  = 6.62607015e-34;
    c0 = 299792458;
    nu = c0 / lambda_m;
    NF = 10^(NF_db / 10);
    G  = 10^(G_db  / 10);
    N_ASE = NF * h * nu * (G - 1);

    a = SNR * n_spans * C_NLI / B^2;
    c = SNR * n_spans * N_ASE * B;
    r = roots([a, 0, -1, c]);
    r = real(r(abs(imag(r)) < 1e-9 * max(abs(r), 1)));
    r = r(r > 0);
    if isempty(r)
        error('energy:transmitter_amplified:noSolution', ...
            'Required SNR exceeds NLI-limited maximum for these link parameters.');
    end
    P_tx = min(r);

    E = P_tx / (B * log2(M) * eta);
end
