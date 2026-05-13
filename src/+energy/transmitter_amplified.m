function E = transmitter_amplified(SNR, B, lambda_m, n_spans, L_span, ...
                                    alpha, beta_2, gamma_nl, NF_db, M, eta)
    % Wall-plug transmitter energy per bit [J/bit] for an amplified link
    % limited by ASE + non-linear interference.
    %
    % Solves eq:snr_nl of report/energy/energy.tex for the transmit power
    % P_tx (a depressed cubic) and combines with eq:tx_energy and
    % eq:tx_energy_true.  Returns the lower-power branch — the higher-
    % power root lies above the NLI-limited optimum.
    %
    % SNR      : linear required SNR (scalar)
    % B        : full WDM bandwidth [Hz]
    % lambda_m : optical carrier wavelength [m]
    % n_spans  : number of identical fibre spans (each with one EDFA)
    % L_span   : span length [km]
    % alpha    : fibre loss [nepers/km]
    % beta_2   : group velocity dispersion [s^2/km]
    % gamma_nl : fibre non-linear coefficient [W^-1 km^-1]
    % NF_db    : amplifier noise figure [dB]
    % M        : modulation order
    % eta      : laser wall-plug efficiency in (0, 1]

    h  = 6.62607015e-34;
    c0 = 299792458;
    nu = c0 / lambda_m;
    NF = 10^(NF_db / 10);

    % Amplifier gain compensates one span of fibre loss exactly.
    G     = exp(alpha * L_span);
    N_ASE = NF * h * nu * (G - 1);

    % Effective length and GN-model NLI coefficient (eq:snr_nl).  Units
    % are consistent in km: gamma [W^-1 km^-1], L_eff [km], alpha [km^-1],
    % beta_2 [s^2/km] ⇒ C_NLI [W^-2 Hz^2].  Argument of the log is
    % dimensionless: (s^2/km)/(km^-1) · Hz^2 = s^2 · Hz^2.
    L_eff = (1 - exp(-alpha * L_span)) / alpha;
    C_NLI = 8 * gamma_nl^2 * L_eff^2 * alpha / (27 * pi * abs(beta_2)) ...
            * log(abs(beta_2) / alpha * pi^2 * B^2);

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
