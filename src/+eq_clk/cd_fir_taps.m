function g = cd_fir_taps(D, L, CLambda, Rs, SpSIn, NTap)
%CD_FIR_TAPS  Chromatic-dispersion compensating FIR (time-domain chirp).
%
%   g = cd_fir_taps(D, L, CLambda, Rs, SpSIn, NTap)
%
%   Returns the (column-vector) sampled chirp impulse response that
%   compensates a fibre with dispersion D [ps/(nm km)] over length L [km]
%   at centre wavelength CLambda [nm] and symbol rate Rs [GBd], sampled at
%   SpSIn Sa/symbol.  The symmetric window has 2*floor((NTap-1)/2)+1 taps:
%
%       g[m] = (Ts/sqrt(i*A)) * exp(i*pi*(m*Ts)^2/A),   m = -K..K
%
%   with A = CLambda^2 * D * L / c, Ts = 1/(SpSIn * Rs).  See equalize_td
%   for the reference implementation; this helper centralises the tap
%   computation so the combined blocks can reuse it.

    c = 299792458;

    D_si       = D * 1e-6;
    L_si       = L * 1e3;
    CLambda_si = CLambda * 1e-9;
    Rs_si      = Rs * 1e9;

    Ts = 1 / (SpSIn * Rs_si);
    A  = CLambda_si^2 * D_si * L_si / c;

    if A == 0
        g = 1;
        return;
    end

    K     = floor((NTap - 1) / 2);
    m     = (-K:K).';
    alpha = Ts / sqrt(1i * A);
    g     = alpha * exp(1i * pi * (m * Ts).^2 / A);
end
