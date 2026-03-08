function Y = add_chromatic_dispersion(X, L, SpS, Rs, D, CWavelength)
%ADD_CHROMATIC_DISPERSION  Apply chromatic dispersion.
%
%   Y = add_chromatic_dispersion(X, L, SpS, Rs, D, CWavelength)
%
%   Inputs
%     X           - input signal [samples x N_pol]
%     L           - fibre length [km]
%     SpS         - samples per symbol
%     Rs          - symbol rate [GBd]
%     D           - dispersion coefficient [ps/(nm*km)]
%     CWavelength - central wavelength [nm]

    c  = 299792458;
    L  = L * 1e3;           % km -> m
    Rs = Rs * 1e9;          % GBd -> Bd
    D  = D * 1e-6;          % ps/(nm*km) -> s/m^2
    CWavelength = CWavelength * 1e-9;  % nm -> m

    w = 2*pi*(-1/2:1/size(X,1):1/2-1/size(X,1)).' * SpS * Rs;
    G = exp(1i*((D*CWavelength^2)/(4*pi*c)) * L * w.^2);
    Y = ifft(ifftshift(G .* fftshift(fft(X))));
end
