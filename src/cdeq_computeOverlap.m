function NOverlap = cdeq_computeOverlap(D, L, CLambda, Rs, SpSIn, NFFT)
%CDEQ_COMPUTEOVERLAP  Compute overlap length for overlap-save CD equalizer.
%
%   NOverlap = cdeq_computeOverlap(D, L, CLambda, Rs, SpSIn, NFFT)
%
%   Inputs (user-friendly units)
%     D       - dispersion coefficient [ps/(nm*km)]
%     L       - fibre length [km]
%     CLambda - central wavelength [nm]
%     Rs      - symbol rate [GBd]
%     SpSIn   - samples per symbol
%     NFFT    - FFT block size

    c = 299792458;

    % Unit conversion
    D       = D * 1e-6;
    L       = L * 1e3;
    CLambda = CLambda * 1e-9;
    Rs      = Rs * 1e9;

    beta2    = -D * CLambda^2 / (2*pi*c);
    NOverlap = ceil(6.67 * abs(beta2) * (Rs^2) * L * SpSIn);
    NOverlap = NOverlap + mod(NOverlap, 2);  % ensure even
    NOverlap = min(NOverlap, floor(NFFT/2));
end
