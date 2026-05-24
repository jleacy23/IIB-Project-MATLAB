function Out = equalize_td(In, D, L, CLambda, Rs, NPol, SpSIn, NTap)
%EQUALIZE_TD  Time-domain (FIR) CD compensation — floating-point reference.
%
%   Out = equalize_td(In, D, L, CLambda, Rs, NPol, SpSIn, NTap)
%
%   Floating-point reference for equalize_td_fxp.  Instead of overlap-save
%   frequency-domain compensation (see equalize.m), this convolves the
%   input with the truncated chromatic dispersion impulse response (a
%   sampled chirp FIR).  The FIR realises the same all-pass transfer
%   function as the frequency-domain equalizer:
%
%       H_CD(f) = exp(-i*pi * CLambda^2 * D * L / c * f^2)
%
%   whose continuous inverse Fourier transform is
%
%       h(t) = (1/sqrt(i*A)) * exp(i*pi * t^2 / A),   A = CLambda^2*D*L/c.
%
%   Sampled at the input rate (Ts = 1/(SpSIn*Rs)) and scaled by Ts so the
%   discrete convolution approximates the continuous one, the taps become
%
%       g[m] = (Ts/sqrt(i*A)) * exp(i*pi * (m*Ts)^2 / A),   m = -K..K.
%
%   Inputs (user-friendly units)
%     In       - input signal [samples x NPol]
%     D        - dispersion coefficient [ps/(nm*km)]
%     L        - fibre length [km]
%     CLambda  - central wavelength [nm]
%     Rs       - symbol rate [GBd]
%     NPol     - number of polarizations (1 or 2)
%     SpSIn    - samples per symbol
%     NTap     - FIR length (use cd_eq.computeOverlap to size to the signal
%                bandwidth).  The symmetric chirp window uses the largest
%                odd value not exceeding NTap, i.e. 2*floor((NTap-1)/2)+1.
%
%   The input is treated as periodic (circular convolution), matching the
%   cyclic extension used by the overlap-save reference and the fixed-point
%   implementation.  Output length equals input length.

    c = 299792458;

    % Unit conversion
    D_si       = D * 1e-6;
    L_si       = L * 1e3;
    CLambda_si = CLambda * 1e-9;
    Rs_si      = Rs * 1e9;

    Ts = 1 / (SpSIn * Rs_si);          % input sampling period [s]
    A  = CLambda_si^2 * D_si * L_si / c;

    NIn = size(In, 1);

    %% Trivial case: no dispersion → identity
    if A == 0
        Out = In;
        return;
    end

    %% CD impulse response (FIR taps)
    %  The FIR length is supplied by the caller (typically via
    %  cd_eq.computeOverlap, which sizes the support to the signal band).
    %  The symmetric chirp window has 2K+1 = closest odd value <= NTap.
    K     = floor((NTap - 1) / 2);
    K     = min(K, floor((NIn - 1) / 2));        % cannot exceed signal length
    m     = (-K:K).';                        % tap lags (column)
    alpha = Ts / sqrt(1i * A);               % unit-gain normalisation
    g     = alpha * exp(1i * pi * (m * Ts).^2 / A);

    %% Circular FIR convolution per polarisation
    %    y[n] = sum_{l=-K..K} g[l] * x[n - l]   (indices wrap mod NIn)
    Out = zeros(NIn, NPol);
    idxBase = (1:NIn).';
    for pol = 1:NPol
        acc = zeros(NIn, 1);
        for t = 1:numel(g)
            l   = m(t);                          % lag
            idx = mod(idxBase - l - 1, NIn) + 1; % circular shift
            acc = acc + g(t) * In(idx, pol);
        end
        Out(:, pol) = acc;
    end
end
