function Out = equalize_td_fxp(In, D, L, CLambda, Rs, NPol, SpSIn, T) %#codegen
%EQUALIZE_TD_FXP  Fixed-point time-domain (FIR) CD compensation.
%
%   Out = equalize_td_fxp(In, D, L, CLambda, Rs, NPol, SpSIn, T)
%
%   Alternate to equalize_fxp: instead of overlap-save frequency-domain
%   compensation, this convolves the input with the truncated chromatic
%   dispersion impulse response (a sampled chirp FIR).  The FIR realises
%   exactly the same all-pass transfer function used by equalize_fxp:
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
%   The Ts/sqrt(i*A) factor keeps the FIR (approximately) unit-gain, so the
%   signal amplitude is preserved through the fixed-point datapath.
%
%   Inputs (user-friendly units)
%     In       - input signal [samples x NPol] (fi or double)
%     D        - dispersion coefficient [ps/(nm*km)]
%     L        - fibre length [km]
%     CLambda  - central wavelength [nm]
%     Rs       - symbol rate [GBd]
%     NPol     - number of polarizations (1 or 2)
%     SpSIn    - samples per symbol
%     T        - (optional) fixed-point types table from
%                equalize_td_fxp_types.  If omitted, uses 'fixed32'.
%
%   The types table T must supply prototype fi objects for:
%     T.x    - input / output signal
%     T.hcd  - FIR tap coefficients
%     T.acc  - accumulator (multiply-accumulate)
%
%   Implementation notes for codegen:
%     - Taps are pre-computed in double and cast once to fi (T.hcd).
%     - The number of taps follows the standard chirp-support criterion
%       (Savory) sized to the signal band B = Rs: the chirp instantaneous
%       frequency reaches the band edge Rs/2 at |m| = K = |A|*Rs^2*SpS/2,
%       giving an odd FIR of 2K+1 taps (matches the report's N_CD).
%     - The input is treated as periodic (circular convolution), matching
%       the cyclic extension used by the overlap-save float reference and
%       avoiding edge transients.  Output length equals input length.
%     - Each polarisation is processed independently.
%     - All arithmetic uses SpecifyPrecision fimath so every product and
%       sum is truncated to the same WL/FL — no bit growth.

    %% Default types table
    if nargin < 8 || isempty(T)
        T = cd_eq.equalize_td_fxp_types('fixed32');
    end

    c = 299792458;

    %% Unit conversion (double — used only for tap computation)
    D_si       = D * 1e-6;
    L_si       = L * 1e3;
    CLambda_si = CLambda * 1e-9;
    Rs_si      = Rs * 1e9;

    Ts = 1 / (SpSIn * Rs_si);          % input sampling period [s]
    A  = CLambda_si^2 * D_si * L_si / c;

    NIn = size(In, 1);

    %% Trivial case: no dispersion → identity
    if A == 0
        Out = cast(In, 'like', T.x);
        return;
    end

    %% CD impulse response (FIR taps), pre-computed in double
    %  The filter is dimensioned to the signal bandwidth B = Rs (Nyquist
    %  pulse-shaping is assumed), not the full sampling Nyquist Fs/2.  The
    %  dispersive memory over that band is the delay spread
    %     Delta_tau/Ts = |A| * Rs^2 * SpS   [samples],
    %  matching the report's N_CD = Delta_tau/T + 1 (tab:cd_taps).  The
    %  impulse-response chirp has instantaneous frequency f_inst = m*Ts/A,
    %  which reaches the signal band edge Rs/2 at |m| = K, giving an odd FIR
    %  of 2K+1 taps.  (An oversampled signal carries no energy beyond Rs/2,
    %  so spanning the chirp out to Fs/2 would only add cost.)
    NspanSamp = abs(A) * Rs_si^2 * SpSIn;        % delay spread Delta_tau/Ts
    K     = floor(NspanSamp / 2);
    K     = min(K, floor((NIn - 1) / 2));        % cannot exceed signal length
    m     = (-K:K).';                        % tap lags (column)
    alpha = Ts / sqrt(1i * A);               % unit-gain normalisation
    g     = alpha * exp(1i * pi * (m * Ts).^2 / A);

    g_fi  = cast(g, 'like', T.hcd);          % cast taps to fixed-point
    nTap  = 2 * K + 1;

    %% Pre-allocate output (accumulator type)
    Out = complex(zeros(NIn, NPol, 'like', T.acc));

    InX = cast(In, 'like', T.x);

    %% ================================================================
    %  Per-polarisation circular FIR convolution
    %    y[n] = sum_{l=-K..K} g[l] * x[n - l]   (indices wrap mod NIn)
    %  ================================================================
    for pol = 1:NPol
        for n = 1:NIn
            acc = complex(cast(0, 'like', T.acc));
            for t = 1:nTap
                % tap t corresponds to lag l = m(t) = (t-1) - K
                idx = mod(n - (t - 1 - K) - 1, NIn) + 1;
                acc = acc + g_fi(t) * InX(idx, pol);
            end
            Out(n, pol) = acc;
        end
    end
end
