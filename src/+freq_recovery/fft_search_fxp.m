function [y, frequency_offset] = fft_search_fxp(x, training, Rs, Nfft, po2Twiddle, CordicIts, max_freq, T, data_aided, D) %#codegen
%FFT_SEARCH_FXP  Fixed-point FFT-based frequency offset estimator.
%
%   [y, frequency_offset] = fft_search_fxp(x, training, Rs, Nfft, po2Twiddle, CordicIts, T)
%
%   The received training sequence is de-rotated against the known
%   sequence in the angle domain (no complex multipliers):
%
%     phi_z(k) = angle(x(k)) - angle(training(k))   via CORDIC
%
%   A unit-amplitude complex sequence is then formed:
%
%     z_unit(k) = cordicrotate( phi_z(k), 1+0j )
%
%   and zero-padded to Nfft (must be a power of 2 >= L).  The fixed-point
%   FFT (fft.fft_fxp) is applied.  The peak bin is found by loop-based
%   magnitude comparison, and a fine frequency correction is computed
%   using the Jacobsen interpolator on the three neighbours.  The
%   frequency estimate is averaged across polarisations and the
%   correction applied sample-by-sample with cordicrotate.
%
%   Inputs
%     x          - input subframe  [Nsym x NPol]  (fi or castable to T.x)
%     training   - known training symbols [L x NPol]
%     Rs         - symbol rate [GBd]  (double scalar)
%     Nfft       - FFT size (integer power of 2, >= L)
%     po2Twiddle - logical: round twiddle factors to powers of 2 for fft_fxp
%     CordicIts  - number of CORDIC iterations
%     max_freq   - maximum allowed frequency offset normalized 
%     T          - fixed-point types table from freq_recovery.fxp_types
%     data_aided - true = training-aided (default), false = blind 4th-power
%     D          - number of data symbols for blind mode (required when
%                  data_aided = false)
%     max_freq   - maximum allowed frequency offset normalized to Rs
%
%   Outputs
%     y                - frequency-corrected subframe [Nsym x NPol], type T.x
%     frequency_offset - estimated frequency offset [Hz]  (double)

    if nargin < 7 || isempty(max_freq)
        max_freq = 1.0;
    end
    if max_freq <= 0
        error('freq_recovery:fft_search_fxp:BadMaxFreq', ...
              'max_freq must be positive.');
    end
    if nargin < 8 || isempty(T)
        T = freq_recovery.fxp_types('fixed16');
    end
    if nargin < 9, data_aided = true; end

    %% ----------------------------------------------------------------
    %  Fixed-point constants
    %% ----------------------------------------------------------------
    ZERO_TH    = cast(0,   'like', T.theta);
    UNIT_RE    = cast(1,   'like', T.acc);
    ZERO_ACC   = cast(0,   'like', T.acc);
    CORDIC_ITS = coder.const(CordicIts);

    %% ----------------------------------------------------------------
    %  Dimensions
    %% ----------------------------------------------------------------
    [L, N_pol] = size(training);
    Nsym       = size(x, 1);
    Nfft_c     = Nfft;
    %% ----------------------------------------------------------------
    %  FFT types — derived from T so the same WL/FL propagates into
    %  the butterfly without a separate hardcoded configuration.
    %% ----------------------------------------------------------------
    T_fft.x   = T.x;
    T_fft.tw  = T.x;
    T_fft.acc = T.acc;

    %% ----------------------------------------------------------------
    %  Pre-compute training phases (training-aided mode only)
    %% ----------------------------------------------------------------
    training_fi = cast(training, 'like', T.x);
    phi_tr = zeros(L, N_pol, 'like', T.theta);
    if data_aided
        for p = 1:N_pol
            for k = 1:L
                phi_tr(k, p) = cast(cordicangle(training_fi(k, p), CORDIC_ITS), ...
                                     'like', T.theta);
            end
        end
    end

    %% ----------------------------------------------------------------
    %  Cast input
    %% ----------------------------------------------------------------
    x_fi = cast(x, 'like', T.x);

    %% ----------------------------------------------------------------
    %  Compute frequency estimate per polarisation
    %% ----------------------------------------------------------------
    f_per_pol = zeros(1, N_pol);   % double Hz

    for p = 1:N_pol

        %% Build z_pad: unit-circle sequence, zero-padded to Nfft
        z_pad = complex(zeros(Nfft_c, 1, 'like', T.acc));

        if data_aided
            %% Training-aided: phi_z(k) = angle(x(k)) - angle(training(k))
            for k = 1:L
                phi_x_k  = cast(cordicangle(x_fi(k, p), CORDIC_ITS), 'like', T.theta);
                phi_z_k  = phi_x_k - phi_tr(k, p);

                unit_in  = complex(UNIT_RE, ZERO_ACC);
                z_pad(k) = cast(cordicrotate(phi_z_k, unit_in, CORDIC_ITS), 'like', T.acc);
            end
        else
            %% Blind: phi_z(k) = 4*angle(x_data(k)),  x_data = x(L+1..L+D)
            for k = 1:D
                phi_x_k  = cast(cordicangle(x_fi(L+k, p), CORDIC_ITS), 'like', T.theta);
                phi_z_k  = cast(mod(4.0 * double(phi_x_k), 2*pi) - pi, 'like', T.theta);

                unit_in  = complex(UNIT_RE, ZERO_ACC);
                z_pad(k) = cast(cordicrotate(phi_z_k, unit_in, CORDIC_ITS), 'like', T.acc);
            end
        end
        % Bins No+1 .. Nfft are already zero (zero-padding)

        %% Fixed-point FFT
        Z = fft.fft_fxp(z_pad, false, po2Twiddle, T_fft);   % [Nfft x 1]

        %% Find peak bin by loop (no vectorised max for codegen clarity)
        peak_mag_sq = -1.0;
        peak_bin_1  = 1;   % 1-indexed

        for k = 1:Nfft_c
            Zk_re  = double(real(Z(k)));
            Zk_im  = double(imag(Z(k)));
            mag_sq = Zk_re*Zk_re + Zk_im*Zk_im;
            if mag_sq > peak_mag_sq
                peak_mag_sq = mag_sq;
                peak_bin_1  = k;
            end
        end

        %% Fine interpolation (Jacobsen estimator) — in double
        %  Bins are 1-indexed here; wrap at spectrum edges
        km1 = mod(peak_bin_1 - 2, Nfft_c) + 1;
        kp1 = mod(peak_bin_1,     Nfft_c) + 1;

        Xm = complex(double(real(Z(km1))), double(imag(Z(km1))));
        Xk = complex(double(real(Z(peak_bin_1))), double(imag(Z(peak_bin_1))));
        Xp = complex(double(real(Z(kp1))), double(imag(Z(kp1))));

        denom = 2.0*Xk - Xm - Xp;
        if abs(denom) > 0.0
            delta = -real((Xp - Xm) / denom);
        else
            delta = 0.0;
        end

        %% Convert peak bin (1-indexed) to signed 0-indexed bin
        %  peak_bin_0 in {0 .. Nfft-1}; frequencies > Nyquist fold to negative
        peak_bin_0 = double(peak_bin_1) - 1.0 + delta;
        if peak_bin_0 >= double(Nfft_c) / 2.0
            peak_bin_0 = peak_bin_0 - double(Nfft_c);
        end

        f_per_pol(p) = peak_bin_0 / double(Nfft_c) * (Rs * 1e9);
    end

    %% ----------------------------------------------------------------
    %  Average across polarisations
    %% ----------------------------------------------------------------
    frequency_offset_Hz = 0.0;
    for p = 1:N_pol
        frequency_offset_Hz = frequency_offset_Hz + f_per_pol(p);
    end
    frequency_offset_Hz = frequency_offset_Hz / double(N_pol);
    if ~data_aided
        frequency_offset_Hz = frequency_offset_Hz / 4.0;
    end

    %% ----------------------------------------------------------------
    %  Phase correction: keep scaled phase in T.theta, then cast back to
    %  double and apply exp(+j*theta) in floating point.
    %  delta_theta already carries the negative sign for derotation.
    %% ----------------------------------------------------------------
    fprintf('Estimated Frequency Offset = %f for FL = %f \n', frequency_offset_Hz, T.x.FractionLength);
    delta_theta = cast(-2.0 * pi * frequency_offset_Hz / (Rs * 1e9 * max_freq), 'like', T.theta);
    y = complex(zeros(Nsym, N_pol, 'like', T.x));
    x_float = double(x_fi);
    theta_wrap = pi / max_freq;

    for p = 1:N_pol
        theta_fi = ZERO_TH;
        for i = 1:Nsym
            theta = double(theta_fi) * max_freq;
            y(i, p) = cast(x_float(i, p) * exp(1j * theta), 'like', T.x);

            % Explicit phase wrap in scaled domain. Do not rely on fi
            % overflow, which wraps at numeric range rather than 2*pi.
            theta_next = double(theta_fi + delta_theta);
            theta_next = mod(theta_next + theta_wrap, 2 * theta_wrap) - theta_wrap;
            theta_fi = cast(theta_next, 'like', T.theta);
        end
    end

    frequency_offset = frequency_offset_Hz;
end
