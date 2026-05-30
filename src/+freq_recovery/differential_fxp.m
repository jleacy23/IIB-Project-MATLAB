function [y, frequency_offset] = differential_fxp(x, training, Rs, CordicIts, T, data_aided, D, max_freq) %#codegen
%DIFFERENTIAL_FXP  Fixed-point simple differential-phase frequency estimator.
%
%   [y, frequency_offset] = differential_fxp(x, training, Rs, CordicIts)
%   [y, frequency_offset] = differential_fxp(x, training, Rs, CordicIts, T)
%   [y, frequency_offset] = differential_fxp(x, training, Rs, CordicIts, T, data_aided, D)
%
%   Fixed-point equivalent of freq_recovery.differential.  Angles are
%   computed via CORDIC; arithmetic uses the type table T from
%   freq_recovery.fxp_types.  All loops are explicit so that MATLAB Coder
%   can generate straight C without dynamic allocation.
%
%   Algorithm
%     1. Build phase of de-rotated observation:
%          phi_z(k,p) = angle(x(k,p)) - angle(training(k,p))   [training-aided]
%          phi_z(k,p) = 4 * angle(x(L+k,p))                    [blind 4th-power]
%     2. Form successive differences:
%          dphi(k,p) = phi_z(k+1,p) - phi_z(k,p)
%     3. Estimate offset per polarisation:
%          f(p) = Rs / (2*pi) * (1/No) * sum_k dphi(k,p)
%     4. Average across polarisations and apply linear phase ramp via CORDIC.
%
%   Inputs
%     x          - input subframe  [Nsym x NPol]  (fi or castable to T.x)
%     training   - known training symbols [L x NPol]
%     Rs         - symbol rate [GBd]  (double scalar)
%     CordicIts  - number of CORDIC iterations  (integer)
%     T          - type table from freq_recovery.fxp_types (default 'fixed16')
%     data_aided - true = training-aided (default), false = blind 4th-power
%     D          - number of data symbols for blind mode
%
%   Outputs
%     y                - frequency-corrected subframe [Nsym x NPol], type T.x
%     frequency_offset - estimated frequency offset [Hz]  (double)

    if nargin < 5 || isempty(T)
        T = freq_recovery.fxp_types('fixed16');
    end
    if nargin < 6, data_aided = true; end
    if nargin < 8 || isempty(max_freq), max_freq = 1.0; end
    if max_freq <= 0
        error('freq_recovery:differential_fxp:BadMaxFreq', ...
              'max_freq must be positive.');
    end

    %% ----------------------------------------------------------------
    %  Fixed-point constants
    %% ----------------------------------------------------------------
    ZERO_TH  = cast(0, 'like', T.theta);
    ZERO_ACC = cast(0, 'like', T.acc);
    PI_TH    = cast(pi,   'like', T.theta);
    TWOPI_TH = cast(2*pi, 'like', T.theta);

    %% ----------------------------------------------------------------
    %  Dimensions
    %% ----------------------------------------------------------------
    [L, N_pol] = size(training);
    Nsym       = size(x, 1);
    if data_aided
        No = L;
    else
        No = D;
    end

    %% ----------------------------------------------------------------
    %  Cast inputs
    %% ----------------------------------------------------------------
    x_fi        = cast(x,        'like', T.x);
    training_fi = cast(training, 'like', T.x);

    %% ----------------------------------------------------------------
    %  Step 1 – Compute phi_z: phase of de-rotated observation sequence
    %           phi_z(k,p) stored as T.theta
    %% ----------------------------------------------------------------
    phi_z = zeros(No, N_pol, 'like', T.theta);

    for p = 1:N_pol
        if data_aided
            %% Training-aided: phi_z(k) = angle(x(k)) - angle(training(k))
            for k = 1:L
                phi_x = cordic.vectoring(real(x_fi(k, p)),        imag(x_fi(k, p)),        CordicIts, T);
                phi_t = cordic.vectoring(real(training_fi(k, p)), imag(training_fi(k, p)), CordicIts, T);
                phi_z(k, p) = phi_x - phi_t;
                if phi_z(k, p) > PI_TH
                    phi_z(k, p) = phi_z(k, p) - TWOPI_TH;
                elseif phi_z(k, p) < -PI_TH
                    phi_z(k, p) = phi_z(k, p) + TWOPI_TH;
                end
            end
        else
            %% Blind: phi_z(k) = 4 * angle(x_data(k))
            for k = 1:D
                phi_x = cordic.vectoring(real(x_fi(L+k, p)), imag(x_fi(L+k, p)), CordicIts, T);
                phi_z(k, p) = cast(mod(4.0 * double(phi_x), 2*pi) - pi, 'like', T.theta);
            end
        end
    end

    %% ----------------------------------------------------------------
    %  Step 2 & 3 – Successive phase differences, accumulate, scale
    %% ----------------------------------------------------------------
    f_per_pol = zeros(1, N_pol);   % double [Hz]
    Rs_Hz     = Rs * 1e9;
    f_coeff   = Rs_Hz / (2.0 * pi);

    for p = 1:N_pol
        acc = ZERO_ACC;
        for k = 1:No - 1
            dphi = phi_z(k + 1, p) - phi_z(k, p);
            if dphi > PI_TH
                dphi = dphi - TWOPI_TH;
            elseif dphi < -PI_TH
                dphi = dphi + TWOPI_TH;
            end
            dphi = cast(dphi, "like", T.acc);
            acc  = acc + dphi;
        end
        % Divide by No (observation length, matching the float version)
        f_per_pol(p) = f_coeff * double(acc) / double(No);
    end

    %% ----------------------------------------------------------------
    %  Blind mode: undo 4th-power scaling
    %% ----------------------------------------------------------------
    if ~data_aided
        for p = 1:N_pol
            f_per_pol(p) = f_per_pol(p) / 4.0;
        end
    end

    %% ----------------------------------------------------------------
    %  Average across polarisations
    %% ----------------------------------------------------------------
    frequency_offset_Hz = 0.0;
    for p = 1:N_pol
        frequency_offset_Hz = frequency_offset_Hz + f_per_pol(p);
    end
    frequency_offset_Hz = frequency_offset_Hz / double(N_pol);

    %% ----------------------------------------------------------------
    %  Phase correction: accumulate the scaled phase ramp in T.theta and
    %  apply the per-symbol rotation exp(+j*theta) via CORDIC.
    %% ----------------------------------------------------------------
    delta_theta = cast(-2.0 * pi * frequency_offset_Hz / (Rs_Hz * max_freq), 'like', T.theta);
    y = complex(zeros(Nsym, N_pol, 'like', T.x));
    theta_wrap = pi / max_freq;

    for p = 1:N_pol
        theta_fi = ZERO_TH;
        for i = 1:Nsym
            theta = double(theta_fi) * max_freq;
            [yr, yi] = cordic.rotate(real(x_fi(i, p)), imag(x_fi(i, p)), theta, CordicIts, T);
            y(i, p) = complex(yr, yi);

            % Add in double so the explicit phase wrap below is not
            % pre-empted by fi-overflow on theta_fi + delta_theta.
            theta_next = double(theta_fi) + double(delta_theta);
            theta_next = mod(theta_next + theta_wrap, 2 * theta_wrap) - theta_wrap;
            theta_fi = cast(theta_next, 'like', T.theta);
        end
    end

    frequency_offset = frequency_offset_Hz;
end
