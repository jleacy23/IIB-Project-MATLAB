function [y, frequency_offset] = differential_fxp(x, training, Rs, CordicIts, T, data_aided, D) %#codegen
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

    %% ----------------------------------------------------------------
    %  Fixed-point constants
    %% ----------------------------------------------------------------
    ZERO_TH  = cast(0, 'like', T.theta);
    ZERO_ACC = cast(0, 'like', T.acc);
    CORDIC_ITS = coder.const(CordicIts);
    PI_TH    = cast(pi,   'like', T.theta);
    TWOPI_TH = cast(2*pi, 'like', T.theta);
    PI_VAL   = cast(pi,   'like', T.theta);
    PI_OVER2 = cast(pi/2, 'like', T.theta);

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
                phi_x = cast(cordicangle(x_fi(k, p), CORDIC_ITS), 'like', T.theta);
                phi_t = cast(cordicangle(training_fi(k, p), CORDIC_ITS), 'like', T.theta);
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
                phi_x = cast(cordicangle(x_fi(L + k, p), CORDIC_ITS), 'like', T.theta);
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
    %  Phase correction: linear ramp applied sample-by-sample via CORDIC
    %% ----------------------------------------------------------------
    delta_theta = cast(-2.0 * pi * frequency_offset_Hz / Rs_Hz, 'like', T.theta);
    y = complex(zeros(Nsym, N_pol, 'like', T.x));

    for p = 1:N_pol
        theta_fi = ZERO_TH;
        for i = 1:Nsym
            theta_d = mod(double(theta_fi) + pi, 2*pi) - pi;
            s_in = x_fi(i, p);
            if theta_d > pi/2
                theta_d = theta_d - pi;  s_in = -s_in;
            elseif theta_d < -pi/2
                theta_d = theta_d + pi;  s_in = -s_in;
            end
            theta_safe = cast(theta_d, 'like', T.theta);
            if theta_safe > PI_OVER2
                theta_safe = theta_safe - PI_VAL;  s_in = -s_in;
            elseif theta_safe < -PI_OVER2
                theta_safe = theta_safe + PI_VAL;  s_in = -s_in;
            end
            y(i, p)  = cast(cordicrotate(theta_safe, s_in, CORDIC_ITS), 'like', T.x);
            theta_fi = theta_fi + delta_theta;
        end
    end

    frequency_offset = frequency_offset_Hz
end
