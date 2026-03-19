function [y, frequency_offset] = tretter_kay_fxp(x, training, Rs, CordicIts, T, data_aided, D, max_freq) %#codegen
%TRETTER_KAY_FXP  Fixed-point Tretter/Kay weighted phase-difference estimator.
%
%   [y, frequency_offset] = tretter_kay_fxp(x, training, Rs, CordicIts)
%   [y, frequency_offset] = tretter_kay_fxp(x, training, Rs, CordicIts, T)
%   [y, frequency_offset] = tretter_kay_fxp(x, training, Rs, CordicIts, T, data_aided, D)
%
%   Fixed-point equivalent of freq_recovery.tretter_kay.  Angles are
%   computed via CORDIC; arithmetic uses the type table T from
%   freq_recovery.fxp_types.  All loops are explicit for MATLAB Coder.
%
%   Algorithm
%     1. Build phase of de-rotated observation (same as differential_fxp).
%     2. Form successive phase differences:
%          dphi(k) = phi_z(k+1) - phi_z(k)
%     3. Apply Kay optimal weights (cast to T.acc):
%          w(k) = 6*k*(No-k) / (No*(No^2-1))
%     4. Weighted sum accumulated in T.acc (fixed-point); scaled to Hz.
%     5. Average across polarisations; apply phase ramp via CORDIC.
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
        error('freq_recovery:tretter_kay_fxp:BadMaxFreq', ...
              'max_freq must be positive.');
    end

    %% ----------------------------------------------------------------
    %  Fixed-point constants
    %% ----------------------------------------------------------------
    ZERO_TH  = cast(0,   'like', T.theta);
    ZERO_ACC = cast(0,   'like', T.acc);
    CORDIC_ITS = coder.const(CordicIts);
    TWOPI_TH = cast(2*pi, 'like', T.theta);
    PI_TH    = cast(pi,   'like', T.theta);

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
    %  Step 2 & 3 – Weighted phase-difference sum
    %
    %  w(k) = 6*k*(No-k) / (No*(No^2-1)),  k = 1 .. No-1
    %  w_k computed in double then cast to T.acc; dphi stays in T.theta;
    %  product and accumulation performed in T.acc.
    %% ----------------------------------------------------------------
    No_d       = double(No);
    w_denom    = No_d * (No_d * No_d - 1.0);   % No*(No^2-1)
    f_per_pol  = zeros(1, N_pol);               % double [Hz]
    Rs_Hz      = Rs * 1e9;
    f_coeff    = Rs_Hz / (2.0 * pi);

    for p = 1:N_pol
        w_acc = ZERO_ACC;   % T.acc fixed-point accumulator
        for k = 1:No - 1
            k_d    = double(k);
            w_k_fi = cast(6.0 * k_d * (No_d - k_d) / w_denom, 'like', T.acc);
            dphi   = phi_z(k + 1, p) - phi_z(k, p);
            if dphi > PI_TH
                dphi = dphi - TWOPI_TH;
            elseif dphi < -PI_TH
                dphi = dphi + TWOPI_TH;
            end
            dphi_acc = cast(dphi, 'like', T.acc);
            w_acc    = w_acc + cast(w_k_fi * dphi_acc, 'like', T.acc);
        end
        f_per_pol(p) = f_coeff * double(w_acc);
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
    %  Phase correction: keep scaled phase in T.theta, then cast back to
    %  double and apply exp(+j*theta) in floating point.
    %% ----------------------------------------------------------------
    delta_theta = cast(-2.0 * pi * frequency_offset_Hz / (Rs_Hz * max_freq), 'like', T.theta);
    y = complex(zeros(Nsym, N_pol, 'like', T.x));
    x_float = double(x_fi);
    theta_wrap = pi / max_freq;

    for p = 1:N_pol
        theta_fi = ZERO_TH;
        for i = 1:Nsym
            theta = double(theta_fi) * max_freq;
            y(i, p) = cast(x_float(i, p) * exp(1j * theta), 'like', T.x);

            theta_next = double(theta_fi + delta_theta);
            theta_next = mod(theta_next + theta_wrap, 2 * theta_wrap) - theta_wrap;
            theta_fi = cast(theta_next, 'like', T.theta);
        end
    end

    frequency_offset = frequency_offset_Hz;
end
