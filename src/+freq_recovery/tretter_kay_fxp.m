function [y, frequency_offset] = tretter_kay_fxp(x, training, Rs, CordicIts, T) %#codegen
%TRETTER_KAY_FXP  Fixed-point Tretter/Kay frequency offset estimator.
%
%   [y, frequency_offset] = tretter_kay_fxp(x, training, Rs, CordicIts, T)
%
%   Multiplier-free derotation: the phase of z(k) = x(k)*conj(training(k))
%   is obtained without forming z(k) explicitly:
%
%     angle(z(k)) = angle(x(k)) - angle(training(k))
%
%   Both angles are computed via CORDIC.  The phase difference
%
%     angle(dz(k)) = angle(z(k+1)) - angle(z(k))
%                  = [angle(x(k+1)) - angle(tr(k+1))] - [angle(x(k)) - angle(tr(k))]
%
%   is then weighted by the Tretter-Kay coefficients and accumulated.
%
%   All phase operations use CORDIC; the output phase correction is
%   applied by accumulating a fixed phase increment and using
%   cordicrotate sample-by-sample.
%
%   Inputs
%     x          - input subframe  [Nsym x NPol]  (fi or castable to T.x)
%     training   - known training symbols [L x NPol]
%     Rs         - symbol rate [GBd]  (double scalar)
%     CordicIts  - number of CORDIC iterations (double scalar)
%     T          - fixed-point types struct from freq_recovery.fxp_types
%
%   Outputs
%     y                - frequency-corrected subframe [Nsym x NPol], type T.x
%     frequency_offset - estimated frequency offset [kHz]  (double)

    if nargin < 5 || isempty(T)
        T = freq_recovery.fxp_types('fixed16');
    end

    %% ----------------------------------------------------------------
    %  Fixed-point constants
    %% ----------------------------------------------------------------
    ZERO_ACC   = cast(0, 'like', T.acc);
    ZERO_TH    = cast(0, 'like', T.theta);
    CORDIC_ITS = coder.const(CordicIts);

    %% ----------------------------------------------------------------
    %  Dimensions
    %% ----------------------------------------------------------------
    [L, N_pol] = size(training);
    Nsym       = size(x, 1);

    %% ----------------------------------------------------------------
    %  Pre-compute Tretter-Kay weights  w(k) = 6k(L-k) / (L(L^2-1))
    %  These are compile-time constants for a fixed training length.
    %% ----------------------------------------------------------------
    w = zeros(L-1, 1, 'like', T.acc);
    denom = double(L) * (double(L)^2 - 1.0);
    for k = 1:L-1
        w(k) = cast(6.0 * double(k) * double(L - k) / denom, 'like', T.acc);
    end

    %% ----------------------------------------------------------------
    %  Pre-compute training phases
    %  phi_tr(k,p) = cordicangle(training(k,p))
    %  cordicangle output FL = input FL - 2; cast immediately to T.theta.
    %% ----------------------------------------------------------------
    training_fi = cast(training, 'like', T.x);
    phi_tr = zeros(L, N_pol, 'like', T.theta);
    for p = 1:N_pol
        for k = 1:L
            phi_tr(k, p) = cast(cordicangle(training_fi(k, p), CORDIC_ITS), ...
                                 'like', T.theta);
        end
    end

    %% ----------------------------------------------------------------
    %  Cast input
    %% ----------------------------------------------------------------
    x_fi = cast(x, 'like', T.x);

    %% ----------------------------------------------------------------
    %  Estimate: weighted phase-difference sum, one value per pol
    %% ----------------------------------------------------------------
    f_angle = zeros(1, N_pol, 'like', T.acc);  % angle-domain offset per pol [rad]

    for p = 1:N_pol

        %% Compute phi_z(k) = cordicangle(x(k)) - phi_tr(k)
        phi_z_prev = ZERO_TH;
        phi_z_curr = ZERO_TH;

        % k = 1: seed phi_z_prev
        phi_x_1  = cast(cordicangle(x_fi(1, p), CORDIC_ITS), 'like', T.theta);
        phi_z_prev = phi_x_1 - phi_tr(1, p);

        %% Accumulate weighted phase differences
        acc = ZERO_ACC;
        for k = 2:L
            phi_x_k  = cast(cordicangle(x_fi(k, p), CORDIC_ITS), 'like', T.theta);
            phi_z_curr = phi_x_k - phi_tr(k, p);

            angle_dz = phi_z_curr - phi_z_prev;       % phase diff, no multiply
            % wrap to [-pi, pi] to avoid large angle issues (e.g. > 2*pi offset)
            if angle_dz > pi
                angle_dz = angle_dz - 2*pi;
            elseif angle_dz < -pi
                angle_dz = angle_dz + 2*pi;
            end
            acc      = acc + cast(w(k-1), 'like', T.acc) * cast(angle_dz, 'like', T.acc);

            phi_z_prev = phi_z_curr;
        end

        f_angle(p) = acc;
    end

    %% ----------------------------------------------------------------
    %  Average offset across polarisations (in angle domain)
    %  Convert to Hz in double
    %    f0_Hz = Rs_Hz / (2*pi) * mean(f_angle)
    %----------------------------------------------------------------
    f_angle_sum = ZERO_ACC;
    for p = 1:N_pol
        f_angle_sum = f_angle_sum + f_angle(p);
    end
    f_angle_mean = double(f_angle_sum) / double(N_pol);

    Rs_Hz          = Rs * 1e9;
    frequency_offset_Hz = Rs_Hz / (2.0 * pi) * f_angle_mean;

    %% ----------------------------------------------------------------
    %  Phase correction: accumulate linear phase ramp and cordicrotate
    %  delta_theta = -2*pi*f0/Rs per sample  (constant increment)
    %% ----------------------------------------------------------------
    delta_theta    = cast(-2.0 * pi * frequency_offset_Hz / Rs_Hz, 'like', T.theta);
    y = complex(zeros(Nsym, N_pol, 'like', T.x));

    for p = 1:N_pol
        theta_fi = ZERO_TH;
        for i = 1:Nsym
            y(i, p) = cast(cordicrotate(theta_fi, x_fi(i, p), CORDIC_ITS), 'like', T.x);
            theta_fi = theta_fi + delta_theta;   % wraps naturally (Wrap overflow)
        end
    end

    frequency_offset = frequency_offset_Hz / 1e3;   % Hz -> kHz
end
