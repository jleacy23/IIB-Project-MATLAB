function [y, frequency_offset] = fitz_fxp(x, training, Rs, N, CordicIts, T) %#codegen
%FITZ_FXP  Fixed-point Fitz (1994) autocorrelation frequency estimator.
%
%   [y, frequency_offset] = fitz_fxp(x, training, Rs, N, CordicIts, T)
%
%   Computes the lag-N autocorrelation of z(k) = x(k)*conj(training(k))
%   without complex multipliers by working in the angle domain:
%
%     angle(z(k)*conj(z(k-N))) = angle(z(k)) - angle(z(k-N))
%                               = [angle(x(k)) - angle(tr(k))]
%                               - [angle(x(k-N)) - angle(tr(k-N))]
%
%   Both angles are computed via CORDIC.  The complex autocorrelation sum
%   is accumulated by projecting each phase-difference term onto a unit
%   circle via cordicrotate, then the angle of the resulting complex sum
%   is extracted with cordicangle.
%
%     R_unit = sum_{k=N+1}^{L} exp( j*(phi_z(k) - phi_z(k-N)) )
%     f0 = Rs / (2*pi*N) * angle(R_unit)
%
%   Inputs
%     x          - input subframe  [Nsym x NPol]  (fi or castable to T.x)
%     training   - known training symbols [L x NPol]
%     Rs         - symbol rate [GBd]  (double scalar)
%     N          - autocorrelation lag  (double scalar, 1 <= N < L)
%     CordicIts  - number of CORDIC iterations
%     T          - fixed-point types struct from freq_recovery.fxp_types
%
%   Outputs
%     y                - frequency-corrected subframe [Nsym x NPol], type T.x
%     frequency_offset - estimated frequency offset [kHz]  (double)

    if nargin < 6 || isempty(T)
        T = freq_recovery.fxp_types('fixed16');
    end

    %% ----------------------------------------------------------------
    %  Fixed-point constants
    %% ----------------------------------------------------------------
    ZERO_TH    = cast(0, 'like', T.theta);
    ZERO_ACC   = cast(0, 'like', T.acc);
    CORDIC_ITS = coder.const(CordicIts);

    %% ----------------------------------------------------------------
    %  Dimensions
    %% ----------------------------------------------------------------
    [L, N_pol] = size(training);
    Nsym       = size(x, 1);

    N_lag = N;   % make lag a compile-time constant

    %% ----------------------------------------------------------------
    %  Cast training and input
    %% ----------------------------------------------------------------
    training_fi = cast(training, 'like', T.x);

    %% ----------------------------------------------------------------
    %  Cast input
    %% ----------------------------------------------------------------
    x_fi = cast(x, 'like', T.x);

    %% ----------------------------------------------------------------
    %  Estimate: autocorrelation-of-phases approach, one angle per pol
    %% ----------------------------------------------------------------
    R_angle = zeros(1, N_pol, 'like', T.theta);

    for p = 1:N_pol

        %% Compute z(k) = x(k) * conj(training(k))  for k = 1..L
        z = complex(zeros(L, 1, 'like', T.acc));
        for k = 1:L
            z(k) = cast(x_fi(k, p) * conj(training_fi(k, p)), 'like', T.acc);
        end

        %% Accumulate autocorrelation  R = sum_{k=N+1}^{L} z(k)*conj(z(k-N))
        %  Each term is one complex multiply; cordicangle is called once on the sum.
        R_re = ZERO_ACC;
        R_im = ZERO_ACC;

        for k = N_lag+1 : L
            term = cast(z(k) * conj(z(k - N_lag)), 'like', T.acc);
            R_re = R_re + real(term);
            R_im = R_im + imag(term);
        end

        R_fi       = complex(cast(R_re, 'like', T.theta), cast(R_im, 'like', T.theta));
        R_angle(p) = cast(cordicangle(R_fi, CORDIC_ITS), 'like', T.theta);
    end

    %% ----------------------------------------------------------------
    %  Average angle across polarisations then convert to Hz in double
    %    f0_Hz = Rs_Hz / (2*pi*N) * mean(R_angle)
    %% ----------------------------------------------------------------
    R_angle_sum = ZERO_TH;
    for p = 1:N_pol
        R_angle_sum = R_angle_sum + R_angle(p);
    end
    R_angle_mean = double(R_angle_sum) / double(N_pol);

    Rs_Hz               = Rs * 1e9;
    frequency_offset_Hz = Rs_Hz / (2.0 * pi * double(N_lag)) * R_angle_mean;

    %% ----------------------------------------------------------------
    %  Phase correction: accumulate linear ramp and cordicrotate
    %% ----------------------------------------------------------------
    delta_theta = cast(-2.0 * pi * frequency_offset_Hz / Rs_Hz, 'like', T.theta);
    y = complex(zeros(Nsym, N_pol, 'like', T.x));

    for p = 1:N_pol
        theta_fi = ZERO_TH;
        for i = 1:Nsym
            y(i, p) = cast(cordicrotate(theta_fi, x_fi(i, p), CORDIC_ITS), 'like', T.x);
            theta_fi = theta_fi + delta_theta;
        end
    end

    frequency_offset = frequency_offset_Hz / 1e3;   % Hz -> kHz
end
