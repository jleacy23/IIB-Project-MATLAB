function [y, frequency_offset] = differential_kay_fxp(x, training, Rs, CordicIts, T, data_aided, D) %#codegen
%DIFFERENTIAL_KAY_FXP  Fixed-point two-stage differential + Tretter-Kay estimator.
%
%   [y, frequency_offset] = differential_kay_fxp(x, training, Rs, CordicIts)
%   [y, frequency_offset] = differential_kay_fxp(x, training, Rs, CordicIts, T)
%   [y, frequency_offset] = differential_kay_fxp(x, training, Rs, CordicIts, T, data_aided, D)
%
%   Fixed-point equivalent of freq_recovery.differential_kay.
%
%   Stage 1  Coarse:  differential_fxp obtains a coarse estimate f_coarse
%            and produces a coarsely corrected output x_coarse.
%
%   Stage 2  Fine:    tretter_kay_fxp runs on x_coarse to estimate the
%            residual offset f_fine.
%
%   The combined estimate f_coarse + f_fine is then applied to the
%   original input x in a single pass so that quantisation from the
%   intermediate correction does not accumulate.
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
    CORDIC_ITS = coder.const(CordicIts);

    Nsym  = size(x, 1);
    N_pol = size(x, 2);
    Rs_Hz = Rs * 1e9;

    x_fi = cast(x, 'like', T.x);

    %% ----------------------------------------------------------------
    %  Stage 1: Coarse estimate via differential_fxp
    %% ----------------------------------------------------------------
    if data_aided
        [x_coarse, f_coarse] = freq_recovery.differential_fxp( ...
            x_fi, training, Rs, CordicIts, T);
    else
        [x_coarse, f_coarse] = freq_recovery.differential_fxp( ...
            x_fi, training, Rs, CordicIts, T, false, D);
    end

    %% ----------------------------------------------------------------
    %  Stage 2: Fine residual estimate via tretter_kay_fxp
    %           Operates on the coarsely corrected signal
    %% ----------------------------------------------------------------
    if data_aided
        [~, f_fine] = freq_recovery.tretter_kay_fxp( ...
            x_coarse, training, Rs, CordicIts, T);
    else
        [~, f_fine] = freq_recovery.tretter_kay_fxp( ...
            x_coarse, training, Rs, CordicIts, T, false, D);
    end

    %% ----------------------------------------------------------------
    %  Combined estimate
    %% ----------------------------------------------------------------
    frequency_offset_Hz = f_coarse + f_fine;

    %% ----------------------------------------------------------------
    %  Apply total correction to the original input in a single pass
    %  (avoids accumulating quantisation error from the two-stage path)
    %% ----------------------------------------------------------------
    delta_theta = cast(-2.0 * pi * frequency_offset_Hz / Rs_Hz, 'like', T.theta);
    y = complex(zeros(Nsym, N_pol, 'like', T.x));

    for p = 1:N_pol
        theta_fi = ZERO_TH;
        for i = 1:Nsym
            y(i, p)  = cast(cordicrotate(theta_fi, x_fi(i, p), CORDIC_ITS), 'like', T.x);
            theta_fi = theta_fi + delta_theta;
        end
    end

    frequency_offset = frequency_offset_Hz;
end
