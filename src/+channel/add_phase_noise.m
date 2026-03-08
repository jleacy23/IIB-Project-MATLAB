function Y = add_phase_noise(X, Rs, Linewidth)
%ADD_PHASE_NOISE  Apply laser phase noise (Wiener process).
%
%   Y = add_phase_noise(X, Rs, Linewidth)
%
%   Inputs
%     X         - input signal [samples x N_pol]
%     Rs        - symbol rate [GBd]
%     Linewidth - laser linewidth [Hz]

    Rs = Rs * 1e9;  % GBd -> Bd

    var_phi   = 2*pi * Linewidth / Rs;
    delta_phi = sqrt(var_phi) * randn(size(X,1), 1);
    phi       = cumsum(delta_phi);
    Y         = X .* exp(1i*phi);
end
