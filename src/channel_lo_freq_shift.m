function Y = channel_lo_freq_shift(X, DeltaF, Rs, SpS)
    % convert units
    deltaF = DeltaF * 1e6;  % MHz -> Hz
    rs     = Rs * 1e9;       % GBd -> Bd
    T = 1 / (rs * SpS);
    DeltaThetaF = 2 * pi * deltaF * T;
    k = repmat((0:size(X,1)-1).', 1, size(X,2));
    Y = X .* exp(1i * DeltaThetaF * k);
end