function tx = rectPulse(symbols, sps)
%RECTPULSE  Rectangular pulse shaping (upsample + hold).
%
%   tx = rectPulse(symbols, sps)

    [Ns, Np] = size(symbols);
    tx = zeros(Ns*sps, Np);

    for p = 1:Np
        up = upsample(symbols(:,p), sps);
        h  = ones(sps,1);
        tx(:,p) = conv(up, h, 'same');
    end
end
