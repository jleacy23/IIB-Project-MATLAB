function tx = qam_rrcPulse(symbols, sps, rolloff, span)
%QAM_RRCPULSE  Root-raised-cosine pulse shaping.
%
%   tx = qam_rrcPulse(symbols, sps, rolloff, span)

    [Ns, Np] = size(symbols);
    h  = rcosdesign(rolloff, span, sps, 'sqrt');
    tx = zeros(Ns*sps, Np);

    for p = 1:Np
        up = upsample(symbols(:,p), sps);
        tx(:,p) = conv(up, h, 'same');
    end
end
