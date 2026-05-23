function tx = nyquistPulse(symbols, sps, rolloff, span)
%NYQUISTPULSE  Raised-cosine (Nyquist) pulse shaping.
%
%   tx = nyquistPulse(symbols, sps, rolloff, span)
%
%   Shapes the symbols with a full raised-cosine filter — the canonical
%   Nyquist pulse, which is inter-symbol-interference free at the symbol
%   instants by itself (unlike the root-raised-cosine of rrcPulse, which is
%   ISI-free only when paired with a matching receive filter).  Use this
%   where a band-limited transmit pulse is wanted but no separate matched
%   filter is applied at the receiver (e.g. equaliser characterisation),
%   matching the report's "Nyquist pulse-shaping" assumption.
%
%   Inputs
%     symbols - symbol stream [Ns x Np]
%     sps     - samples per symbol
%     rolloff - raised-cosine roll-off factor (0..1)
%     span    - filter span in symbols
%
%   Output
%     tx      - shaped signal [Ns*sps x Np]
%
%   The filter is scaled to unit peak so the symbol-rate samples recover the
%   transmitted symbol amplitudes (sampling tx(1:sps:end) is ISI-free).

    [Ns, Np] = size(symbols);
    h  = rcosdesign(rolloff, span, sps, 'normal');
    h  = h / max(h);                 % unit peak -> symbol amplitude preserved
    tx = zeros(Ns*sps, Np);

    for p = 1:Np
        up = upsample(symbols(:,p), sps);
        tx(:,p) = conv(up, h, 'same');
    end
end
