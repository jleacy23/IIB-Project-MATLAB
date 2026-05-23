function rx = matched_filter(y, sps, pulse, rolloff, span)
%MATCHED_FILTER  Apply a matched filter to a received signal.
%
%   rx = matched_filter(y, sps, pulse)
%   rx = matched_filter(y, sps, pulse, rolloff, span)
%
%   Inputs
%     y       - received signal [samples x NPol]
%     sps     - samples per symbol
%     pulse   - 'rect', 'rrc', or 'nyquist'
%     rolloff - roll-off factor (required when pulse = 'rrc' or 'nyquist')
%     span    - filter span in symbols (required when pulse = 'rrc' or 'nyquist')
%
%   Output
%     rx      - filtered signal [samples x NPol], same size as y
%
%   Notes
%     'rrc'     — root-raised-cosine; matched to an rrcPulse transmitter so
%                 the combined response is a (Nyquist) raised cosine.
%     'nyquist' — full raised-cosine receive filter, matched to a
%                 nyquistPulse transmitter.  Since that transmit pulse is
%                 already ISI-free, this receive filter is optional (it
%                 maximises SNR / rejects out-of-band noise but is not
%                 required to satisfy the Nyquist criterion).

    Np = size(y, 2);

    if strcmp(pulse, 'rect')
        h = ones(sps, 1) / sps;
    elseif strcmp(pulse, 'rrc')
        if nargin < 5
            error('matched_filter:missingArgs', ...
                'rolloff and span are required for RRC matched filter.');
        end
        h = rcosdesign(rolloff, span, sps, 'sqrt').';
        h = h / sum(abs(h).^2);   % normalise for unit energy
    elseif strcmp(pulse, 'nyquist')
        if nargin < 5
            error('matched_filter:missingArgs', ...
                'rolloff and span are required for Nyquist matched filter.');
        end
        h = rcosdesign(rolloff, span, sps, 'normal').';
        h = h / sum(abs(h).^2);   % normalise for unit energy
    else
        error('matched_filter:badPulse', ...
            'Unsupported pulse type ''%s''. Use ''rect'', ''rrc'', or ''nyquist''.', pulse);
    end

    rx = zeros(size(y));
    for p = 1:Np
        rx(:, p) = conv(y(:, p), h, 'same');
    end
end
