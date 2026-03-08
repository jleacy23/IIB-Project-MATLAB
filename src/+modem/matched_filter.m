function rx = matched_filter(y, sps, pulse, rolloff, span)
%MATCHED_FILTER  Apply a matched filter to a received signal.
%
%   rx = matched_filter(y, sps, pulse)
%   rx = matched_filter(y, sps, pulse, rolloff, span)
%
%   Inputs
%     y       - received signal [samples x NPol]
%     sps     - samples per symbol
%     pulse   - 'rect' or 'rrc'
%     rolloff - RRC roll-off factor (required when pulse = 'rrc')
%     span    - RRC filter span in symbols (required when pulse = 'rrc')
%
%   Output
%     rx      - filtered signal [samples x NPol], same size as y

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
    else
        error('matched_filter:badPulse', ...
            'Unsupported pulse type ''%s''. Use ''rect'' or ''rrc''.', pulse);
    end

    rx = zeros(size(y));
    for p = 1:Np
        rx(:, p) = conv(y(:, p), h, 'same');
    end
end
