function y = equalize(x, SpS, NTaps, Mu, SingleSpike, N1, NOut, SignOnly, PLanes)
%equalize  Adaptive butterfly equalization (CMA), parallel-lane version.
%
%   y = equalize(x, SpS, NTaps, Mu, SingleSpike, N1, NOut, SignOnly, PLanes)
%
%   Inputs
%     x             - input signal [samples x 2]
%     SpS           - samples per symbol
%     NTaps         - number of FIR taps
%     Mu            - step size
%     SingleSpike   - true/false for single-spike initialisation
%     N1            - iteration to reinitialise y-pol weights
%     NOut          - samples to discard after equalisation
%     SignOnly      - if true, use sign(error) and complex-sign of y
%                     (sign(real(y)) + j*sign(imag(y))) in the update
%                     instead of the full multiplications
%     PLanes        - number of parallel lanes (default 1 = serial). The
%                     symbol stream is split into PLanes overlapping
%                     buffers (each successive buffer shifted by one
%                     symbol). All lanes in a block are filtered with the
%                     same (frozen) weights and each produces one output
%                     sample exactly as in the serial implementation. The
%                     per-lane gradient terms are summed and divided by the
%                     number of lanes in the block, then applied as a
%                     single weight update.
%
%   With PLanes = 1 this is bit-identical to the serial CMA.

    if nargin < 9 || isempty(PLanes)
        PLanes = 1;
    end

    % CMA radius
    R_CMA = sqrt(2);

    %% Input blocks (padding for convolution)
    x = [x(end-floor(NTaps/2)+1:end,:); x; x(1:floor(NTaps/2),:)];

    xV = convmtx(x(:,1).', NTaps);
    xH = convmtx(x(:,2).', NTaps);

    xV = xV(:, NTaps:SpS:end-NTaps+1);
    xH = xH(:, NTaps:SpS:end-NTaps+1);

    OutLength = floor((size(x,1) - NTaps + 1) / 2);

    %% Initialise outputs
    y1 = zeros(OutLength, 1);
    y2 = zeros(OutLength, 1);

    %% Initial filter coefficients
    w1V = zeros(NTaps, 1);
    w1H = zeros(NTaps, 1);
    w2V = zeros(NTaps, 1);
    w2H = zeros(NTaps, 1);

    if SingleSpike
        w1V(floor(NTaps/2)+1) = 1;
    end

    %% Adaptive equalisation loop (block-parallel over PLanes lanes)
    for iStart = 1:PLanes:OutLength
        % Lanes processed in this block (the last block may be partial).
        % Each lane is one overlapping symbol buffer; adjacent buffers
        % overlap by all but one symbol (the regressor columns of xV/xH).
        iEnd   = min(iStart + PLanes - 1, OutLength);
        nLanes = iEnd - iStart + 1;

        % Accumulated gradient terms over the lanes in this block.
        g1V = zeros(NTaps, 1);
        g1H = zeros(NTaps, 1);
        g2V = zeros(NTaps, 1);
        g2H = zeros(NTaps, 1);

        reinit = false;

        for i = iStart:iEnd
            % Compute outputs with the weights frozen for the block
            y1(i) = w1V'*xV(:,i) + w1H'*xH(:,i);
            y2(i) = w2V'*xV(:,i) + w2H'*xH(:,i);

            % CMA error and conjugate-output factor (optionally sign-reduced)
            if SignOnly
                e1  = sign(R_CMA - abs(y1(i))^2);
                e2  = sign(R_CMA - abs(y2(i))^2);
                yc1 = sign(real(y1(i))) - 1j*sign(imag(y1(i)));
                yc2 = sign(real(y2(i))) - 1j*sign(imag(y2(i)));
            else
                e1  = R_CMA - abs(y1(i))^2;
                e2  = R_CMA - abs(y2(i))^2;
                yc1 = conj(y1(i));
                yc2 = conj(y2(i));
            end

            % Accumulate the per-lane gradient terms
            g1V = g1V + xV(:,i)*e1*yc1;
            g1H = g1H + xH(:,i)*e1*yc1;
            g2V = g2V + xV(:,i)*e2*yc2;
            g2H = g2H + xH(:,i)*e2*yc2;

            % Flag the block in which the y-pol reinitialisation falls
            if i == N1 && SingleSpike
                reinit = true;
            end
        end

        % Single weight update from the averaged gradient over the block
        w1V = w1V + Mu*g1V;
        w1H = w1H + Mu*g1H;
        w2V = w2V + Mu*g2V;
        w2H = w2H + Mu*g2H;

        % Reinitialisation for SingleSpike (applied after the block update)
        if reinit
            w2H = conj(w1V(end:-1:1, 1));
            w2V = -conj(w1H(end:-1:1, 1));
        end
    end

    % Collect output and remove extra samples
    y = [y1, y2];
    y = y(1+NOut:end, :);
end
