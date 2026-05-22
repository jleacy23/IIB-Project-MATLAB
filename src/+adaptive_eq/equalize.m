function y = equalize(x, SpS, NTaps, Mu, SingleSpike, N1, NOut, SignOnly, PLanes, Mode, Pilots, BlockLen)
%equalize  Adaptive butterfly equalization (CMA / pilot-aided), block-parallel.
%
%   y = equalize(x, SpS, NTaps, Mu, SingleSpike, N1, NOut, SignOnly, ...
%                PLanes, Mode, Pilots, BlockLen)
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
%     PLanes        - number of parallel lanes (default 1 = serial). Kept
%                     for backward compatibility; when BlockLen is omitted
%                     it sets the weight-update block length (see BlockLen).
%     Mode          - update mode (default 0):
%                       0 = CMA          (blind constant-modulus)
%                       1 = pilot-aided  (data-aided LMS on pilots)
%     Pilots        - [NBlocks x 2] known pilot symbols, one per block, at
%                     the equalizer output scale (e.g. the +/-3+/-3j CPON
%                     pilots from modem.modulate, tiled per subframe).
%                     Pilots(b,:) is the pilot at the first symbol of
%                     block b.  Empty (default) => no pilots, pure CMA over
%                     every symbol (original behaviour).
%     BlockLen      - weight-update block length in symbols (default =
%                     PLanes).  For CPON set BlockLen = 32: weights are held
%                     across the 32 symbols of a block and updated once at
%                     the block end ("weights shared by the whole block").
%                     The pilot sits at the first symbol of each block.
%
%   CPON adaptation
%     The symbol stream is processed in blocks of BlockLen (= 32 for CPON).
%     All symbols in a block are filtered with the same frozen weights;
%     a single weight update is applied at the end of the block.
%       * Mode 0 (CMA): the pilot symbol does NOT contribute to the weight
%         update; only the data symbols of the block drive the CMA gradient.
%       * Mode 1 (pilot-aided): the error is computed only from the pilot of
%         each polarisation, e = Pilots(b,:) - y(pilot), and the resulting
%         LMS gradient updates the block-shared weights.
%     Note: only the first symbol of every block is treated as a pilot (per
%     the CPON rule p = 1 + 32k); the TS2..TS11 training symbols in block 1
%     are not special-cased here.
%
%   With Mode = 0, empty Pilots and BlockLen = PLanes this is bit-identical
%   to the previous parallel-lane CMA (and to the serial CMA at PLanes = 1).

    if nargin < 9 || isempty(PLanes)
        PLanes = 1;
    end
    if nargin < 10 || isempty(Mode)
        Mode = 0;
    end
    if nargin < 11
        Pilots = [];
    end
    if nargin < 12 || isempty(BlockLen)
        BlockLen = PLanes;
    end

    usePilots     = ~isempty(Pilots);
    NBlocksPilots = size(Pilots, 1);

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

    %% Adaptive equalisation loop (one BlockLen block per iteration)
    for iStart = 1:BlockLen:OutLength
        % Symbols processed in this block (the last block may be partial).
        % All are filtered with the same frozen weights; a single update is
        % applied at the block end so the weights are shared by the block.
        iEnd = min(iStart + BlockLen - 1, OutLength);

        % CPON block index (1-based) for the pilot lookup.  iStart-1 is an
        % exact multiple of BlockLen, so b is integer.
        b = floor((iStart - 1) / BlockLen) + 1;

        % Accumulated gradient terms over the symbols in this block.
        g1V = zeros(NTaps, 1);
        g1H = zeros(NTaps, 1);
        g2V = zeros(NTaps, 1);
        g2H = zeros(NTaps, 1);

        reinit = false;

        for i = iStart:iEnd
            % Compute outputs with the weights frozen for the block
            y1(i) = w1V'*xV(:,i) + w1H'*xH(:,i);
            y2(i) = w2V'*xV(:,i) + w2H'*xH(:,i);

            % The pilot is the first symbol of each block (CPON: p = 1+32k).
            isPilot = usePilots && (i == iStart) && (b <= NBlocksPilots);

            if Mode == 1
                % --- Pilot-aided LMS: only the pilot drives the update ---
                if isPilot
                    e1 = Pilots(b, 1) - y1(i);   % data-aided error
                    e2 = Pilots(b, 2) - y2(i);
                    if SignOnly
                        f1 = sign(real(e1)) - 1j*sign(imag(e1));
                        f2 = sign(real(e2)) - 1j*sign(imag(e2));
                    else
                        f1 = conj(e1);
                        f2 = conj(e2);
                    end
                    g1V = g1V + xV(:,i)*f1;
                    g1H = g1H + xH(:,i)*f1;
                    g2V = g2V + xV(:,i)*f2;
                    g2H = g2H + xH(:,i)*f2;
                end
            else
                % --- CMA: every symbol except the pilot contributes ---
                if ~isPilot
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
                    g1V = g1V + xV(:,i)*e1*yc1;
                    g1H = g1H + xH(:,i)*e1*yc1;
                    g2V = g2V + xV(:,i)*e2*yc2;
                    g2H = g2H + xH(:,i)*e2*yc2;
                end
            end

            % Flag the block in which the y-pol reinitialisation falls
            if i == N1 && SingleSpike
                reinit = true;
            end
        end

        % Single weight update from the summed block gradient
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
