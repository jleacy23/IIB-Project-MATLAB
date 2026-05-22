function y = equalize_fxp(x, SpS, NTaps, Mu, SingleSpike, N1, NOut, SignOnly, UpdateStep, T, PLanes, Mode, Pilots, BlockLen) %#codegen
%equalize_fxp  Fixed-point adaptive butterfly equalization (CMA / pilot-aided).
%
%   y = equalize_fxp(x, SpS, NTaps, Mu, SingleSpike, N1, NOut, SignOnly, ...
%                    UpdateStep, T, PLanes, Mode, Pilots, BlockLen)
%
%   Inputs
%     x             - input signal [samples x 2] (fi or double)
%     SpS           - samples per symbol
%     NTaps         - number of FIR taps
%     Mu            - step size (fi or double)
%     SingleSpike   - true/false for single-spike initialisation
%     N1            - iteration to reinitialise y-pol weights
%     NOut          - samples to discard after equalisation
%     SignOnly      - if true, use sign(error) and complex-sign of y
%                     (sign(real(y)) + j*sign(imag(y))) in the update
%                     instead of the full multiplications
%     UpdateStep    - (optional) number of output samples between weight
%                     updates.  1 (default) updates every block.
%     T             - (optional) fixed-point types table from
%                     equalize_fxp_types.  If omitted, calls
%                     equalize_fxp_types('fixed16').
%     PLanes        - (optional) number of parallel lanes (default 1).
%                     Retained for backward compatibility; when BlockLen is
%                     omitted it sets the weight-update block length.
%     Mode          - (optional) update mode (default 0):
%                       0 = CMA          (blind constant-modulus)
%                       1 = pilot-aided  (data-aided LMS on pilots)
%     Pilots        - (optional) [NBlocks x 2] known pilot symbols, one per
%                     block, at the equalizer output scale (e.g. the CPON
%                     +/-3+/-3j pilots tiled per subframe).  Pilots(b,:) is
%                     the pilot at the first symbol of block b.  Empty
%                     (default) => no pilots, pure CMA over every symbol.
%     BlockLen      - (optional) weight-update block length in symbols
%                     (default = PLanes).  For CPON set BlockLen = 32:
%                     weights are held across the 32 symbols of a block and
%                     a single update is applied at the block end.  The
%                     pilot sits at the first symbol of each block.
%
%   CPON adaptation (see docs/cpon_framing_structure.md)
%     Symbols are processed in blocks of BlockLen (= 32 for CPON), all
%     filtered with the same frozen weights, with one weight update per
%     block (gated by UpdateStep):
%       * Mode 0 (CMA): the pilot symbol does NOT contribute to the update;
%         only the data symbols drive the CMA gradient.
%       * Mode 1 (pilot-aided): the error e = Pilots(b,:) - y(pilot) is
%         taken from the pilot of each polarisation, and its LMS gradient
%         updates the block-shared weights.
%     Only the first symbol of every block is treated as a pilot.  The
%     pilot reference is carried at T.y precision.
%
%   With Mode = 0, empty Pilots and BlockLen = PLanes this is bit-identical
%   to the previous parallel-lane fixed-point CMA (including UpdateStep).
%
%   The types table T must supply prototype fi objects for:
%     T.x      - input signal type
%     T.w      - filter coefficient type (needs enough FL for mu*grad)
%     T.y      - equalizer output type
%     T.acc    - accumulator type (inner product)
%     T.err    - error signal type
%     T.grad   - unscaled gradient type  x*err*conj(y)  before mu scaling
%     T.R_CMA  - CMA radius type
%
%   Weight update
%     The per-tap gradient x*err*conj(y), accumulated over the lanes of a
%     block in T.grad fixed-point precision, is converted to double, scaled
%     by double(Mu), and cast back to T.w.  This prevents very small mu
%     values from being rounded to zero inside the fixed-point datapath.
%
%   Best practices applied (per MathWorks Fixed-Point Designer manual):
%     - Data type definitions are separated from the algorithm via a types
%       table (cast/zeros ...'like'...).
%     - Subscripted assignment (:)= is used everywhere inside the loop to
%       prevent bit growth and preserve declared types.
%     - fimath uses SpecifyPrecision for products and sums so every
%       arithmetic result is truncated to the same WL/FL — no rescaling,
%       mimicking a uniform fixed-point datapath (FPGA / ASIC).
%     - convmtx (not fixed-point friendly) is replaced with explicit
%       indexing into a tap-delay line.

    %% Defaults for optional trailing arguments
    if nargin < 9 || isempty(UpdateStep)
        UpdateStep = 1;
    end
    if nargin < 10 || isempty(T)
        T = equalize_fxp_types('fixed16');
    end
    if nargin < 11 || isempty(PLanes)
        PLanes = 1;
    end
    if nargin < 12 || isempty(Mode)
        Mode = 0;
    end
    if nargin < 13
        Pilots = [];
    end
    if nargin < 14 || isempty(BlockLen)
        BlockLen = PLanes;
    end

    %% Pilot handling.  Carry the pilot reference at T.y precision so the
    %  data-aided error e = pilot - y is formed in the output datapath type.
    usePilots     = ~isempty(Pilots);
    NBlocksPilots = size(Pilots, 1);
    Pilots_fi     = cast(Pilots, 'like', T.y);

    %% Cast CMA radius
    R_CMA = cast(sqrt(2), 'like', T.R_CMA);

    %% Step size kept as double so very small values are not rounded to zero
    mu_dbl = double(Mu);

    %% Circular-pad input and cast to fixed-point
    halfTaps = floor(NTaps/2);
    xPad = [x(end-halfTaps+1:end, :); x; x(1:halfTaps, :)];
    xPad = cast(xPad, 'like', T.x);

    Nsamples = size(xPad, 1);
    OutLength = floor((Nsamples - NTaps + 1) / SpS);

    %% Build tap-delay-line index matrix (replaces convmtx)
    %  Each column j holds the NTaps indices into xPad that correspond to
    %  the j-th output sample.
    tapIdx = zeros(NTaps, OutLength);
    for j = 1:OutLength
        startIdx = (j-1)*SpS + 1;
        tapIdx(:, j) = (startIdx : startIdx + NTaps - 1).';
    end

    %% Initialise outputs (complex — equaliser outputs are complex-valued)
    y1 = complex(zeros(OutLength, 1, 'like', T.y));
    y2 = complex(zeros(OutLength, 1, 'like', T.y));

    %% Initial filter coefficients (complex — updated with complex gradients)
    w1V = complex(zeros(NTaps, 1, 'like', T.w));
    w1H = complex(zeros(NTaps, 1, 'like', T.w));
    w2V = complex(zeros(NTaps, 1, 'like', T.w));
    w2H = complex(zeros(NTaps, 1, 'like', T.w));

    if SingleSpike
        w1V(halfTaps + 1) = cast(1, 'like', T.w);
    end

    %% Pre-allocate temporaries (typed once, reused via (:)= )
    acc1 = cast(complex(0, 0), 'like', T.acc);
    acc2 = cast(complex(0, 0), 'like', T.acc);
    err1 = cast(0, 'like', T.err);
    err2 = cast(0, 'like', T.err);
    yc1  = complex(zeros(1, 1, 'like', T.y));
    yc2  = complex(zeros(1, 1, 'like', T.y));
    % Pilot-aided LMS error e = pilot - y and its gradient factor (T.y).
    ep1  = complex(zeros(1, 1, 'like', T.y));
    ep2  = complex(zeros(1, 1, 'like', T.y));
    fac1 = complex(zeros(1, 1, 'like', T.y));
    fac2 = complex(zeros(1, 1, 'like', T.y));

    %% Pre-allocate per-tap gradient accumulators (T.grad precision).
    %  These sum the gradient terms over the lanes of a block before the
    %  single (double-scaled) weight update.
    g1V  = complex(zeros(NTaps, 1, 'like', T.grad));
    g1H  = complex(zeros(NTaps, 1, 'like', T.grad));
    g2V  = complex(zeros(NTaps, 1, 'like', T.grad));
    g2H  = complex(zeros(NTaps, 1, 'like', T.grad));

    %% ====================================================================
    %  Adaptive equalisation loop (one BlockLen block per iteration)
    %  ====================================================================
    for iStart = 1:BlockLen:OutLength
        % Symbols processed in this block (the last block may be partial),
        % all filtered with the same frozen weights.
        iEnd = min(iStart + BlockLen - 1, OutLength);

        % CPON block index (1-based) for the pilot lookup.  iStart-1 is an
        % exact multiple of BlockLen, so b is integer.
        b = floor((iStart - 1) / BlockLen) + 1;

        % Zero the per-tap gradient accumulators for this block.
        g1V(:) = complex(0, 0);
        g1H(:) = complex(0, 0);
        g2V(:) = complex(0, 0);
        g2H(:) = complex(0, 0);

        doReinit = false;

        for i = iStart:iEnd
            % --- Extract tap vectors for this lane ---
            xv_i = xPad(tapIdx(:, i), 1);   % vertical   pol taps
            xh_i = xPad(tapIdx(:, i), 2);   % horizontal pol taps

            % --- Butterfly outputs (inner products, frozen weights) ---
            acc1(:) = complex(0, 0);
            acc2(:) = complex(0, 0);
            for k = 1:NTaps
                acc1(:) = acc1 + conj(w1V(k)) * xv_i(k) + conj(w1H(k)) * xh_i(k);
                acc2(:) = acc2 + conj(w2V(k)) * xv_i(k) + conj(w2H(k)) * xh_i(k);
            end
            y1(i) = acc1;
            y2(i) = acc2;

            % The pilot is the first symbol of each block (CPON: p = 1+32k).
            isPilot = usePilots && (i == iStart) && (b <= NBlocksPilots);

            if Mode == 1
                % --- Pilot-aided LMS: only the pilot drives the update ---
                if isPilot
                    ep1(:) = Pilots_fi(b, 1) - y1(i);   % data-aided error
                    ep2(:) = Pilots_fi(b, 2) - y2(i);
                    if SignOnly
                        fac1(:) = complex(sign(real(ep1)), -sign(imag(ep1)));
                        fac2(:) = complex(sign(real(ep2)), -sign(imag(ep2)));
                    else
                        fac1(:) = conj(ep1);
                        fac2(:) = conj(ep2);
                    end
                    for k = 1:NTaps
                        g1V(k) = g1V(k) + xv_i(k) * fac1;
                        g1H(k) = g1H(k) + xh_i(k) * fac1;
                        g2V(k) = g2V(k) + xv_i(k) * fac2;
                        g2H(k) = g2H(k) + xh_i(k) * fac2;
                    end
                end
            else
                % --- CMA: every symbol except the pilot contributes ---
                if ~isPilot
                    if SignOnly
                        err1(:) = sign(R_CMA - abs(y1(i))^2);
                        err2(:) = sign(R_CMA - abs(y2(i))^2);
                        yc1(:) = complex(sign(real(y1(i))), -sign(imag(y1(i))));
                        yc2(:) = complex(sign(real(y2(i))), -sign(imag(y2(i))));
                    else
                        err1(:) = (R_CMA - abs(y1(i))^2);
                        err2(:) = (R_CMA - abs(y2(i))^2);
                        yc1(:) = conj(y1(i));
                        yc2(:) = conj(y2(i));
                    end
                    for k = 1:NTaps
                        g1V(k) = g1V(k) + xv_i(k) * err1 * yc1;
                        g1H(k) = g1H(k) + xh_i(k) * err1 * yc1;
                        g2V(k) = g2V(k) + xv_i(k) * err2 * yc2;
                        g2H(k) = g2H(k) + xh_i(k) * err2 * yc2;
                    end
                end
            end

            % Flag the block in which the y-pol reinitialisation falls
            if i == N1 && SingleSpike
                doReinit = true;
            end
        end

        % --- Single weight update for the block (gated by UpdateStep) ---
        %  Gradient summed over the block (CMA: data symbols; pilot mode:
        %  the pilot) in T.grad, then scaled by mu_dbl in double to avoid
        %  small mu rounding to zero in fxp.  Gating on iEnd reproduces the
        %  serial mod(i,UpdateStep) behaviour when BlockLen == 1.
        if mod(iEnd, UpdateStep) == 0
            for k = 1:NTaps
                w1V(k) = w1V(k) + cast(mu_dbl * double(g1V(k)), 'like', T.w);
                w1H(k) = w1H(k) + cast(mu_dbl * double(g1H(k)), 'like', T.w);
                w2V(k) = w2V(k) + cast(mu_dbl * double(g2V(k)), 'like', T.w);
                w2H(k) = w2H(k) + cast(mu_dbl * double(g2H(k)), 'like', T.w);
            end
        end

        % --- Reinitialisation for SingleSpike (after the block update) ---
        if doReinit
            w2H(:) = conj(w1V(end:-1:1, 1));
            w2V(:) = -conj(w1H(end:-1:1, 1));
        end
    end

    %% Collect output and remove transient
    y = [y1, y2];
    y = y(1+NOut:end, :);
end
