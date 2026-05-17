function y = equalize_fxp(x, SpS, NTaps, Mu, SingleSpike, N1, NOut, SignOnly, UpdateStep, T, PLanes) %#codegen
%equalize_fxp  Fixed-point adaptive butterfly equalization (CMA), parallel-lane.
%
%   y = equalize_fxp(x, SpS, NTaps, Mu, SingleSpike, N1, NOut, SignOnly, UpdateStep, T, PLanes)
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
%     PLanes        - (optional) number of parallel lanes (default 1 =
%                     serial).  The symbol stream is split into PLanes
%                     overlapping buffers (each successive buffer shifted
%                     by one symbol).  All lanes in a block are filtered
%                     with the same (frozen) weights and each produces one
%                     output sample exactly as in the serial case.  The
%                     per-tap gradient terms are summed over the lanes in
%                     the block (no averaging) and applied as a single
%                     weight update at the end of the block.
%
%   With PLanes = 1 this is bit-identical to the serial fixed-point CMA
%   (including the UpdateStep gating).
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

    %% Pre-allocate per-tap gradient accumulators (T.grad precision).
    %  These sum the gradient terms over the lanes of a block before the
    %  single (double-scaled) weight update.
    g1V  = complex(zeros(NTaps, 1, 'like', T.grad));
    g1H  = complex(zeros(NTaps, 1, 'like', T.grad));
    g2V  = complex(zeros(NTaps, 1, 'like', T.grad));
    g2H  = complex(zeros(NTaps, 1, 'like', T.grad));

    %% ====================================================================
    %  Adaptive equalisation loop (block-parallel over PLanes lanes)
    %  ====================================================================
    for iStart = 1:PLanes:OutLength
        % Lanes processed in this block (the last block may be partial).
        iEnd = min(iStart + PLanes - 1, OutLength);

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

            % --- CMA error and conjugate-output factor ---
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

            % --- Accumulate per-tap gradient over the lanes (T.grad) ---
            for k = 1:NTaps
                g1V(k) = g1V(k) + xv_i(k) * err1 * yc1;
                g1H(k) = g1H(k) + xh_i(k) * err1 * yc1;
                g2V(k) = g2V(k) + xv_i(k) * err2 * yc2;
                g2H(k) = g2H(k) + xh_i(k) * err2 * yc2;
            end

            % Flag the block in which the y-pol reinitialisation falls
            if i == N1 && SingleSpike
                doReinit = true;
            end
        end

        % --- Single CMA weight update for the block (gated by UpdateStep) ---
        %  Gradient summed over the block lanes in T.grad, then scaled by
        %  mu_dbl in double to avoid small mu rounding to zero in fxp.
        %  Gating on iEnd reproduces the serial mod(i,UpdateStep) behaviour
        %  when PLanes == 1.
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
