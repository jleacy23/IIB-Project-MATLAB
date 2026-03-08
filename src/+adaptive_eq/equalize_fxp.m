function y = equalize_fxp(x, SpS, Eq, NTaps, Mu, SingleSpike, N1, N2, NOut, T) %#codegen
%equalize_fxp  Fixed-point adaptive butterfly equalization (CMA / RDE / CMA+RDE).
%
%   y = equalize_fxp(x, SpS, Eq, NTaps, Mu, SingleSpike, N1, N2, NOut, T)
%
%   Inputs
%     x             - input signal [samples x 2] (fi or double)
%     SpS           - samples per symbol
%     Eq            - algorithm: 'CMA', 'RDE', or 'CMA+RDE'
%     NTaps         - number of FIR taps
%     Mu            - step size (fi or double)
%     SingleSpike   - true/false for single-spike initialisation
%     N1            - iteration to reinitialise y-pol weights
%     N2            - iteration to switch CMA->RDE ([] if unused)
%     NOut          - samples to discard after equalisation
%     T             - (optional) fixed-point types table from
%                     equalize_fxp_types.  If omitted, calls
%                     equalize_fxp_types('fixed16').
%
%   The types table T must supply prototype fi objects for:
%     T.x      - input signal type
%     T.w      - filter coefficient type
%     T.y      - equalizer output type
%     T.acc    - accumulator type (inner product)
%     T.err    - error signal type
%     T.mu     - step-size type
%     T.R_CMA  - CMA radius type
%     T.R_RDE  - RDE radii type
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

    %% Default types table
    if nargin < 10 || isempty(T)
        T = equalize_fxp_types('fixed16');
    end

    %% Algorithm / mode flags (integer logic – no fi needed)
    CMAFlag  = false;
    RDEFlag  = false;
    CMAtoRDE = false;

    if strcmp(Eq, 'CMA')
        CMAFlag = true;
    elseif strcmp(Eq, 'RDE')
        RDEFlag = true;
    elseif strcmp(Eq, 'CMA+RDE')
        CMAFlag  = true;
        CMAtoRDE = true;
    else
        error('Unsupported equalizer type');
    end

    %% Cast constants to fixed-point types
    % CMA radius
    if CMAFlag
        if ~CMAtoRDE
            R_CMA = cast(sqrt(2), 'like', T.R_CMA);
        else
            R_CMA = cast(1.32, 'like', T.R_CMA);
        end
    else
        R_CMA = cast(0, 'like', T.R_CMA);   % placeholder
    end

    % RDE radii
    if CMAtoRDE || RDEFlag
        R_RDE = cast([1/sqrt(5), 1, 3/sqrt(5)], 'like', T.R_RDE);
    else
        R_RDE = cast([0 0 0], 'like', T.R_RDE);   % placeholder
    end

    %% Cast step size
    mu_fxp = cast(Mu, 'like', T.mu);

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

    %% ====================================================================
    %  Adaptive equalisation loop
    %  ====================================================================
    for i = 1:OutLength
        % --- Extract tap vectors for this iteration ---
        xv_i = xPad(tapIdx(:, i), 1);   % vertical   pol taps
        xh_i = xPad(tapIdx(:, i), 2);   % horizontal pol taps

        % --- Compute butterfly outputs (inner products) ---
        acc1(:) = complex(0, 0);
        acc2(:) = complex(0, 0);
        for k = 1:NTaps
            acc1(:) = acc1 + conj(w1V(k)) * xv_i(k) + conj(w1H(k)) * xh_i(k);
            acc2(:) = acc2 + conj(w2V(k)) * xv_i(k) + conj(w2H(k)) * xh_i(k);
        end
        y1(i) = acc1;
        y2(i) = acc2;

        % --- Coefficient update ---
        if CMAFlag
            % CMA error: e = (R - |y|^2)
            err1(:) = (R_CMA - abs(y1(i))^2);
            err2(:) = (R_CMA - abs(y2(i))^2);
            for k = 1:NTaps
                w1V(k) = w1V(k) + mu_fxp * xv_i(k) * err1 * conj(y1(i));
                w1H(k) = w1H(k) + mu_fxp * xh_i(k) * err1 * conj(y1(i));
                w2V(k) = w2V(k) + mu_fxp * xv_i(k) * err2 * conj(y2(i));
                w2H(k) = w2H(k) + mu_fxp * xh_i(k) * err2 * conj(y2(i));
            end

            % Switch CMA -> RDE
            if CMAtoRDE && i == N2
                CMAFlag = false;
                RDEFlag = true;
            end

        elseif RDEFlag
            % RDE: find closest ring radius
            [~, r1] = min(abs(R_RDE - cast(abs(y1(i)), 'like', T.R_RDE)));
            [~, r2] = min(abs(R_RDE - cast(abs(y2(i)), 'like', T.R_RDE)));

            err1(:) = (R_RDE(r1)^2 - abs(y1(i))^2);
            err2(:) = (R_RDE(r2)^2 - abs(y2(i))^2);
            for k = 1:NTaps
                w1V(k) = w1V(k) + mu_fxp * xv_i(k) * err1 * conj(y1(i));
                w1H(k) = w1H(k) + mu_fxp * xh_i(k) * err1 * conj(y1(i));
                w2V(k) = w2V(k) + mu_fxp * xv_i(k) * err2 * conj(y2(i));
                w2H(k) = w2H(k) + mu_fxp * xh_i(k) * err2 * conj(y2(i));
            end
        end

        % --- Reinitialisation for SingleSpike ---
        if i == N1 && SingleSpike
            w2H(:) = conj(w1V(end:-1:1, 1));
            w2V(:) = -conj(w1H(end:-1:1, 1));
        end
    end

    %% Collect output and remove transient
    y = [y1, y2];
    y = y(1+NOut:end, :);
end
