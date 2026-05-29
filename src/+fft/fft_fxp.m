function X = fft_fxp(x, inverse, po2Twiddle, T) %#codegen
%FFT_FXP  Fixed-point radix-2 FFT/IFFT for code generation.
%
%   X = fft_fxp(x, inverse, po2Twiddle, T)
%
%   Performs an in-place Cooley-Tukey decimation-in-time radix-2 FFT (or
%   IFFT when inverse == true) on a fixed-size complex column vector.
%
%   Inputs
%     x          - complex column vector [N x 1], N must be a power of 2
%     inverse    - logical scalar; true -> IFFT (conjugate twiddles)
%     po2Twiddle - logical scalar; true -> round each twiddle-factor
%                  component (real, imag) to the nearest signed power of 2
%                  so that every twiddle multiplication becomes a bit-shift
%     T          - (optional) fixed-point types table from fft_fxp_types.
%                  If omitted, defaults to fft_fxp_types('fixed32').
%
%   Output
%     X          - complex column vector [N x 1], type T.acc
%
%   Normalisation / precision strategy
%     The radix-2 butterflies run UNSCALED in a WIDE internal accumulator
%     (wideAccLike), whose word length is a fixed design constant with
%     enough integer headroom for the ~N-fold magnitude growth of an
%     unscaled transform and generous fractional bits.  The transform is
%     normalised and quantised to the output type T.acc ONLY at the end:
%       * forward FFT: divide by N (single exact power-of-2 shift) then
%         cast to T.acc, so the output is 1/N-normalised (== fft(x)/N);
%       * inverse FFT: no scaling, cast to T.acc (== ifft(X)*N).
%     The forward<->inverse round trip is therefore the identity.  Keeping
%     the accumulator wide avoids the per-stage truncation noise (and DC
%     bias) of a scaled FFT, so the only quantisation is the T.acc output
%     register precision and the T.tw twiddle-ROM precision.  Float configs
%     (T.acc double/single) run natively in that type.
%
%   Types table fields
%     T.x    - input signal prototype (initial cast)
%     T.tw   - twiddle-factor (ROM) prototype
%     T.acc  - output prototype (the desired quantised precision)

    %% Defaults
    if nargin < 4 || isempty(T)
        T = fft.fft_fxp_types('fixed32');
    end

    N         = size(x, 1);
    numStages = round(log2(double(N)));

    %% Wide internal accumulator (fixed-point) or native float type
    accWide = wideAccLike(T.acc);
    accIsFi = isfi(accWide);

    %% ----------------------------------------------------------------
    %  Bit-reverse permutation: x -> Xw  (cast into the wide accumulator)
    %  ----------------------------------------------------------------
    Xw = complex(zeros(N, 1, 'like', accWide));
    for i = 0:N-1
        rev    = bitrev_idx(i, numStages);
        Xw(rev + 1) = cast(x(i + 1), 'like', accWide);
    end

    %% ----------------------------------------------------------------
    %  Cooley-Tukey butterfly stages (UNSCALED, wide accumulator)
    %  ----------------------------------------------------------------
    for s = 1:numStages
        halfLen   = 2^(s - 1);           % half-butterfly span
        fullLen   = 2^s;                  % full span
        numGroups = N / fullLen;          % groups this stage

        for k = 0:halfLen - 1
            % --- twiddle angle ---
            theta = -2 * pi * double(k) / double(fullLen);
            if inverse
                theta = -theta;           % conjugate for IFFT
            end

            % --- twiddle components (quantised to the T.tw ROM type) ---
            wr = cos(theta);
            wi = sin(theta);
            if po2Twiddle
                wr = roundPow2(wr);
                wi = roundPow2(wi);
            end
            W  = complex(cast(wr, 'like', T.tw), ...
                         cast(wi, 'like', T.tw));
            % Promote the (T.tw-precision) twiddle into the wide accumulator
            % so the MAC accumulates at full width (operands share fimath).
            Wc = cast(W, 'like', accWide);

            % --- butterfly across groups ---
            for g = 0:numGroups - 1
                idx_top = g * fullLen + k + 1;
                idx_bot = idx_top + halfLen;

                u = Xw(idx_top);
                t = Wc * Xw(idx_bot);

                Xw(idx_top) = u + t;       % unscaled
                Xw(idx_bot) = u - t;
            end
        end
    end

    %% ----------------------------------------------------------------
    %  Normalise (1/N on forward only) and quantise to the output type
    %  ----------------------------------------------------------------
    X = complex(zeros(N, 1, 'like', T.acc));
    if inverse
        for i = 1:N
            X(i) = cast(Xw(i), 'like', T.acc);
        end
    elseif accIsFi
        M = cast(numStages, 'int32');
        for i = 1:N
            % 1/N via a single arithmetic right shift in the wide
            % accumulator (exact, power of two), then quantise to T.acc.
            X(i) = cast(bitshift(Xw(i), -M), 'like', T.acc);
        end
    else
        for i = 1:N
            X(i) = cast(Xw(i) / N, 'like', T.acc);
        end
    end
end

%% ====================================================================
%  Local functions
%  ====================================================================

function w = wideAccLike(proto)
%WIDEACCLIKE  Wide fixed-point accumulator prototype for the FFT internals.
%   Word length is a FIXED design constant (not the swept output precision):
%   ~20 integer bits hold the unscaled ~N-fold magnitude growth and 28
%   fractional bits keep the butterfly MACs near-lossless.  For
%   floating-point output types the FFT simply runs in that type.
    if isfi(proto)
        WL = 48; FL = 28;
        F  = fimath( ...
            'RoundingMethod',       'Floor', ...
            'OverflowAction',       'Wrap',  ...
            'ProductMode',          'SpecifyPrecision', ...
            'ProductWordLength',     WL, ...
            'ProductFractionLength', FL, ...
            'SumMode',              'SpecifyPrecision', ...
            'SumWordLength',         WL, ...
            'SumFractionLength',     FL);
        w = fi([], 1, WL, FL, F);
    else
        w = proto;
    end
end

function rev = bitrev_idx(idx, nbits)
%BITREV_IDX  Bit-reverse a 0-based index with nbits bits.
%   Uses pure scalar arithmetic (mod / floor) so that no integer-only
%   bitwise intrinsics are needed — fully codegen-compatible.
    rev = 0;
    val = idx;
    for b = 1:nbits          %#ok<FXUP>
        rev = rev * 2 + mod(val, 2);
        val = floor(val / 2);
    end
end

function y = roundPow2(x)
%ROUNDPOW2  Snap x to the nearest signed power of two.
%   Values with |x| < 2^{-16} are mapped to zero (catches floating-point
%   near-zeros such as cos(pi/2) ~ 6e-17).
    if abs(x) < pow2(-16)
        y = 0;
    else
        y = sign(x) * pow2(round(log2(abs(x))));
    end
end
