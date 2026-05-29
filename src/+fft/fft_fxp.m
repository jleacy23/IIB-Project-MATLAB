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
%     The FFT runs at a SINGLE precision, T.acc: the input is cast to T.acc
%     immediately, the radix-2 butterflies run in T.acc, and the output is
%     T.acc.  (In the uniform static configuration T.x == T.acc.)
%       * forward FFT: each radix-2 stage divides its butterfly outputs by
%         2 and writes them back to the accumulator.  Over the log2(N)
%         stages this applies the full 1/N normalisation (== fft(x)/N) while
%         keeping the magnitude bounded (no N-fold growth);
%       * inverse FFT: NO per-stage division (== ifft(X)*N), so the
%         forward<->inverse round trip is unit gain.
%     The only quantisation is the T.acc datapath precision and the T.tw
%     twiddle-ROM precision.  Float configs run natively in that type.
%
%   Types table fields
%     T.x    - upstream I/O prototype (uniform static config: == T.acc)
%     T.tw   - twiddle-factor (ROM) prototype
%     T.acc  - single FFT datapath precision (input cast, butterflies, output)

    %% Defaults
    if nargin < 4 || isempty(T)
        T = fft.fft_fxp_types('fixed32');
    end

    N         = size(x, 1);
    numStages = round(log2(double(N)));

    %% Internal accumulator: the caller-supplied T.acc type (set in the
    %  types table / test).  Float configs run natively in that type.
    accWide = T.acc;
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
    %  Cooley-Tukey butterfly stages.  Forward: each stage divides by 2
    %  (block-floating-point style).  Inverse: unscaled.
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

                top = u + t;
                bot = u - t;

                % Forward pass: divide each stage by 2 (cumulative 1/N over
                % the log2(N) stages) and write back to the accumulator.
                % Inverse pass: no scaling, so the round trip is unit gain.
                if ~inverse
                    if accIsFi
                        top = bitshift(top, -1);
                        bot = bitshift(bot, -1);
                    else
                        top = top / 2;
                        bot = bot / 2;
                    end
                end

                Xw(idx_top) = top;
                Xw(idx_bot) = bot;
            end
        end
    end

    %% ----------------------------------------------------------------
    %  Quantise to the output type T.acc (the single FFT precision).
    %  Forward scaling (1/N) has already been applied per-stage above;
    %  inverse is unscaled.
    %  ----------------------------------------------------------------
    X = complex(zeros(N, 1, 'like', T.acc));
    for i = 1:N
        X(i) = cast(Xw(i), 'like', T.acc);
    end
end

%% ====================================================================
%  Local functions
%  ====================================================================

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
