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
%     inverse    - logical scalar; true -> IFFT (conjugate twiddles + 1/N)
%     po2Twiddle - logical scalar; true -> round each twiddle-factor
%                  component (real, imag) to the nearest signed power of 2
%                  so that every twiddle multiplication becomes a bit-shift
%     T          - (optional) fixed-point types table from fft_fxp_types.
%                  If omitted, defaults to fft_fxp_types('fixed32').
%
%   Output
%     X          - complex column vector [N x 1]
%
%   The accumulator type T.acc is used for all intermediate butterfly
%   values and the output.  The fimath attached to T.acc controls
%   product / sum word lengths.  No normalisation is applied between
%   butterfly stages; for the IFFT the 1/N scaling is applied once at
%   the output.
%
%   Types table fields
%     T.x    - input signal prototype (initial cast)
%     T.tw   - twiddle-factor prototype
%     T.acc  - accumulator / output prototype (butterfly sums & products)

    %% Defaults
    if nargin < 4 || isempty(T)
        T = fft_fxp_types('fixed32');
    end

    N         = size(x, 1);
    numStages = round(log2(double(N)));

    %% ----------------------------------------------------------------
    %  Bit-reverse permutation: x -> X  (cast into accumulator type)
    %  ----------------------------------------------------------------
    X = complex(zeros(N, 1, 'like', T.acc));
    for i = 0:N-1
        rev    = bitrev_idx(i, numStages);
        X(rev + 1) = cast(x(i + 1), 'like', T.acc);
    end

    %% ----------------------------------------------------------------
    %  Cooley-Tukey butterfly stages
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

            % --- twiddle components ---
            wr = cos(theta);
            wi = sin(theta);
            if po2Twiddle
                wr = roundPow2(wr);
                wi = roundPow2(wi);
            end

            W = complex(cast(wr, 'like', T.tw), ...
                        cast(wi, 'like', T.tw));

            % --- butterfly across groups ---
            for g = 0:numGroups - 1
                idx_top = g * fullLen + k + 1;
                idx_bot = idx_top + halfLen;

                u = X(idx_top);
                t = W * X(idx_bot);

                X(idx_top) = u + t;
                X(idx_bot) = u - t;
            end
        end
    end

    %% ----------------------------------------------------------------
    %  IFFT: scale output by 1/N
    %  ----------------------------------------------------------------
    if inverse
        M = cast(log2(double(N)), 'int32');
        for i = 1:N
            X(i) = bitshift(X(i), -M);  % equivalent to X(i) / N with arithmetic shift
        end
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
