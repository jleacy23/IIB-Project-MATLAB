function X = fft_flp(x, inverse, po2Twiddle)
%FFT_FLP  Floating-point radix-2 FFT/IFFT with optional power-of-2 twiddles.
%
%   X = fft_flp(x)
%   X = fft_flp(x, inverse)
%   X = fft_flp(x, inverse, po2Twiddle)
%
%   Cooley-Tukey decimation-in-time radix-2 FFT in double precision.
%   Matches MATLAB's built-in fft / ifft (to working precision) when
%   po2Twiddle is false; when true, each twiddle-factor component
%   (real and imag) is snapped to the nearest signed power of two, so the
%   corresponding multiplications reduce to bit-shifts in a hardware
%   implementation.  Use this floating-point variant in algorithm-
%   development paths; the bit-true fixed-point version is fft.fft_fxp.
%
%   Inputs
%     x          - input signal, size [N x ...].  N must be a power of 2.
%                  The transform is taken along the first dimension; any
%                  trailing dimensions are processed independently
%                  (matches fft's column-wise behaviour on matrices and
%                  arrays).
%     inverse    - logical; true -> IFFT (conjugate twiddles, no
%                  scaling).  The forward FFT applies 1/2 scaling between
%                  each butterfly stage, spreading the full 1/N
%                  normalisation across the log2(N) stages.  Default false.
%     po2Twiddle - logical; true -> snap each twiddle component to the
%                  nearest signed power of two.  Default false.
%
%   Output
%     X          - transformed signal, same size as x.

    if nargin < 2 || isempty(inverse)
        inverse = false;
    end
    if nargin < 3 || isempty(po2Twiddle)
        po2Twiddle = false;
    end

    inSize = size(x);
    N = inSize(1);
    if N < 2 || bitand(N, N - 1) ~= 0
        error('fft_flp:lenNotPo2', ...
              'Input length along dim 1 must be a power of 2 (got %d).', N);
    end
    numStages = round(log2(N));

    M = prod(inSize(2:end));
    x = reshape(x, N, M);

    %% Bit-reverse permutation
    X = complex(zeros(N, M));
    for i = 0:N-1
        rev = bitrev_idx(i, numStages);
        X(rev + 1, :) = x(i + 1, :);
    end

    %% Cooley-Tukey butterfly stages
    for s = 1:numStages
        halfLen   = 2^(s - 1);
        fullLen   = 2^s;
        numGroups = N / fullLen;

        for k = 0:halfLen - 1
            theta = -2 * pi * k / fullLen;
            if inverse
                theta = -theta;
            end

            wr = cos(theta);
            wi = sin(theta);
            if po2Twiddle
                wr = roundPow2(wr);
                wi = roundPow2(wi);
            end
            W = complex(wr, wi);

            for g = 0:numGroups - 1
                idx_top = g * fullLen + k + 1;
                idx_bot = idx_top + halfLen;

                u = X(idx_top, :);
                t = W .* X(idx_bot, :);

                if inverse
                    % IFFT: no scaling
                    X(idx_top, :) = u + t;
                    X(idx_bot, :) = u - t;
                else
                    % Forward FFT: 1/2 inter-stage scaling so the full
                    % 1/N normalisation is spread across the log2(N)
                    % stages (mirrors the fixed-point fft_fxp path).
                    X(idx_top, :) = (u + t) / 2;
                    X(idx_bot, :) = (u - t) / 2;
                end
            end
        end
    end

    X = reshape(X, inSize);
end


%% ====================================================================
%  Local functions
%  ====================================================================

function rev = bitrev_idx(idx, nbits)
%BITREV_IDX  Bit-reverse a 0-based index with nbits bits.
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
