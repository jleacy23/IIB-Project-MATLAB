function Out = recovery_godard_fxp(In, NSymb, N, beta, po2Twiddle, T) %#codegen
%recovery_godard_fxp  Fixed-point feedforward Modified Godard clock recovery.
%
%   Out = recovery_godard_fxp(In, NSymb, N, beta, po2Twiddle, T)
%
%   See clk_recovery.recovery_godard for the floating-point reference.
%   Both the timing estimate and the timing correction are performed in
%   the frequency domain on the same FFT block of N samples.
%
%   Per block:
%       τ̂_T   = atan2(Σ Im, Σ Re) / (2π)        (Josten 2017 eq. 5,
%                                                  arg/(2π) feedforward)
%       τ̂_s   = τ̂_T · η                         (input-sample delay)
%       R_corr(k) = R(k) · exp(-j 2π k_idx τ̂_s / N)
%       Out_blk   = IFFT( R_corr )
%
%   Inputs
%     In         - input signal at 2 Sa/symbol (column vector, one pol)
%     NSymb      - number of transmitted symbols (output limited to NSymb*2)
%     N          - FFT block size (power of 2)
%     beta       - pulse-shaping roll-off factor (0 < beta <= 1)
%     po2Twiddle - logical; round FFT twiddles to nearest signed power of 2
%     T          - (optional) fixed-point types table from
%                  recovery_godard_fxp_types.  Defaults to 'fixed32'.
%
%   Output
%     Out        - clock-recovered signal (column vector, 2 Sa/symbol)
%
%   Codegen notes
%     - Uses fft.fft_fxp for both forward and inverse transforms.
%     - atan2, cos and sin are evaluated in double via the standard
%       fi-to-double escape (matches the pattern in pilots_only_fxp /
%       fft_search_fxp).  Only the result is cast back to fi.
%     - All bulk arithmetic (FFT butterflies, MG sum, freq-domain
%       multiply) is done in fi with SpecifyPrecision fimath.

    if nargin < 6 || isempty(T)
        T = clk_recovery.recovery_godard_fxp_types('fixed32');
    end

    eta = 2;                                            % input SpS

    % MG bin shift and excess-bandwidth summation bounds
    shift = int32(round((1 - 1/eta) * double(N)));
    kLo   = int32(round((1 - beta) / (2*eta) * double(N)) + 1);
    kHi   = int32(round((1 + beta) / (2*eta) * double(N)));

    % FFT sub-types table
    Tfft.x   = T.x;
    Tfft.tw  = T.tw;
    Tfft.acc = T.acc;

    NN      = int32(N);
    halfN   = int32(N / 2);
    invN    = 1 / double(N);
    LIn     = int32(length(In));
    nBlocks = int32(floor(double(LIn) / double(N)));

    maxOut = int32(NSymb * 2);
    Out    = complex(zeros(maxOut, 1, 'like', T.x));

    for iBlk = int32(1):nBlocks
        blkStart = (iBlk - int32(1)) * NN + int32(1);

        % --- Extract block, cast to FFT input type ---
        InB = complex(zeros(NN, 1, 'like', T.x));
        for j = int32(1):NN
            InB(j) = cast(In(blkStart + j - int32(1)), 'like', T.x);
        end

        % --- Forward FFT ---
        R = fft.fft_fxp(InB, false, po2Twiddle, Tfft);

        % --- MG sum: Σ R(k) · conj(R(k+shift)) over excess-BW bins ---
        sum_re = cast(0, 'like', T.acc);
        sum_im = cast(0, 'like', T.acc);
        for k = kLo:kHi
            a_re = real(R(k));
            a_im = imag(R(k));
            c_re = real(R(k + shift));
            c_im = imag(R(k + shift));
            sum_re(:) = sum_re + (a_re * c_re + a_im * c_im);
            sum_im(:) = sum_im + (a_im * c_re - a_re * c_im);
        end

        % --- Feedforward MG estimate (units of T) → input-sample delay ---
        tauT    = atan2(double(sum_im), double(sum_re)) / (2*pi);
        tauSamp = tauT * eta;
        dPhi    = -2 * pi * tauSamp * invN;          % phase slope per bin

        % --- Frequency-domain phase-ramp correction ---
        R_corr = complex(zeros(NN, 1, 'like', T.acc));
        for k = int32(1):NN
            if k <= halfN
                kBin = double(k) - 1;
            else
                kBin = double(k) - 1 - double(N);
            end
            phi  = dPhi * kBin;
            cc   = cast(cos(phi), 'like', T.tw);
            ss   = cast(sin(phi), 'like', T.tw);
            r_re = real(R(k));
            r_im = imag(R(k));
            R_corr(k) = complex(r_re * cc - r_im * ss, ...
                                r_re * ss + r_im * cc);
        end

        % --- Inverse FFT ---
        outBlk = fft.fft_fxp(R_corr, true, po2Twiddle, Tfft);

        % --- Write to output buffer (capped at maxOut) ---
        for j = int32(1):NN
            outIdx = blkStart + j - int32(1);
            if outIdx > maxOut
                break;
            end
            Out(outIdx) = cast(outBlk(j), 'like', T.x);
        end
    end

    % --- Pass-through tail samples that don't fill a complete block ---
    tailStart = nBlocks * NN + int32(1);
    for j = tailStart:LIn
        if j > maxOut
            break;
        end
        Out(j) = cast(In(j), 'like', T.x);
    end
end
