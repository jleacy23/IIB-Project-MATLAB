function Out = recovery_godard_fxp(In, NSymb, N, beta, G, po2Twiddle, T) %#codegen
%recovery_godard_fxp  Fixed-point feedforward Modified Godard clock recovery.
%
%   Out = recovery_godard_fxp(In, NSymb, N, beta, G, po2Twiddle, T)
%
%   See clk_recovery.recovery_godard for the floating-point reference.
%   Both the timing estimate and the timing correction are performed in
%   the frequency domain on the same FFT block of N samples.
%
%   Per block:
%       τ̂_T,w = atan2(Σ Im, Σ Re) / (2π)        (wrapped, in (-0.5, 0.5]·T)
%       τ̂_T   = unwrap(τ̂_T,w  vs  previous block)
%       τ̂_s   = τ̂_T · η                          (cumulative input-sample delay)
%       R_corr(k) = R(k) · exp(-j 2π k_idx τ̂_s / N)
%       Out_blk   = IFFT( R_corr )
%
%   Overlap-save: input is zero-padded with G samples each side; blocks
%   step by M = N − 2G input samples; G samples are discarded at each
%   end of every IFFT block.  G must be >= the maximum expected |τ̂_s|
%   (cumulative SFO drift across the record, in input samples).
%
%   Inputs
%     In         - input signal at 2 Sa/symbol (column vector, one pol)
%     NSymb      - number of transmitted symbols (output limited to NSymb*2)
%     N          - FFT block size (power of 2)
%     beta       - pulse-shaping roll-off factor (0 < beta <= 1)
%     G          - overlap-save guard length per block edge (< N/2)
%     po2Twiddle - logical; round FFT twiddles to nearest signed power of 2
%     T          - (optional) fixed-point types table from
%                  recovery_godard_fxp_types.  Defaults to 'fixed32'.
%
%   Output
%     Out        - clock-recovered signal (column vector, 2 Sa/symbol)
%
%   Codegen notes
%     - Uses fft.fft_fxp for both forward and inverse transforms.
%     - atan2, round, cos and sin are evaluated in double via the
%       standard fi-to-double escape — only the results are cast back
%       to fi.
%     - All bulk arithmetic (FFT butterflies, MG sum, freq-domain
%       multiply) is in fi with SpecifyPrecision fimath.

    if nargin < 7 || isempty(T)
        T = clk_recovery.recovery_godard_fxp_types('fixed32');
    end

    eta = 2;

    % MG bin shift and excess-bandwidth summation bounds
    shift = int32(round((1 - 1/eta) * double(N)));
    kLo   = int32(round((1 - beta) / (2*eta) * double(N)) + 1);
    kHi   = int32(round((1 + beta) / (2*eta) * double(N)));

    % FFT sub-types table
    Tfft.x   = T.x;
    Tfft.tw  = T.tw;
    Tfft.acc = T.acc;

    NN    = int32(N);
    GG    = int32(G);
    MM    = NN - int32(2) * GG;
    halfN = int32(N / 2);
    invN  = 1 / double(N);
    LIn   = int32(length(In));

    % Padded input length and block count
    LPad    = LIn + int32(2) * GG;
    nBlocks = int32(floor(double(LPad - NN) / double(MM))) + int32(1);

    % Build padded input buffer (zeros at both ends)
    InPad = complex(zeros(LPad, 1, 'like', T.x));
    for j = int32(1):LIn
        InPad(GG + j) = cast(In(j), 'like', T.x);
    end

    maxOut    = int32(NSymb * 2);
    Out       = complex(zeros(maxOut, 1, 'like', T.x));
    tauT_prev = 0;                           % running unwrapped estimate (T units)

    for b = int32(1):nBlocks
        inStart = (b - int32(1)) * MM + int32(1);
        if inStart + NN - int32(1) > LPad
            break;
        end

        % --- Extract block ---
        InB = complex(zeros(NN, 1, 'like', T.x));
        for j = int32(1):NN
            InB(j) = InPad(inStart + j - int32(1));
        end

        % --- Forward FFT ---
        R = fft.fft_fxp(InB, false, po2Twiddle, Tfft);

        % --- MG sum over excess-BW bins ---
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

        % --- Feedforward MG estimate + unwrap ---
        tauT_w    = atan2(double(sum_im), double(sum_re)) / (2*pi);
        tauT      = tauT_w + round(tauT_prev - tauT_w);
        tauT_prev = tauT;
        tauSamp   = tauT * eta;
        dPhi      = -2 * pi * tauSamp * invN;        % phase slope per bin

        % --- Phase-ramp correction ---
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

        % --- Overlap-save: keep middle M samples ---
        outStart = (b - int32(1)) * MM + int32(1);
        for j = int32(1):MM
            outIdx = outStart + j - int32(1);
            if outIdx > maxOut
                break;
            end
            if outIdx > LIn
                break;
            end
            srcIdx = GG + j;
            Out(outIdx) = cast(outBlk(srcIdx), 'like', T.x);
        end
    end
end
