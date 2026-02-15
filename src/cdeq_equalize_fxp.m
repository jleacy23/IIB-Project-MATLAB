function Out = cdeq_equalize_fxp(In, D, L, CLambda, Rs, NPol, SpSIn, NFFT, po2Twiddle, T) %#codegen
%CDEQ_EQUALIZE_FXP  Fixed-point overlap-save frequency-domain CD compensation.
%
%   Out = cdeq_equalize_fxp(In, D, L, CLambda, Rs, NPol, SpSIn, NFFT,
%                            po2Twiddle, T)
%
%   Inputs (user-friendly units)
%     In         - input signal [samples x NPol] (fi or double)
%     D          - dispersion coefficient [ps/(nm*km)]
%     L          - fibre length [km]
%     CLambda    - central wavelength [nm]
%     Rs         - symbol rate [GBd]
%     NPol       - number of polarizations (1 or 2)
%     SpSIn      - samples per symbol
%     NFFT       - FFT block size (must be a power of 2)
%     po2Twiddle - logical; when true, FFT twiddle factors are rounded
%                  to the nearest signed power of 2 (shifts only)
%     T          - (optional) fixed-point types table from
%                  cdeq_equalize_fxp_types.  If omitted, uses 'fixed32'.
%
%   The types table T must supply prototype fi objects for:
%     T.x    - input / output signal
%     T.tw   - FFT twiddle factors (passed to fft_fxp)
%     T.hcd  - CD frequency-response coefficients
%     T.acc  - accumulator (FFT butterflies + freq-domain multiply)
%
%   Implementation notes for codegen:
%     - Uses fft_fxp (radix-2 Cooley-Tukey) for FFT/IFFT.
%     - The CD transfer function H_CD is pre-computed in double,
%       reordered from fftshift to natural FFT bin order via ifftshift,
%       then cast to fi — eliminating fftshift/ifftshift in the loop.
%     - Instead of an explicit overlap buffer, the input is prepended
%       with NOverlap zeros and a sliding window of stride stepLen
%       naturally captures the overlap-save behaviour.
%     - Each polarisation is processed independently to avoid 3-D fi
%       arrays, which simplifies the codegen datapath.
%     - All arithmetic uses SpecifyPrecision fimath so every product
%       and sum is truncated to the same WL/FL — no bit growth.

    %% Default types table
    if nargin < 10 || isempty(T)
        T = cdeq_equalize_fxp_types('fixed32');
    end

    c = 299792458;

    %% Unit conversion (double — used only for H_CD computation)
    D_si       = D * 1e-6;
    L_si       = L * 1e3;
    CLambda_si = CLambda * 1e-9;
    Rs_si      = Rs * 1e9;

    %% Compute overlap length (double helper, same as float version)
    NOverlap = cdeq_computeOverlap(D, L, CLambda, Rs, SpSIn, NFFT);
    stepLen  = NFFT - NOverlap;
    halfOv   = NOverlap / 2;          % NOverlap is guaranteed even

    %% CD frequency response — pre-computed in double, natural FFT order
    %  HCD is computed in fftshift (centred) order then reordered via
    %  ifftshift so that no fftshift/ifftshift is needed inside the loop.
    n  = (0:NFFT-1)' - NFFT/2;        % centred frequency indices
    fN = SpSIn * Rs_si / 2;           % Nyquist frequency

    HCD_centered = exp(-1i * pi * CLambda_si^2 * D_si * L_si / c * ...
                       (n * 2 * fN / NFFT).^2);

    HCD_nat = ifftshift(HCD_centered);       % natural FFT bin order
    HCD_fi  = cast(HCD_nat, 'like', T.hcd);  % cast to fixed-point

    %% Build FFT types sub-table for fft_fxp calls
    Tfft.x   = T.x;
    Tfft.tw  = T.tw;
    Tfft.acc = T.acc;

    %% Input extension (same logic as cdeq_equalize)
    NIn    = size(In, 1);
    AuxLen = NIn / stepLen;

    if AuxLen ~= ceil(AuxLen)
        NExtra = ceil(AuxLen) * stepLen - NIn;
    else
        NExtra = NOverlap;
    end

    halfEx = NExtra / 2;

    InPad = cast([In(end - halfEx + 1 : end, :); ...
                  In; ...
                  In(1 : halfEx, :)], 'like', T.x);

    NPadded = size(InPad, 1);
    nBlocks = NPadded / stepLen;

    %% Prepend NOverlap zeros so the overlap-save sliding window works
    %  without an explicit overlap buffer.  Block i starts at index
    %  (i-1)*stepLen + 1 in InPad2 and spans NFFT samples.
    InPad2 = [complex(zeros(NOverlap, NPol, 'like', T.x)); InPad];

    %% Pre-allocate output (accumulator type)
    OutPad = complex(zeros(NPadded, NPol, 'like', T.acc));

    %% ================================================================
    %  Per-polarisation overlap-save processing
    %  ================================================================
    for pol = 1:NPol
        for i = 1:nBlocks
            wStart = (i - 1) * stepLen + 1;

            % --- Extract NFFT-length block (overlap + new data) ------
            InB = InPad2(wStart : wStart + NFFT - 1, pol);

            % --- Forward FFT ----------------------------------------
            X = fft_fxp(InB, false, po2Twiddle, Tfft);

            % --- Frequency-domain CD compensation -------------------
            Y = complex(zeros(NFFT, 1, 'like', T.acc));
            for k = 1:NFFT
                Y(k) = X(k) * HCD_fi(k);
            end

            % --- Inverse FFT ----------------------------------------
            outFDE = fft_fxp(Y, true, po2Twiddle, Tfft);

            % --- Keep only the valid (non-aliased) samples -----------
            oStart = (i - 1) * stepLen + 1;
            for j = 1:stepLen
                OutPad(oStart + j - 1, pol) = outFDE(halfOv + j);
            end
        end
    end

    %% Remove extension samples
    DInit = 1 + (NExtra + NOverlap) / 2;
    DFin  = (NExtra - NOverlap) / 2;
    Out   = OutPad(DInit : end - DFin, :);
end
