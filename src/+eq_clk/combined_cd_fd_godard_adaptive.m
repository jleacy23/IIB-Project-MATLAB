function [y, cfoBinsApplied] = combined_cd_fd_godard_adaptive(In, SpS, ...
        NFFT, NOverlap, D, L, CLambda, Rs, Rolloff, ki, kp, NSymb, ...
        AdaptOpts, cfoEnable, po2Twiddle)
%COMBINED_CD_FD_GODARD_ADAPTIVE  Frequency-domain CD + RRC matched filter
%   (overlap-save) with optional one-shot coarse CFO correction and
%   Godard timing recovery sharing the same FFT and overlap, followed
%   by a butterfly CMA equaliser.
%
%   y = combined_cd_fd_godard_adaptive(In, SpS, NFFT, NOverlap, D, L, ...
%           CLambda, Rs, Rolloff, ki, kp, NSymb, AdaptOpts)
%   y = combined_cd_fd_godard_adaptive(..., AdaptOpts, cfoEnable)
%   y = combined_cd_fd_godard_adaptive(..., cfoEnable, po2Twiddle)
%
%   Inputs
%     In          - input signal [samples x 2]
%     SpS         - samples per symbol
%     NFFT        - FFT block size (shared between CD/MF and Godard)
%     NOverlap    - overlap-save overlap length (even, < NFFT)
%     D, L, CLambda, Rs - dispersion / signal parameters
%     Rolloff     - RRC roll-off (also defines the Godard excess-band bins)
%     ki, kp      - Godard PI loop-filter gains (units of tau in samples)
%     NSymb       - number of transmitted symbols
%     AdaptOpts   - adaptive equaliser settings struct (NTaps, Mu, ...) -
%                   see eq_clk.apply_adaptive_eq.
%     cfoEnable   - (optional) when truthy (nonzero), apply a one-shot
%                   coarse CFO correction (eq_clk.coarse_cfo_fd):
%                   centroid of the first NFFT samples' FFT drives a
%                   continuous-valued time-domain phasor across the
%                   whole input.  Default false.
%     po2Twiddle  - (optional) when true, the per-block FFT/IFFT use
%                   fft.fft_flp with twiddle factors snapped to the
%                   nearest signed power of two (for hardware shift-only
%                   multiplications).  Default false.
%
%   The Godard Modified-Godard timing metric is evaluated on the
%   already-corrected spectrum (CD + MF + current phase ramp).  Its imag
%   part drives a per-block PI loop filter whose accumulated tau is split
%   NCO-style into an integer part (folded into the per-block read pointer)
%   and a bounded fractional residual (applied as the next-block phase
%   ramp), so the ramp slope never grows past +/-half a sample.  The
%   processing is a streaming zero-padded overlap-save (no cyclic frame
%   extension), so it no longer requires the total drift to be an integer
%   number of samples.

    if nargin < 14 || isempty(cfoEnable)
        cfoEnable = false;
    end
    if nargin < 15 || isempty(po2Twiddle)
        po2Twiddle = false;
    end

    NPol = size(In, 2);

    %% One-shot coarse CFO correction (time-domain phasor) ---------
    cfoBinsApplied = 0;
    if cfoEnable
        [In, cfoBinsApplied] = eq_clk.coarse_cfo_fd(In, NFFT);
    end

    %% Static frequency masks (natural FFT order) --------------------
    HCDshift = eq_clk.cd_fd_response(D, L, CLambda, Rs, SpS, NFFT);
    HMFshift = eq_clk.rrc_fd_response(Rolloff, NFFT, SpS);
    Hstatic  = ifftshift(HCDshift .* HMFshift);   % natural FFT order

    %% Godard band parameters (natural order, recovery_godard convention)
    eta   = SpS;
    shift = round((1 - 1/eta) * NFFT);
    kLo   = round((1 - Rolloff) / (2*eta) * NFFT) + 1;
    kHi   = round((1 + Rolloff) / (2*eta) * NFFT);
    k_idx = [0:NFFT/2-1, -NFFT/2:-1].';

    %% Zero-padded streaming overlap-save framing (no cyclic extension) --
    %  Mirrors the working clk_recovery.recovery_godard: a guard of halfOv
    %  zeros on each end, blocks read at stride stepLen, halfOv discarded
    %  per edge.  The accumulated timing is split NCO-style inside the loop
    %  into an integer part (folded into the per-block read pointer) and a
    %  bounded fractional residual (the phase ramp).
    NIn     = size(In, 1);
    stepLen = NFFT - NOverlap;          % M: valid output samples per block
    halfOv  = NOverlap / 2;             % G: overlap-save guard per edge

    InPad   = [zeros(halfOv, NPol); In; zeros(halfOv, NPol)];
    LPad    = size(InPad, 1);
    nBlk    = floor((LPad - NFFT) / stepLen) + 1;

    %% Overlap-save loop with embedded Godard PI ---------------------
    Out = zeros(nBlk * stepLen, NPol);

    LF_I    = 0;
    tauSamp = 0;

    for i = 1:nBlk
        % NCO-style integer/fractional split of the timing estimate: the
        % integer part folds into the read pointer, only the bounded
        % fractional residual mu drives the phase ramp.
        dInt = round(tauSamp);
        mu   = tauSamp - dInt;           % fractional residual in [-0.5, 0.5)
        ramp = exp(-1j * 2*pi * k_idx * mu / NFFT);

        % Read window with the integer timing folded into the read pointer;
        % indices outside the padded input are zero-filled.
        rdStart = (i - 1) * stepLen + 1 - dInt;
        idx     = (rdStart : rdStart + NFFT - 1).';
        valid   = idx >= 1 & idx <= LPad;
        InB     = zeros(NFFT, 1, NPol);
        for p = 1:NPol
            InB(valid, 1, p) = InPad(idx(valid), p);
        end

        % Natural-order FFT (matches recovery_godard convention)
        R = fft.fft_flp(InB, false, po2Twiddle);

        % CD + matched filter + current fractional phase ramp
        R_corr = (R .* Hstatic) .* ramp;

        % Godard metric on the corrected spectrum (per-pol then summed)
        S = 0;
        for p = 1:NPol
            Rp = R_corr(:, 1, p);
            S  = S + sum(Rp(kLo:kHi) .* conj(Rp(kLo+shift:kHi+shift)));
        end
        e = imag(S);

        % PI update for the next block
        LF_I    = LF_I + ki * e;
        tauSamp = kp * e + LF_I;

        % IFFT and overlap-save save
        OutFDE = fft.fft_flp(R_corr, true, po2Twiddle);
        oStart = (i - 1) * stepLen + 1;
        Out(oStart : oStart + stepLen - 1, :) = ...
            reshape(OutFDE(halfOv+1 : halfOv+stepLen, 1, :), stepLen, NPol);
    end

    %% Trim to the original input span (output aligns from sample 1) -----
    if size(Out, 1) > NIn
        z = Out(1:NIn, :);
    else
        z = Out;
    end

    %% Adaptive butterfly CMA ---------------------------------------
    if NPol == 1
        z = [z, z];
        singlePol = true;
    else
        singlePol = false;
    end
    y = eq_clk.apply_adaptive_eq(z, SpS, AdaptOpts);
    if singlePol
        y = y(:, 1);
    end

    if ~isempty(NSymb) && NSymb > 0
        NOutMax = NSymb;
        if size(y, 1) > NOutMax
            y = y(1:NOutMax, :);
        end
    end
end
