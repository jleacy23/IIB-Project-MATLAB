function [y, cfoBinsApplied] = combined_cd_fd_godard_adaptive_fxp(...
        In, SpS, NFFT, NOverlap, D, L, CLambda, Rs, Rolloff, ki, kp, ...
        NSymb, AdaptOpts, cfoEnable, po2Twiddle, T) %#codegen
%COMBINED_CD_FD_GODARD_ADAPTIVE_FXP  Fixed-point CD-FD + RRC matched filter
%   (overlap-save) with optional coarse CFO correction and embedded
%   Modified-Godard timing recovery sharing the same FFT and overlap,
%   followed by butterfly CMA.
%
%   [y, cfoBinsApplied] = combined_cd_fd_godard_adaptive_fxp(In, SpS, ...
%       NFFT, NOverlap, D, L, CLambda, Rs, Rolloff, ki, kp, NSymb, ...
%       AdaptOpts, cfoEnable, po2Twiddle, T)
%
%   The composite types table T has the following sub-tables (see
%   eq_clk.combined_cd_fd_godard_adaptive_fxp_types):
%
%     T.Static  - frequency-domain CD + matched filter
%     T.Godard  - embedded Godard timing recovery (metric, loop filter,
%                 phase-ramp twiddle)
%     T.AdaptEq - adaptive butterfly equaliser
%
%   The optional coarse CFO correction is applied in floating point
%   (cast fi -> double, eq_clk.coarse_cfo_fd in float, cast back to
%   T.Static.x).
%
%   The Godard PI loop runs in-line with the overlap-save iteration: each
%   block's corrected spectrum (CD + MF + current phase ramp) feeds the
%   Modified-Godard timing metric S; imag(S) drives the per-block PI loop
%   whose accumulated tau is applied as the next block's phase ramp.

    if nargin < 14 || isempty(cfoEnable)
        cfoEnable = false;
    end
    if nargin < 15 || isempty(po2Twiddle)
        po2Twiddle = false;
    end
    if nargin < 16 || isempty(T)
        T = eq_clk.combined_cd_fd_godard_adaptive_fxp_types('fixed32');
    end

    NPol = size(In, 2);

    %% ================================================================
    %  1. Coarse CFO correction (floating point)
    %  ================================================================
    cfoBinsApplied = 0;
    if cfoEnable
        InDouble = double(In);
        [InCFO_d, cfoBinsApplied] = eq_clk.coarse_cfo_fd(InDouble, NFFT);
        InS = cast(InCFO_d, 'like', T.Static.x);
    else
        InS = cast(In, 'like', T.Static.x);
    end

    %% ================================================================
    %  2. Static frequency masks (natural FFT order)
    %  ================================================================
    HCDshift = eq_clk.cd_fd_response(D, L, CLambda, Rs, SpS, NFFT);
    HMFshift = eq_clk.rrc_fd_response(Rolloff, NFFT, SpS);
    Hstatic  = ifftshift(HCDshift .* HMFshift);
    Hstatic_fi = cast(Hstatic, 'like', T.Static.hcd);

    Tfft.x   = T.Static.x;
    Tfft.tw  = T.Static.tw;
    Tfft.acc = T.Static.acc;

    %% Godard band parameters (natural order, recovery_godard convention)
    eta   = SpS;
    shift = round((1 - 1/eta) * NFFT);
    kLo   = round((1 - Rolloff) / (2*eta) * NFFT) + 1;
    kHi   = round((1 + Rolloff) / (2*eta) * NFFT);
    k_idx = [0:NFFT/2-1, -NFFT/2:-1].';

    %% Input cyclic extension to integer block count
    NIn     = size(InS, 1);
    stepLen = NFFT - NOverlap;
    AuxLen  = NIn / stepLen;
    if AuxLen ~= ceil(AuxLen)
        NExtra = ceil(AuxLen) * stepLen - NIn;
    else
        NExtra = NOverlap;
    end
    halfEx = NExtra / 2;
    halfOv = NOverlap / 2;

    InPad = cast([InS(end - halfEx + 1 : end, :); ...
                  InS; ...
                  InS(1 : halfEx, :)], 'like', T.Static.x);

    NPadded = size(InPad, 1);
    InPad2  = [complex(zeros(NOverlap, NPol, 'like', T.Static.x)); InPad];
    nBlocks = NPadded / stepLen;

    OutPad = complex(zeros(NPadded, NPol, 'like', T.Static.acc));

    %% ================================================================
    %  3. Overlap-save loop with embedded Godard PI
    %  ================================================================
    LF_I    = cast(0, 'like', T.Godard.lf);
    tauSamp = cast(0, 'like', T.Godard.lf);

    ki_fi = cast(ki, 'like', T.Godard.lf);
    kp_fi = cast(kp, 'like', T.Godard.lf);

    % FFT buffers and the corrected spectrum (kept across pols within a
    % block so the Godard metric can sum over polarisations).
    R_corr_all = complex(zeros(NFFT, NPol, 'like', T.Static.acc));

    for i = 1:nBlocks
        % --- Build the per-block phase ramp from current tauSamp -------
        %  Computed in double, then cast to T.Godard.tw (unit-magnitude
        %  complex twiddle).  Matches the fft_search_fxp pattern of
        %  computing exp(.) in double and casting back.
        tau_d  = double(tauSamp);
        ramp_d = exp(-1j * 2*pi * k_idx * tau_d / NFFT);
        ramp   = cast(ramp_d, 'like', T.Godard.tw);

        % --- Per-polarisation FFT, CD/MF mask, phase ramp --------------
        for pol = 1:NPol
            wStart = (i - 1) * stepLen + 1;
            InB    = InPad2(wStart : wStart + NFFT - 1, pol);

            X = fft.fft_fxp(InB, false, po2Twiddle, Tfft);

            % Apply CD + matched filter + current phase ramp
            for k = 1:NFFT
                Rfilt   = X(k) * Hstatic_fi(k);
                R_corr_all(k, pol) = cast(Rfilt * cast(ramp(k), 'like', T.Static.acc), ...
                                          'like', T.Static.acc);
            end

            % --- Inverse FFT and overlap-save save ---------------------
            outFDE = fft.fft_fxp(R_corr_all(:, pol), true, po2Twiddle, Tfft);

            oStart = (i - 1) * stepLen + 1;
            for j = 1:stepLen
                OutPad(oStart + j - 1, pol) = outFDE(halfOv + j);
            end
        end

        % --- Godard metric on the corrected spectrum (per-pol sum) -----
        S = complex(cast(0, 'like', T.Godard.metric));
        for pol = 1:NPol
            for kk = kLo:kHi
                a   = cast(R_corr_all(kk,         pol), 'like', T.Godard.metric);
                b   = cast(R_corr_all(kk + shift, pol), 'like', T.Godard.metric);
                S(:) = S + a * conj(b);
            end
        end
        e = cast(imag(S), 'like', T.Godard.ek);

        % --- PI update for the next block ------------------------------
        e_lf    = cast(e, 'like', T.Godard.lf);
        LF_I(:) = LF_I + ki_fi * e_lf;
        tauSamp(:) = kp_fi * e_lf + LF_I;
    end

    %% Remove cyclic extension
    DInit = 1 + (NExtra + NOverlap) / 2;
    DFin  = (NExtra - NOverlap) / 2;
    z     = OutPad(DInit : end - DFin, :);

    %% ================================================================
    %  4. Adaptive butterfly equaliser
    %  ================================================================
    if NPol == 1
        zEq = complex(zeros(size(z, 1), 2, 'like', T.AdaptEq.x));
        for i = 1:size(z, 1)
            zEq(i, 1) = cast(z(i, 1), 'like', T.AdaptEq.x);
            zEq(i, 2) = cast(z(i, 1), 'like', T.AdaptEq.x);
        end
        singlePol = true;
    else
        zEq = cast(z, 'like', T.AdaptEq.x);
        singlePol = false;
    end

    yEq = adaptive_eq.equalize_fxp(zEq, SpS, ...
        AdaptOpts.NTaps, AdaptOpts.Mu, AdaptOpts.SingleSpike, ...
        AdaptOpts.N1, AdaptOpts.NOut, AdaptOpts.SignOnly, ...
        AdaptOpts.UpdateStep, T.AdaptEq, AdaptOpts.PLanes, ...
        AdaptOpts.Mode, AdaptOpts.Pilots, AdaptOpts.BlockLen, ...
        AdaptOpts.SubframeBlocks);

    if singlePol
        y = yEq(:, 1);
    else
        y = yEq;
    end

    %% Trim to NSymb if requested
    if ~isempty(NSymb) && NSymb > 0
        if size(y, 1) > NSymb
            y = y(1:NSymb, :);
        end
    end
end
