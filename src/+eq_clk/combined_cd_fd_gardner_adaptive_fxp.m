function [y, cfoBinsApplied] = combined_cd_fd_gardner_adaptive_fxp(...
        In, SpS, NFFT, NOverlap, D, L, CLambda, Rs, Rolloff, ki, kp, ...
        NSymb, NLanes, AdaptOpts, cfoEnable, po2Twiddle, T) %#codegen
%COMBINED_CD_FD_GARDNER_ADAPTIVE_FXP  Fixed-point CD-FD + Gardner DPLL +
%   butterfly CMA, with optional coarse CFO correction.
%
%   [y, cfoBinsApplied] = combined_cd_fd_gardner_adaptive_fxp(In, SpS, ...
%       NFFT, NOverlap, D, L, CLambda, Rs, Rolloff, ki, kp, NSymb, ...
%       NLanes, AdaptOpts, cfoEnable, po2Twiddle, T)
%
%   Drop-in fixed-point version of eq_clk.combined_cd_fd_gardner_adaptive.
%   Each pipeline section uses its own bit widths via the composite types
%   table T (eq_clk.combined_cd_fd_gardner_adaptive_fxp_types):
%
%     T.Static  - frequency-domain CD + matched filter (overlap-save)
%     T.Clk     - Gardner DPLL clock recovery
%     T.AdaptEq - adaptive butterfly equaliser
%
%   The optional coarse CFO correction is applied in floating point
%   (cast fi -> double, eq_clk.coarse_cfo_fd in float, cast back to
%   T.Static.x) — no dedicated fxp helper is used for it.
%
%   AdaptOpts is a struct of explicit fields (NTaps, Mu, SingleSpike, N1,
%   NOut, SignOnly, UpdateStep, PLanes, Mode, Pilots, BlockLen,
%   SubframeBlocks).  For codegen, the build script supplies a typed
%   prototype.  Pilots must be a fi matrix at T.AdaptEq.y precision
%   (use [] for pure CMA — pass a 0x2 fi matrix in codegen).
%
%   Required: SpS = 2 (Gardner DPLL).

    if SpS ~= 2
        error('combined_cd_fd_gardner_adaptive_fxp:SpS', ...
              'Gardner DPLL requires SpS = 2.');
    end
    if nargin < 15 || isempty(cfoEnable)
        cfoEnable = false;
    end
    if nargin < 16 || isempty(po2Twiddle)
        po2Twiddle = false;
    end
    if nargin < 17 || isempty(T)
        T = eq_clk.combined_cd_fd_gardner_adaptive_fxp_types('fixed32');
    end

    NPol = size(In, 2);

    %% ================================================================
    %  1. Coarse CFO correction (one-shot, time-domain phasor)
    %  ================================================================
    %  Runs in floating point.  Cast fi input to double, call the
    %  floating-point coarse_cfo_fd, then cast the corrected signal back
    %  to T.Static.x for the downstream fxp data path.
    cfoBinsApplied = 0;
    if cfoEnable
        InDouble = double(In);
        [InCFO_d, cfoBinsApplied] = eq_clk.coarse_cfo_fd(InDouble, NFFT);
        InS = cast(InCFO_d, 'like', T.Static.x);
    else
        InS = cast(In, 'like', T.Static.x);
    end

    %% ================================================================
    %  2. Static CD + matched filter (overlap-save, frequency domain)
    %  ================================================================
    HCDshift = eq_clk.cd_fd_response(D, L, CLambda, Rs, SpS, NFFT);
    HMFshift = eq_clk.rrc_fd_response(Rolloff, NFFT, SpS);
    Hstatic  = ifftshift(HCDshift .* HMFshift);     % natural FFT order
    Hstatic_fi = cast(Hstatic, 'like', T.Static.hcd);

    Tfft.x   = T.Static.x;
    Tfft.tw  = T.Static.tw;
    Tfft.acc = T.Static.acc;   % single FFT precision (uniform static config)

    %% Input cyclic extension to integer block count
    NIn    = size(InS, 1);
    stepLen = NFFT - NOverlap;
    AuxLen = NIn / stepLen;
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

    %% Prepend NOverlap zeros so a sliding window captures overlap-save
    %  without an explicit overlap buffer (matches cd_eq.equalize_fxp).
    InPad2 = [complex(zeros(NOverlap, NPol, 'like', T.Static.x)); InPad];
    nBlocks = NPadded / stepLen;

    OutPad = complex(zeros(NPadded, NPol, 'like', T.Static.acc));

    %% Per-polarisation overlap-save processing (avoid 3-D fi arrays)
    for pol = 1:NPol
        for i = 1:nBlocks
            wStart = (i - 1) * stepLen + 1;
            InB = InPad2(wStart : wStart + NFFT - 1, pol);

            X = fft.fft_fxp(InB, false, po2Twiddle, Tfft);

            Y = complex(zeros(NFFT, 1, 'like', T.Static.acc));
            for k = 1:NFFT
                Y(k) = X(k) * Hstatic_fi(k);
            end

            outFDE = fft.fft_fxp(Y, true, po2Twiddle, Tfft);

            oStart = (i - 1) * stepLen + 1;
            for j = 1:stepLen
                OutPad(oStart + j - 1, pol) = outFDE(halfOv + j);
            end
        end
    end

    %% Remove cyclic extension
    DInit = 1 + (NExtra + NOverlap) / 2;
    DFin  = (NExtra - NOverlap) / 2;
    z     = OutPad(DInit : end - DFin, :);

    %% ================================================================
    %  3. Gardner timing recovery (per polarisation)
    %  ================================================================
    NPolZ = size(z, 2);
    % Per-pol recovery, then truncate all pols to the shortest output so
    % the clock-recovery output is a tidy matrix for the adaptive EQ.
    zV_in = cast(z(:, 1), 'like', T.Clk.x);
    zV    = clk_recovery.recovery_fxp(zV_in, NSymb, ki, kp, NLanes, T.Clk);
    minLen = size(zV, 1);

    if NPolZ == 2
        zH_in = cast(z(:, 2), 'like', T.Clk.x);
        zH    = clk_recovery.recovery_fxp(zH_in, NSymb, ki, kp, NLanes, T.Clk);
        if size(zH, 1) < minLen
            minLen = size(zH, 1);
        end
        zClk = complex(zeros(minLen, 2, 'like', T.Clk.x));
        for i = 1:minLen
            zClk(i, 1) = zV(i);
            zClk(i, 2) = zH(i);
        end
    else
        zClk = complex(zeros(minLen, 1, 'like', T.Clk.x));
        for i = 1:minLen
            zClk(i, 1) = zV(i);
        end
    end

    %% ================================================================
    %  4. Adaptive butterfly equaliser
    %  ================================================================
    %  The adaptive_eq.equalize_fxp signature expects a 2-pol input.  For
    %  single-pol input, duplicate the column (matches the float wrapper).
    if NPol == 1
        zEq = complex(zeros(size(zClk, 1), 2, 'like', T.AdaptEq.x));
        for i = 1:size(zClk, 1)
            zEq(i, 1) = cast(zClk(i, 1), 'like', T.AdaptEq.x);
            zEq(i, 2) = cast(zClk(i, 1), 'like', T.AdaptEq.x);
        end
        singlePol = true;
    else
        zEq = cast(zClk, 'like', T.AdaptEq.x);
        singlePol = false;
    end

    % CMA radius scale (post/pre input-normalisation energy ratio).  Optional
    % field: defaults to 1 (unscaled R) when the caller does not supply it, so
    % AdaptOpts structs built before RScale existed still compile.
    if isfield(AdaptOpts, 'RScale')
        rScaleEq = double(AdaptOpts.RScale);
    else
        rScaleEq = 1;
    end

    yEq = adaptive_eq.equalize_fxp(zEq, SpS, ...
        AdaptOpts.NTaps, AdaptOpts.Mu, AdaptOpts.SingleSpike, ...
        AdaptOpts.N1, AdaptOpts.NOut, AdaptOpts.SignOnly, ...
        AdaptOpts.UpdateStep, T.AdaptEq, AdaptOpts.PLanes, ...
        AdaptOpts.Mode, AdaptOpts.Pilots, AdaptOpts.BlockLen, ...
        AdaptOpts.SubframeBlocks, rScaleEq);

    if singlePol
        y = yEq(:, 1);
    else
        y = yEq;
    end
end
