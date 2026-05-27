function [y, cfoBinsApplied] = combined_cd_fd_gardner_adaptive(In, SpS, ...
        NFFT, NOverlap, D, L, CLambda, Rs, Rolloff, ki, kp, NSymb, ...
        NLanes, AdaptOpts, cfoEnable, po2Twiddle)
%COMBINED_CD_FD_GARDNER_ADAPTIVE  Frequency-domain CD + RRC matched filter
%   (overlap-save) with optional one-shot coarse CFO correction, then
%   time-domain Gardner DPLL timing recovery, then butterfly CMA
%   equaliser.
%
%   y = combined_cd_fd_gardner_adaptive(In, SpS, NFFT, NOverlap, D, L, ...
%           CLambda, Rs, Rolloff, ki, kp, NSymb, NLanes, AdaptOpts)
%   y = combined_cd_fd_gardner_adaptive(..., AdaptOpts, cfoEnable)
%   y = combined_cd_fd_gardner_adaptive(..., cfoEnable, po2Twiddle)
%
%   Inputs
%     In          - input signal [samples x 2]
%     SpS         - samples per symbol (must be 2 for Gardner)
%     NFFT        - overlap-save FFT block size
%     NOverlap    - overlap-save overlap length (even, < NFFT)
%     D, L, CLambda, Rs - dispersion / signal parameters
%     Rolloff     - RRC matched-filter roll-off
%     ki, kp      - Gardner DPLL loop-filter gains
%     NSymb       - number of transmitted symbols
%     NLanes      - clk_recovery parallel lanes per block
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

    if SpS ~= 2
        error('combined_cd_fd_gardner_adaptive:SpS', ...
              'Gardner DPLL requires SpS = 2.');
    end
    if nargin < 15 || isempty(cfoEnable)
        cfoEnable = false;
    end
    if nargin < 16 || isempty(po2Twiddle)
        po2Twiddle = false;
    end

    NPol = size(In, 2);

    %% One-shot coarse CFO correction (time-domain phasor) ---------
    cfoBinsApplied = 0;
    if cfoEnable
        [In, cfoBinsApplied] = eq_clk.coarse_cfo_fd(In, NFFT);
    end

    %% Combined frequency-domain CD + matched-filter mask -----------
    HCDshift = eq_clk.cd_fd_response(D, L, CLambda, Rs, SpS, NFFT);
    HMFshift = eq_clk.rrc_fd_response(Rolloff, NFFT, SpS);
    Hstatic  = ifftshift(HCDshift .* HMFshift);    % natural FFT order

    %% Input cyclic extension to integer block count ----------------
    AuxLen = size(In,1) / (NFFT - NOverlap);
    if AuxLen ~= ceil(AuxLen)
        NExtra = ceil(AuxLen)*(NFFT - NOverlap) - size(In,1);
    else
        NExtra = NOverlap;
    end
    In = [In(end-NExtra/2+1:end,:); In; In(1:NExtra/2,:)];

    BlocksV = reshape(In(:,1), NFFT - NOverlap, ...
                      size(In,1)/(NFFT - NOverlap));
    if NPol == 2
        BlocksH = reshape(In(:,2), NFFT - NOverlap, ...
                          size(In,1)/(NFFT - NOverlap));
        Blocks = cat(3, BlocksV, BlocksH);
    else
        Blocks = BlocksV;
    end

    nBlk = size(Blocks, 2);

    %% Overlap-save loop --------------------------------------------
    Out     = zeros(size(Blocks));
    Overlap = zeros(NOverlap, 1, NPol);

    for i = 1:nBlk
        InB = [Overlap; Blocks(:,i,:)];

        % FFT and CD + matched filter
        R     = fft.fft_flp(InB, false, po2Twiddle);
        Rfilt = R .* Hstatic;

        % IFFT and overlap-save save
        OutFDE  = fft.fft_flp(Rfilt, true, po2Twiddle);
        Overlap = InB(end-NOverlap+1:end, 1, :);
        OutB    = OutFDE(NOverlap/2+1:end-NOverlap/2, 1, :);
        Out(:,i,:) = OutB;
    end

    %% Reassemble output --------------------------------------------
    OutV = reshape(Out(:,:,1), [], 1);
    if NPol == 2
        OutH = reshape(Out(:,:,2), [], 1);
        z = [OutV OutH];
    else
        z = OutV;
    end

    DInit = 1 + (NExtra + NOverlap)/2;
    DFin  = (NExtra - NOverlap)/2;
    z = z(DInit:end-DFin, :);

    %% Gardner timing recovery (per polarisation) -------------------
    zClk = run_gardner_per_pol(z, NSymb, ki, kp, NLanes);

    %% Adaptive butterfly CMA --------------------------------------
    if NPol == 1
        zClk = [zClk, zClk];
        singlePol = true;
    else
        singlePol = false;
    end
    y = eq_clk.apply_adaptive_eq(zClk, SpS, AdaptOpts);
    if singlePol
        y = y(:, 1);
    end
end


% =====================================================================
function z = run_gardner_per_pol(x, NSymb, ki, kp, NLanes)
    NPol = size(x, 2);
    outCols = cell(1, NPol);
    minLen = inf;
    for p = 1:NPol
        out = clk_recovery.recovery(x(:, p), 'NRZ', NSymb, ki, kp, NLanes);
        outCols{p} = out;
        minLen = min(minLen, length(out));
    end
    z = zeros(minLen, NPol);
    for p = 1:NPol
        z(:, p) = outCols{p}(1:minLen);
    end
end
