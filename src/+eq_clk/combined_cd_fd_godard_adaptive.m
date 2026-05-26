function y = combined_cd_fd_godard_adaptive(In, SpS, NFFT, NOverlap, ...
        D, L, CLambda, Rs, Rolloff, ki, kp, NSymb, AdaptOpts)
%COMBINED_CD_FD_GODARD_ADAPTIVE  Frequency-domain CD + RRC matched filter
%   (overlap-save) with Godard timing recovery sharing the same FFT and
%   overlap, followed by a butterfly CMA equaliser.
%
%   y = combined_cd_fd_godard_adaptive(In, SpS, NFFT, NOverlap, D, L, ...
%           CLambda, Rs, Rolloff, ki, kp, NSymb, AdaptOpts)
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
%
%   The Godard Modified-Godard timing metric is evaluated on the
%   already-corrected spectrum (CD + MF + current phase ramp).  Its imag
%   part drives a per-block PI loop filter whose accumulated tau is
%   applied as the next-block phase ramp.

    NPol = size(In, 2);

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

    %% Overlap-save loop with embedded Godard PI ---------------------
    Out     = zeros(size(Blocks));
    Overlap = zeros(NOverlap, 1, NPol);

    LF_I    = 0;
    tauSamp = 0;

    for i = 1:nBlk
        InB = [Overlap; Blocks(:,i,:)];

        % Natural-order FFT (matches recovery_godard convention)
        R = fft(InB);

        % Apply CD + matched filter
        Rfilt = R .* Hstatic;

        % Apply current cumulative timing phase ramp
        ramp   = exp(-1j * 2*pi * k_idx * tauSamp / NFFT);
        R_corr = Rfilt .* ramp;

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
        OutFDE  = ifft(R_corr);
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
