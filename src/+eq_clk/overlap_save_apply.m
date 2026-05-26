function Out = overlap_save_apply(In, H, NFFT, NOverlap)
%OVERLAP_SAVE_APPLY  Apply a static freq-domain mask via overlap-save.
%
%   Out = overlap_save_apply(In, H, NFFT, NOverlap)
%
%   In        - input [samples x NPol]
%   H         - frequency response (NFFT x 1) on the fftshifted grid
%   NFFT      - block size
%   NOverlap  - overlap length (even, < NFFT)
%
%   Mirrors the cyclic-extension overlap-save pattern of cd_eq.equalize so
%   the output length matches the input length.

    NPol = size(In, 2);
    if size(H, 3) == 1 && NPol == 2
        H = cat(3, H, H);
    end

    %% Input extension (cyclic) ---------------------------------------
    AuxLen = size(In,1) / (NFFT - NOverlap);
    if AuxLen ~= ceil(AuxLen)
        NExtra = ceil(AuxLen)*(NFFT - NOverlap) - size(In,1);
        In = [In(end-NExtra/2+1:end,:); In; In(1:NExtra/2,:)];
    else
        NExtra = NOverlap;
        In = [In(end-NExtra/2+1:end,:); In; In(1:NExtra/2,:)];
    end

    %% Block formation -----------------------------------------------
    BlocksV = reshape(In(:,1), NFFT - NOverlap, ...
                      size(In,1)/(NFFT - NOverlap));
    if NPol == 2
        BlocksH = reshape(In(:,2), NFFT - NOverlap, ...
                          size(In,1)/(NFFT - NOverlap));
        Blocks = cat(3, BlocksV, BlocksH);
    else
        Blocks = BlocksV;
    end

    %% Processing ----------------------------------------------------
    Out     = zeros(size(Blocks));
    Overlap = zeros(NOverlap, 1, NPol);

    for i = 1:size(Blocks,2)
        InB     = [Overlap; Blocks(:,i,:)];
        InBFreq = fftshift(fft(InB));
        OutFreq = InBFreq .* H;
        OutFDE  = ifft(ifftshift(OutFreq));
        Overlap = InB(end-NOverlap+1:end, 1, :);
        OutB    = OutFDE(NOverlap/2+1:end-NOverlap/2, 1, :);
        Out(:,i,:) = OutB;
    end

    %% Reassemble ----------------------------------------------------
    OutV = reshape(Out(:,:,1), [], 1);
    if NPol == 2
        OutH = reshape(Out(:,:,2), [], 1);
        Out  = [OutV OutH];
    else
        Out = OutV;
    end

    DInit = 1 + (NExtra + NOverlap)/2;
    DFin  = (NExtra - NOverlap)/2;
    Out   = Out(DInit:end-DFin, :);
end
