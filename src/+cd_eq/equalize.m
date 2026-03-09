function Out = equalize(In, D, L, CLambda, Rs, NPol, SpSIn, NFFT)
%EQUALIZE  Overlap-save frequency-domain CD compensation.
%
%   Out = equalize(In, D, L, CLambda, Rs, NPol, SpSIn, NFFT)
%
%   Inputs (user-friendly units)
%     In            - input signal [samples x NPol]
%     D             - dispersion coefficient [ps/(nm*km)]
%     L             - fibre length [km]
%     CLambda       - central wavelength [nm]
%     Rs            - symbol rate [GBd]
%     NPol          - number of polarizations
%     SpSIn         - samples per symbol
%     NFFT          - FFT block size

    c = 299792458;

    % Unit conversion
    D_si       = D * 1e-6;
    L_si       = L * 1e3;
    CLambda_si = CLambda * 1e-9;
    Rs_si      = Rs * 1e9;

    % Compute overlap
    NOverlap = cd_eq.computeOverlap(D, L, CLambda, Rs, SpSIn, NFFT);

    %% Frequency response
    n  = (-NFFT/2:NFFT/2-1)';
    fN = SpSIn * Rs_si / 2;

    HCD = exp(-1i*pi*CLambda_si^2*D_si*L_si/c * ...
              (n*2*fN/NFFT).^2);

    if NPol == 2
        HCD = cat(3, HCD, HCD);
    end

    %% Input extension
    AuxLen = size(In,1) / (NFFT - NOverlap);

    if AuxLen ~= ceil(AuxLen)
        NExtra = ceil(AuxLen)*(NFFT - NOverlap) - size(In,1);
        In = [In(end-NExtra/2+1:end,:); In; In(1:NExtra/2,:)];
    else
        NExtra = NOverlap;
        In = [In(end-NExtra/2+1:end,:); In; In(1:NExtra/2,:)];
    end

    %% Block formation
    BlocksV = reshape(In(:,1), NFFT - NOverlap, ...
                      size(In,1)/(NFFT - NOverlap));

    if NPol == 2
        BlocksH = reshape(In(:,2), NFFT - NOverlap, ...
                          size(In,1)/(NFFT - NOverlap));
        Blocks = cat(3, BlocksV, BlocksH);
    else
        Blocks = BlocksV;
    end

    %% Processing
    Out     = zeros(size(Blocks));
    Overlap = zeros(NOverlap, 1, NPol);

    for i = 1:size(Blocks,2)
        InB = [Overlap; Blocks(:,i,:)];

        InBFreq = fftshift(fft(InB));

        OutFDEFreq = InBFreq .* HCD;

        OutFDE = ifft(ifftshift(OutFDEFreq));

        Overlap = InB(end-NOverlap+1:end, 1, :);
        OutB    = OutFDE(NOverlap/2+1:end-NOverlap/2, 1, :);
        Out(:,i,:) = OutB;
    end

    %% Reassemble output
    OutV = reshape(Out(:,:,1), [], 1);

    if NPol == 2
        OutH = reshape(Out(:,:,2), [], 1);
        Out  = [OutV OutH];
    else
        Out = OutV;
    end

    %% Remove extra samples
    DInit = 1 + (NExtra + NOverlap)/2;
    DFin  = (NExtra - NOverlap)/2;
    Out   = Out(DInit:end-DFin, :);
end
