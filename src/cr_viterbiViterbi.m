function v = cr_viterbiViterbi(x, NPol, NTaps, VVFilter, Pilots, P, L, CSThreshold, UsePilots)
%CR_VITERBIVITERBI  Viterbi-Viterbi carrier phase estimation & correction.
%
%   v = cr_viterbiViterbi(x, NPol, NTaps, VVFilter)
%
%   Inputs
%     x             - input signal [samples x NPol]
%     NPol          - number of polarizations
%     NTaps         - number of past/future symbols for phase estimation
%     VVFilter      - Viterbi-Viterbi filter coefficients
%     Pilots        - pilot symbols [num pilots x NPol]
%     P             - number of pilots per block
%     L             - block length (symbols)
%     CSThreshold   - phase deviation to declare a cycle slip (radians)
%     UsePilots     - flag to enable pilot-based cycle slip correction

    L_filt = 2 * NTaps + 1;
    ThetaML4 = zeros(size(x,1), NPol);

    NBlocks = ceil(size(x,1) / L);
    PhiRef = zeros(NBlocks, NPol);

    % correlate pilots with received symbols to get reference phase (filters out additive noise)
    for pol = 1:NPol
        for b = 1:NBlocks
            PilotBlock = Pilots((b-1)*P+1 : min(b*P, size(Pilots,1)), pol);
            xBlock = x((b-1)*L+1 : min(((b-1)*L+P), size(x,1)), pol);
            PhiRef(b, pol) = angle(PilotBlock' * xBlock);
        end
    end 

    

    for pol = 1:NPol
        xBlocks = [zeros(floor(L_filt/2), 1); x(:,pol); zeros(floor(L_filt/2), 1)];
        xBlocks = convmtx(xBlocks.', L_filt);
        xBlocks = flipud(xBlocks(:, L_filt:end-L_filt+1));

        xBlocks4 = xBlocks.^4;
        ThetaML4(:,pol) = angle(VVFilter.' * xBlocks4);
    end

    % Phase correction
    ThetaML = ThetaML4 / 4 - pi/4;

    % Unwrap phase with pilot-aided CS correction
    ThetaBlocks = reshape(ThetaML, L, NBlocks, NPol);
    ThetaBlocksPU = zeros(size(ThetaBlocks));
    for pol = 1:NPol
        ThetaPrev = 0;
        for b = 1:NBlocks
            ThetaBlock = ThetaBlocks(:, b, pol);
            PhiRefBlock = PhiRef(b, pol);
            for i = 1:length(ThetaBlock)
                n = floor(1/2 + (ThetaPrev - ThetaBlock(i)) / (pi/2));
                ThetaBlocksPU (i, b, pol) = ThetaBlock(i) + n * (pi/2);
                if UsePilots
                    ThetaBlocksPU (i, b, pol) = ThetaBlocksPU(i, b, pol) - pi/2 * round((ThetaBlocksPU(i, b, pol) - PhiRefBlock) / CSThreshold);
                end
                ThetaPrev = ThetaBlocksPU(i, b, pol);
            end
        end
    end

    ThetaPU = reshape(ThetaBlocksPU, size(ThetaML));
    v = x .* exp(-1j*ThetaPU);
end
