function [v, ThetaPU] = cr_viterbiViterbi(x, NPol, VVFilter, L, Pilots, PilotThreshold, UsePilots, BlockBased)
%CR_VITERBIVITERBI  Viterbi-Viterbi carrier phase estimation & correction.
%
%   v = cr_viterbiViterbi(x, NPol, NTaps, VVFilter)
%
%   Inputs
%     x             - input signal [samples x NPol]
%     NPol          - number of polarizations
%     VVFilter      - Viterbi-Viterbi filter coefficients
%     L             - block length
%     Pilots        - Pilot symbols at the start of every block [Pilot length x 1]
%     PilotThreshold - threshold to reverse a cycle slip
%     UsePilots      - whether to use pilots for phase estimation
%     BlockBased     - whether to apply phase correction on a block of symbols
    L_filt = length(VVFilter);
    ThetaML4 = zeros(size(x,1), NPol);

    % Reference phase from pilots
    NBlocks = ceil(size(x,1) / L);
    PhiRef = zeros(NBlocks, NPol);
    P = length(Pilots);

    % Correlate received symbols with pilots and find argument
    for pol = 1:NPol
        for b = 1:NBlocks
            blockStart = (b-1)*L + 1;
            blockEnd = (b-1)*L + P;
            block = x(blockStart:blockEnd, pol);
            corr = sum(conj(Pilots) .* block);
            PhiRef(b, pol) = angle(corr);
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

    if BlockBased
        ThetaML = repelem(ThetaML(1:L:end, :), L, 1);
        % crop to original length
        ThetaML = ThetaML(1:size(x,1), :);
    end

    % Phase unwrapping
    N = size(ThetaML, 1);
    ThetaPU = zeros(N, NPol);
    for pol = 1:NPol
        ThetaPrev = 0;
        for i = 1:N
            n = floor(0.5 + (ThetaPrev - ThetaML(i, pol)) / (pi/2));
            ThetaPU(i, pol) = ThetaML(i, pol) + n * (pi/2);
            if UsePilots
                BlockIdx = ceil(i / L);
                n = 0;
                PhaseDiff = ThetaPU(i, pol) - PhiRef(BlockIdx, pol);
                if PhaseDiff > pi /2
                    n = 1;
                elseif PhaseDiff < -pi/2
                    n = -1;
                end
                % ThetaPU(i, pol) = ThetaPU(i, pol) - pi/2 * round((ThetaPU(i, pol) - PhiRef(BlockIdx, pol)) / PilotThreshold);
                ThetaPU(i, pol) = ThetaPU(i, pol) - n * pi/2;
            end
            ThetaPrev = ThetaPU(i, pol);
        end
    end


    v = x .* exp(-1j*ThetaPU);
end
