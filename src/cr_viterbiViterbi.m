function [v, ThetaPU] = cr_viterbiViterbi(x, NPol, NTaps, VVFilter)
%CR_VITERBIVITERBI  Viterbi-Viterbi carrier phase estimation & correction.
%
%   v = cr_viterbiViterbi(x, NPol, NTaps, VVFilter)
%
%   Inputs
%     x             - input signal [samples x NPol]
%     NPol          - number of polarizations
%     NTaps         - number of past/future symbols for phase estimation
%     VVFilter      - Viterbi-Viterbi filter coefficients

    L_filt = 2 * NTaps + 1;
    ThetaML4 = zeros(size(x,1), NPol);

    for pol = 1:NPol
        xBlocks = [zeros(floor(L_filt/2), 1); x(:,pol); zeros(floor(L_filt/2), 1)];
        xBlocks = convmtx(xBlocks.', L_filt);
        xBlocks = flipud(xBlocks(:, L_filt:end-L_filt+1));

        xBlocks4 = xBlocks.^4;
        ThetaML4(:,pol) = angle(VVFilter.' * xBlocks4);
    end

    % Phase correction
    ThetaML = ThetaML4 / 4 - pi/4;

    % Phase unwrapping
    N = size(ThetaML, 1);
    ThetaPU = zeros(N, NPol);
    for pol = 1:NPol
        ThetaPrev = 0;
        for i = 1:N
            n = floor(0.5 + (ThetaPrev - ThetaML(i, pol)) / (pi/2));
            ThetaPU(i, pol) = ThetaML(i, pol) + n * (pi/2);
            ThetaPrev = ThetaPU(i, pol);
        end
    end

    v = x .* exp(-1j*ThetaPU);
end
