function v = cr_viterbiViterbi(x, NPol, NTaps, VVFilter)
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
    ThetaML = ThetaML4 / 4 - pi/4; % Scale back down and rotate by pi/4 to align with constellation

    % Unwrap phase
    ThetaPU = zeros(size(ThetaML));
    ThetaPrev = zeros(1, NPol);

    for i = 1:size(ThetaML, 1)
        n = floor(1/2 + (ThetaPrev - ThetaML(i,:)) / (pi/2));
        ThetaPU(i,:) = ThetaML(i,:) + n * (pi/2);
        ThetaPrev = ThetaPU(i,:);
    end

    v = x .* exp(-1j*ThetaPU);
end
