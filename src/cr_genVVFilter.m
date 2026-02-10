function w = cr_genVVFilter(Linewidth, Rs, SNR, SymbolEnergy, NPol, NTaps)
%CR_GENVVFILTER  Generate Wiener-optimised Viterbi-Viterbi filter.
%
%   w = cr_genVVFilter(Linewidth, Rs, SNR, SymbolEnergy, NPol, NTaps)
%
%   Inputs
%     Linewidth    - laser linewidth [Hz]
%     Rs           - symbol rate [GBd]
%     SNR          - signal-to-noise ratio [dB]
%     SymbolEnergy - average symbol energy
%     NPol         - number of polarizations
%     NTaps        - number of past/future symbols for phase estimation

    Rs_si = Rs * 1e9;  % GBd -> Bd

    L_filt = 2 * NTaps + 1;
    Ts     = 1 / Rs_si;
    VarDeltaPhi = 2 * pi * Linewidth * Ts;

    % Additive noise variance
    SNRLin = 10^(SNR/10) * 2 * 125e9 / (NPol * Rs_si);
    VarEta = SymbolEnergy / (2 * SNRLin);

    % K matrix
    KAux = zeros(NTaps);
    K    = zeros(L_filt);
    for i = 0:NTaps
        for j = 0:NTaps
            KAux(i+1, j+1) = min(i, j);
        end
    end
    K(1:NTaps+1, 1:NTaps+1) = rot90(KAux, 2);
    K(NTaps+1:L_filt, NTaps+1:L_filt) = KAux;

    I = eye(L_filt);
    C = SymbolEnergy^4 * 16 * VarDeltaPhi * K + ...
        SymbolEnergy^3 * 16 * VarEta * I;
    w = (ones(L_filt,1)' / C).';
    w = w / max(w);
end
