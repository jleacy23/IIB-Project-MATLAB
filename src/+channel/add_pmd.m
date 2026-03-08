function Y = add_pmd(X, L, SpS, Rs, DGDSpec, N_pmd)
%ADD_PMD  Apply polarisation-mode dispersion.
%
%   Y = add_pmd(X, L, SpS, Rs, DGDSpec, N_pmd)
%
%   Inputs
%     X       - input signal [samples x 2]
%     L       - fibre length [km]
%     SpS     - samples per symbol
%     Rs      - symbol rate [GBd]
%     DGDSpec - PMD coefficient [ps/sqrt(km)]
%     N_pmd   - number of PMD stages

    L  = L * 1e3;   % km -> m
    Rs = Rs * 1e9;   % GBd -> Bd

    SDTau = sqrt(3*pi/8) * DGDSpec;
    Tau   = (SDTau * sqrt(L*1e-3) / sqrt(N_pmd)) * 1e-12;

    w = 2*pi*fftshift(-1/2:1/size(X,1):1/2-1/size(X,1)).' * SpS * Rs;

    % Random unitary matrices V and U for mode coupling
    V = zeros(2,2,N_pmd);
    U = zeros(2,2,N_pmd);
    for i = 1:N_pmd
        [V(:,:,i), ~, U(:,:,i)] = svd(randn(2) + 1i*randn(2));
    end

    Freq_EV = fft(X(:,1));
    Freq_EH = fft(X(:,2));

    for i = 1:N_pmd
        U_herm = U(:,:,i)';
        E_1 = U_herm(1,1)*Freq_EV + U_herm(1,2)*Freq_EH;
        E_2 = U_herm(2,1)*Freq_EV + U_herm(2,2)*Freq_EH;

        E_1 = E_1 .* exp( 1i*w*Tau/2);
        E_2 = E_2 .* exp(-1i*w*Tau/2);

        Freq_EV = V(1,1,i)*E_1 + V(1,2,i)*E_2;
        Freq_EH = V(2,1,i)*E_1 + V(2,2,i)*E_2;
    end

    Y(:,1) = ifft(Freq_EV);
    Y(:,2) = ifft(Freq_EH);
end
