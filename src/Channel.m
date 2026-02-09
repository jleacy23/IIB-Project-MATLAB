classdef Channel
    %CHANNEL Summary of this class goes here
    %   Detailed explanation goes here

    properties
        L       %Fiber length [km]
        SNR     % Signal-to-Noise ratio in dB
        SpS     % samples per symbol
        Rs      % Symbol rate [GBd]
        D       % Dispersion coefficient [ps/nm/km]
        CWavelength % Central wavelength [nm]
        c = 299792458 % speed of light [m/s]
        DGDSpec % PMD coefficient
        N_pmd   % Number of PMD stages
        Linewidth % Laser linewidth [Hz]
    end

    methods
        function obj = Channel(L, SNR, SpS, Rs, D, CWavelength, DGDSpec, N_pmd, Linewidth)
            obj.L = L * 1e3;
            obj.SNR = SNR;
            obj.SpS = SpS;
            obj.Rs = Rs * 1e9;
            obj.D = D * 1e-6;
            obj.CWavelength = CWavelength * 1e-9;
            obj.DGDSpec = DGDSpec;
            obj.N_pmd = N_pmd;
            obj.Linewidth = Linewidth;
        end

        function Y = add_awgn(obj, X)
            Y = awgn(X, obj.SNR, 'measured');
        end

        function Y = add_chromatic_dispersion(obj, X)
            w = 2*pi*(-1/2:1/size(X,1):1/2-1/size(X,1)).'*obj.SpS*obj.Rs;
            G = exp(1i*((obj.D*obj.CWavelength^2)/(4*pi*obj.c))*obj.L*w.^2);
            Y = ifft(ifftshift(G.*fftshift(fft(X))));
        end

        function Y = add_pmd(obj, X)
            SDTau = sqrt(3*pi/8) * obj.DGDSpec;
            Tau = (SDTau*sqrt(obj.L*1e-3)/sqrt(obj.N_pmd))*1e-12;

            w = 2*pi*fftshift(-1/2:1/size(X,1):1/2-1/size(X,1)).'*obj.SpS*obj.Rs;

            % Random unitary matrices V and U that describes mode coupling:
            for i =1:obj.N_pmd
                [V(:,:,i),~,U(:,:,i)] = svd(randn(2) + 1i*randn(2));
            end

            Freq_EV = fft(X(:,1));
            Freq_EH = fft(X(:,2));

            for i = 1:obj.N_pmd
                U_herm = U(:,:,i)';
                E_1 = U_herm(1,1)*Freq_EV + U_herm(1,2)*Freq_EH;
                E_2 = U_herm(2,1)*Freq_EV + U_herm(2,2)*Freq_EH;

                E_1 = E_1.*exp(1i*w*Tau/2);
                E_2 = E_2.*exp(-1i*w*Tau/2);

                Freq_EV = V(1,1,i)*E_1 + V(1,2,i)*E_2;
                Freq_EH = V(2,1,i)*E_1 + V(2,2,i)*E_2;
            end

            Y(:,1) = ifft(Freq_EV);
            Y(:,2) = ifft(Freq_EH);
        end

        function Y = add_phase_noise(obj, X)
            var_phi = 2*pi*obj.Linewidth/obj.Rs;
            delta_phi = sqrt(var_phi)*randn(size(X,1),1);
            phi = cumsum(delta_phi);
            Y = X.*exp(1i*phi);
        end
    end
end
