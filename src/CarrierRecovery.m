classdef CarrierRecovery
    properties
        Linewidth  % Laser linewidth [Hz]
        Rs         % Symbol rate [GBd]
        SNR        % SNR in dB
        SymbolEnergy
        NPol       % Number of polarizations
        NTaps      % Number of past and future symbols used for phase estimation
        WL        % Fixed-point word length
        FL        % Fixed-point fraction length
        fi_t      % Fimath object for fixed-point arithmetic
        VVFilter  % Viterbi-Viterbi filter coefficients
        CordicIts  % Number of CORDIC iterations
        UseFixedPoint % Flag to indicate whether to use fixed-point arithmetic
    end

    methods
        function obj = CarrierRecovery(Linewidth, Rs, SNR, SymbolEnergy, NPol, NTaps, WL, FL, CordicIts, UseFixedPoint)
            obj.Linewidth = Linewidth;
            obj.Rs = Rs * 1e9;
            obj.SNR = SNR;
            obj.SymbolEnergy = SymbolEnergy;
            obj.NPol = NPol;
            obj.NTaps = NTaps;
            obj.WL = WL;
            obj.FL = FL;
            obj.CordicIts = CordicIts;
            obj.UseFixedPoint = UseFixedPoint;
            obj.fi_t = fimath('RoundingMethod','Nearest', ...
                'OverflowAction','Saturate', ...
                'ProductMode','SpecifyPrecision', ...
                'ProductWordLength',WL, ...
                'ProductFractionLength',FL, ...
                'SumMode','SpecifyPrecision', ...
                'SumWordLength',WL, ...
                'SumFractionLength',FL);
            obj.VVFilter = obj.genVVFilter(obj.UseFixedPoint);
        end

        function w = genVVFilter(obj, UseFixedPoint)
            L = 2 * obj.NTaps + 1;
            Ts = 1 / obj.Rs;
            VarDeltaPhi = 2 * pi * obj.Linewidth * Ts;

            % additive noise variance
            SNRLin = 10^(obj.SNR/10) * 2 * 125.e9 / (obj.NPol * obj.Rs); % linear SNR per polarization
            VarEta = obj.SymbolEnergy / (2 * SNRLin);

            % K Matrix
            KAux = zeros(obj.NTaps);
            K = zeros(L);
            for i = 0:obj.NTaps
                for j = 0:obj.NTaps
                    KAux(i+1,j+1) = min(i,j);
                end
            end
            K(1:obj.NTaps+1,1:obj.NTaps+1) = rot90(KAux,2);
            K(obj.NTaps+1:L,obj.NTaps+1:L) = KAux;

            I = eye(L);
            C = obj.SymbolEnergy^4 * 16 * VarDeltaPhi * K + obj.SymbolEnergy^3 * 16 * VarEta * I;
            w = (ones(L,1)'/(C)).';
            w = w / max(w); % normalize
            if UseFixedPoint
                w = fi(w, 1, obj.WL, obj.FL, 'fimath', obj.fi_t);
            end
        end

        function v = ViterbiViterbi(obj, x, UseFixedPoint)
            L = 2 * obj.NTaps + 1;
            ThetaML4 = zeros(size(x,1), obj.NPol);
            if UseFixedPoint
                ThetaML4 = fi(ThetaML4, 1, obj.WL, obj.FL, 'fimath', obj.fi_t);
            end

            for pol = 1:obj.NPol
                xBlocks = [zeros(floor(L/2),1); x(:,pol); zeros(floor(L/2),1)];
                xBlocks = convmtx(xBlocks.', L);
                xBlocks = flipud(xBlocks(:,L:end-L+1));

                if UseFixedPoint
                    xBlocks = fi(xBlocks, 1, obj.WL, obj.FL, 'fimath', obj.fi_t);
                end

                xBlocks4 = xBlocks.^4;
                if UseFixedPoint
                    filtered = obj.VVFilter.' * xBlocks4;
                    ThetaML4(:,pol) = cordicatan2(imag(filtered), real(filtered), obj.CordicIts);
                else
                    ThetaML4(:,pol) = angle(obj.VVFilter.' * xBlocks4);
                end
            end
            clearvars xBlocks xBlocks4;

            % TODO: add pilot correction, unwrap will be needed here
            ThetaML = ThetaML4 / 4 - pi/4;
            if UseFixedPoint
                ThetaML = fi(ThetaML, 1, obj.WL, obj.FL, 'fimath', obj.fi_t);
            end

            if UseFixedPoint
                x = fi(x, 1, obj.WL, obj.FL, 'fimath', obj.fi_t);
                v = cordicrotate(-ThetaML, x, obj.CordicIts);
            else
                v = x .* exp(-1j*ThetaML);
            end
        end
    end
end






