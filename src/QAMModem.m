classdef QAMModem
    properties
        M
        N_pol
        bitsPerSymbol
    end

    methods
        %% Constructor
        function obj = QAMModem(M, N_pol)
            if mod(log2(M),1) ~= 0
                error('M must be power of 2.');
            end
            obj.M = M;
            obj.N_pol = N_pol;
            obj.bitsPerSymbol = log2(M);
        end

        %% Generate random bits
        function bits = randomBits(~, Nbits)
            bits = randi([0 1], Nbits, 1);
        end

        %% Bits -> QAM symbols
        function symbols = modulate(obj, bits)
            k = obj.bitsPerSymbol;
            bits = bits(:);

            Ns_total = floor(length(bits) / k);
            Ns_pol   = floor(Ns_total / obj.N_pol);

            if Ns_pol == 0
                error('Not enough bits.');
            end

            bits = bits(1:Ns_pol * obj.N_pol * k);

            bits = reshape(bits, k, []).';
            symIdx = bi2de(bits, 'left-msb');

            syms = qammod(symIdx, obj.M, ...
                'UnitAveragePower', true, ...
                'InputType', 'integer');

            symbols = reshape(syms, Ns_pol, obj.N_pol);
        end

        %% Hard symbol decision
        function decidedSymbols = decideSymbols(obj, rxSymbols)
            % rxSymbols: [Ns x N_pol]

            [Ns, Np] = size(rxSymbols);

            if Np ~= obj.N_pol
                error('Polarization count mismatch.');
            end

            % Demod to nearest constellation point
            idx = qamdemod(rxSymbols(:), obj.M, ...
                'UnitAveragePower', true, ...
                'OutputType', 'integer');

            decidedSymbols = qammod(idx, obj.M, ...
                'UnitAveragePower', true, ...
                'InputType', 'integer');

            decidedSymbols = reshape(decidedSymbols, Ns, Np);
        end

        %% Symbols -> bits
        function bits = symbolsToBits(obj, symbols)
            k = obj.bitsPerSymbol;

            idx = qamdemod(symbols(:), obj.M, ...
                'UnitAveragePower', true, ...
                'OutputType', 'integer');

            bitsMat = de2bi(idx, k, 'left-msb');

            bits = reshape(bitsMat.', [], 1);
        end

        %% Rectangular pulse shaping
        function tx = rectPulse(~, symbols, sps)
            [Ns, Np] = size(symbols);

            tx = zeros(Ns*sps, Np);

            for p = 1:Np
                up = upsample(symbols(:,p), sps);
                h = ones(sps,1);
                tx(:,p) = conv(up, h, 'same');
            end
        end

        %% Root Raised Cosine shaping
        function tx = rrcPulse(~, symbols, sps, rolloff, span)
            [Ns, Np] = size(symbols);

            h = rcosdesign(rolloff, span, sps, 'sqrt');

            tx = zeros(Ns*sps, Np);

            for p = 1:Np
                up = upsample(symbols(:,p), sps);
                tx(:,p) = conv(up, h, 'same');
            end
        end
    end
end
