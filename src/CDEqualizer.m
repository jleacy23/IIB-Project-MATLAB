classdef CDEqualizer
    properties
        D           % Dispersion [ps/(nm*km)]
        L           % Fiber length [km]
        CLambda     % Central wavelength [nm]
        Rs          % Symbol rate [GBd]
        NPol        % Number of polarizations
        SpSIn       % Samples per symbol
        NFFT        % FFT size
        NOverlap    % Overlap size (computed separately)
        WL          % fixed point word length
        FL          % fixed point fraction length
        c = 299792458; % speed of light
        fi_t        % fimath types
    end

    methods
        %% Constructor
        function obj = CDEqualizer(D,L,CLambda,Rs,NPol,SpSIn,NFFT, WL, FL)
            obj.D       = D * 1e-6;
            obj.L       = L * 1e3;
            obj.CLambda = CLambda * 1e-9;
            obj.Rs      = Rs * 1e9;
            obj.NPol    = NPol;
            obj.SpSIn   = SpSIn;
            obj.NFFT    = NFFT;
            obj.WL      = WL;
            obj.FL      = FL;
            obj.fi_t = fimath('RoundingMethod','Nearest', ...
                'OverflowAction','Saturate', ...
                'ProductMode','SpecifyPrecision', ...
                'ProductWordLength',WL, ...
                'ProductFractionLength',FL, ...
                'SumMode','SpecifyPrecision', ...
                'SumWordLength',WL, ...
                'SumFractionLength',FL);

            if nargin < 10 || isempty(NOverlap)
                obj.NOverlap = obj.computeOverlap();
            else
                obj.NOverlap = NOverlap + mod(NOverlap,2);
            end
        end

        %% Compute overlap (separate method)
        function NOverlap = computeOverlap(obj)
            beta2 = -obj.D * obj.CLambda^2 / (2*pi*obj.c);
            NOverlap = ceil(6.67 * abs(beta2) * (obj.Rs^2) * obj.L * obj.SpSIn);
            NOverlap = NOverlap + mod(NOverlap,2); % Ensure even
            %constrain to sensible limits
            NOverlap = min(NOverlap, floor(obj.NFFT/2));
        end

        %% Overlap-save CD compensation
        function Out = equalize(obj, In, UseFixedPoint)

            if isempty(obj.NOverlap)
                error('Overlap not computed. Call computeOverlap first.');
            end

            if UseFixedPoint
                fft_fi = dsp.FFT("OutputDataType", "Same as input");
                ifft_fi = dsp.IFFT("OutputDataType", "Same as input");
            end

            %% Frequency response
            n = (-obj.NFFT/2:obj.NFFT/2-1)';
            fN = obj.SpSIn * obj.Rs / 2;

            HCD = exp(-1i*pi*obj.CLambda^2*obj.D*obj.L/obj.c * ...
                     (n*2*fN/obj.NFFT).^2);

            if obj.NPol == 2
                HCD = cat(3,HCD,HCD);
            end

            if UseFixedPoint
                HCD = fi(HCD,1, obj.WL, obj.FL, 'fimath', obj.fi_t);
            end

            %% Input extension
            AuxLen = size(In,1)/(obj.NFFT-obj.NOverlap);

            if AuxLen ~= ceil(AuxLen)
                NExtra = ceil(AuxLen)*(obj.NFFT-obj.NOverlap) ...
                         - size(In,1);
                In = [In(end-NExtra/2+1:end,:);
                      In;
                      In(1:NExtra/2,:)];
            else
                NExtra = obj.NOverlap;
                In = [In(end-NExtra/2+1:end,:);
                      In;
                      In(1:NExtra/2,:)];
            end

            %% Block formation
            BlocksV = reshape(In(:,1), ...
                        obj.NFFT-obj.NOverlap, ...
                        size(In,1)/(obj.NFFT-obj.NOverlap));

            if obj.NPol == 2
                BlocksH = reshape(In(:,2), ...
                            obj.NFFT-obj.NOverlap, ...
                            size(In,1)/(obj.NFFT-obj.NOverlap));

                Blocks = cat(3,BlocksV,BlocksH);
            else
                Blocks = BlocksV;
            end

            %% Processing
            Out = zeros(size(Blocks));
            Overlap = zeros(obj.NOverlap,1,obj.NPol);

            for i = 1:size(Blocks,2)

                InB = [Overlap; Blocks(:,i,:)];

                if UseFixedPoint
                    InB = fi(InB,1, obj.WL, obj.FL, 'fimath', obj.fi_t);
                end
                
                if UseFixedPoint
                    InBFreq = fftshift(fft_fi(InB));
                else                    
                    InBFreq = fftshift(fft(InB));
                end
                OutFDEFreq = InBFreq .* HCD;

                if UseFixedPoint
                    OutFDE = ifft_fi(ifftshift(OutFDEFreq));
                else
                    OutFDE = ifft(ifftshift(OutFDEFreq));
                end

                


                Overlap = InB(end-obj.NOverlap+1:end,1,:);

                OutB = OutFDE( ...
                        obj.NOverlap/2+1:end-obj.NOverlap/2,1,:);

                Out(:,i,:) = OutB;
            end

            %% Reassemble output
            OutV = reshape(Out(:,:,1),[],1);

            if obj.NPol == 2
                OutH = reshape(Out(:,:,2),[],1);
                Out = [OutV OutH];
            else
                Out = OutV;
            end

            %% Remove extra samples
            DInit = 1 + (NExtra + obj.NOverlap)/2;
            DFin  = (NExtra - obj.NOverlap)/2;

            Out = Out(DInit:end-DFin,:);
        end
    end
end

