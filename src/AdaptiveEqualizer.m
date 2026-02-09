classdef AdaptiveEqualizer
    properties
        % Equalizer parameters
        Eq          % 'CMA', 'RDE', 'CMA+RDE'
        NTaps       % Number of taps in the butterfly filter
        Mu          % Step-size for coefficient calculation
        SingleSpike % true/false for single-spike initialization
        N1          % Pre-initialization iterations (for SingleSpike)
        N2          % Iterations before switching from CMA to RDE
        NOut        % Number of samples to discard after equalization
        WL          % fixed point word length
        FL          % fixed point fraction length
        fi_t        % fimath types
    end

    methods
        %% Constructor
        function obj = AdaptiveEqualizer(ParamDE, WL, FL)
            obj.Eq          = ParamDE.Eq;
            obj.NTaps       = ParamDE.NTaps;
            obj.Mu          = ParamDE.Mu;
            obj.SingleSpike = ParamDE.SingleSpike;
            obj.N1          = ParamDE.N1;
            obj.NOut        = ParamDE.NOut;
            obj.WL          = WL;
            obj.FL          = FL;
            obj.fi_t = fimath('RoundingMethod','Nearest', ...
                'OverflowAction','Saturate', ...
                'ProductMode','SpecifyPrecision', ...
                'ProductWordLength',WL, ...
                'ProductFractionLength',FL, ...
                'SumMode','SpecifyPrecision', ...
                'SumWordLength',WL, ...
                'SumFractionLength',FL);

            if isfield(ParamDE,'N2')
                obj.N2 = ParamDE.N2;
            else
                obj.N2 = [];
            end
        end

        %% Main adaptive equalization method
        function y = equalize(obj, x, SpS, UseFixedPoint)
            % Adaptive equalization for dual-polarization input
            %
            % x: input signal [samples x 2]
            % SpS: samples per symbol

            % Flags
            CMAFlag    = false;
            RDEFlag    = false;
            CMAtoRDE  = false;

            if strcmp(obj.Eq,'CMA')
                CMAFlag = true;
            elseif strcmp(obj.Eq,'RDE')
                RDEFlag = true;
            elseif strcmp(obj.Eq,'CMA+RDE')
                CMAFlag    = true;
                CMAtoRDE  = true;
            else
                error('Unsupported equalizer type');
            end

            % CMA radius
            if CMAFlag
                if ~CMAtoRDE
                    R_CMA = 1;
                else
                    R_CMA = 1.32;
                end
            end

            % RDE radii
            if CMAtoRDE || RDEFlag
                R_RDE = [1/sqrt(5), 1, 3/sqrt(5)];
            end

            %% Input blocks (padding for convolution)
            x = [x(end-floor(obj.NTaps/2)+1:end,:); x ; x(1:floor(obj.NTaps/2),:)];

            xV = convmtx(x(:,1).', obj.NTaps);
            xH = convmtx(x(:,2).', obj.NTaps);

            if UseFixedPoint
                xV = fi(xV, 1, obj.WL, obj.FL, 'fimath', obj.fi_t);
                xH = fi(xH, 1, obj.WL, obj.FL, 'fimath', obj.fi_t);
            end

            xV = xV(:, obj.NTaps:SpS:end-obj.NTaps+1);
            xH = xH(:, obj.NTaps:SpS:end-obj.NTaps+1);

            OutLength = floor((size(x,1)-obj.NTaps+1)/2);

            %% Initialize outputs
            y1 = zeros(OutLength,1);
            y2 = zeros(OutLength,1);

            if UseFixedPoint
                y1 = fi(y1, 1, obj.WL, obj.FL, 'fimath', obj.fi_t);
                y2 = fi(y2, 1, obj.WL, obj.FL, 'fimath', obj.fi_t);
            end

            %% Initial filter coefficients
            w1V = zeros(obj.NTaps,1);
            w1H = zeros(obj.NTaps,1);
            w2V = zeros(obj.NTaps,1);
            w2H = zeros(obj.NTaps,1);

            if UseFixedPoint
                w1V = fi(w1V, 1, obj.WL, obj.FL, 'fimath', obj.fi_t);
                w1H = fi(w1H, 1, obj.WL, obj.FL, 'fimath', obj.fi_t);
                w2V = fi(w2V, 1, obj.WL, obj.FL, 'fimath', obj.fi_t);
                w2H = fi(w2H, 1, obj.WL, obj.FL, 'fimath', obj.fi_t);
            end

            if obj.SingleSpike
                w1V(floor(obj.NTaps/2)+1) = 1;
            end

            %% Adaptive equalization loop
            for i = 1:OutLength
                % Compute outputs
                y1(i) = w1V'*xV(:,i) + w1H'*xH(:,i);
                y2(i) = w2V'*xV(:,i) + w2H'*xH(:,i);

                % Update coefficients
                if CMAFlag
                    [w1V,w1H,w2V,w2H] = obj.CMA(xV(:,i),xH(:,i),y1(i),y2(i),...
                                                 w1V,w1H,w2V,w2H,R_CMA,obj.Mu);

                    % Switch CMA->RDE
                    if CMAtoRDE && i == obj.N2
                        CMAFlag = false;
                        RDEFlag = true;
                    end

                elseif RDEFlag
                    [w1V,w1H,w2V,w2H] = obj.RDE(xV(:,i),xH(:,i),y1(i),y2(i),...
                                                 w1V,w1H,w2V,w2H,R_RDE,obj.Mu);
                end

                % Reinitialization for SingleSpike
                if i == obj.N1 && obj.SingleSpike
                    w2H = conj(w1V(end:-1:1,1));
                    w2V = -conj(w1H(end:-1:1,1));
                end
            end

            % Collect output and remove extra samples
            y = [y1, y2];
            y = y(1+obj.NOut:end,:);
        end

        function [w1V,w1H,w2V,w2H] = CMA(~, xV,xH,y1,y2,w1V,w1H,w2V,w2H,R,Mu)
            w1V = w1V + Mu*xV*(R-abs(y1).^2)*conj(y1);
            w1H = w1H + Mu*xH*(R-abs(y1).^2)*conj(y1);
            w2V = w2V + Mu*xV*(R-abs(y2).^2)*conj(y2);
            w2H = w2H + Mu*xH*(R-abs(y2).^2)*conj(y2);
        end

        function [w1V, w1H, w2V, w2H] = RDE(~, xV,xH,y1,y2,w1V,w1H,w2V,w2H,R,Mu)
            [~,r1]= min(abs(R-abs(y1))) ; [~,r2] = min(abs(R-abs(y2)));
            w1V = w1V + Mu*xV*(R(r1)^2-abs(y1).^2)*conj(y1);
            w1H = w1H + Mu*xH*(R(r1)^2-abs(y1).^2)*conj(y1);
            w2V = w2V + Mu*xV*(R(r2)^2-abs(y2).^2)*conj(y2);
            w2H = w2H + Mu*xH*(R(r2)^2-abs(y2).^2)*conj(y2);
        end
    end
end