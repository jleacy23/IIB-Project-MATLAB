function y = adeq_equalize(x, SpS, Eq, NTaps, Mu, SingleSpike, N1, N2, NOut)
%ADEQ_EQUALIZE  Adaptive butterfly equalization (CMA / RDE / CMA+RDE).
%
%   y = adeq_equalize(x, SpS, Eq, NTaps, Mu, SingleSpike, N1, N2, NOut)
%
%   Inputs
%     x             - input signal [samples x 2]
%     SpS           - samples per symbol
%     Eq            - algorithm: 'CMA', 'RDE', or 'CMA+RDE'
%     NTaps         - number of FIR taps
%     Mu            - step size
%     SingleSpike   - true/false for single-spike initialisation
%     N1            - iteration to reinitialise y-pol weights
%     N2            - iteration to switch CMA->RDE ([] if unused)
%     NOut          - samples to discard after equalisation

    % Flags
    CMAFlag  = false;
    RDEFlag  = false;
    CMAtoRDE = false;

    if strcmp(Eq, 'CMA')
        CMAFlag = true;
    elseif strcmp(Eq, 'RDE')
        RDEFlag = true;
    elseif strcmp(Eq, 'CMA+RDE')
        CMAFlag  = true;
        CMAtoRDE = true;
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
    x = [x(end-floor(NTaps/2)+1:end,:); x; x(1:floor(NTaps/2),:)];

    xV = convmtx(x(:,1).', NTaps);
    xH = convmtx(x(:,2).', NTaps);

    xV = xV(:, NTaps:SpS:end-NTaps+1);
    xH = xH(:, NTaps:SpS:end-NTaps+1);

    OutLength = floor((size(x,1) - NTaps + 1) / 2);

    %% Initialise outputs
    y1 = zeros(OutLength, 1);
    y2 = zeros(OutLength, 1);

    %% Initial filter coefficients
    w1V = zeros(NTaps, 1);
    w1H = zeros(NTaps, 1);
    w2V = zeros(NTaps, 1);
    w2H = zeros(NTaps, 1);

    if SingleSpike
        w1V(floor(NTaps/2)+1) = 1;
    end

    %% Adaptive equalisation loop
    for i = 1:OutLength
        % Compute outputs
        y1(i) = w1V'*xV(:,i) + w1H'*xH(:,i);
        y2(i) = w2V'*xV(:,i) + w2H'*xH(:,i);

        % Update coefficients
        if CMAFlag
            [w1V,w1H,w2V,w2H] = adeq_cma_update( ...
                xV(:,i), xH(:,i), y1(i), y2(i), ...
                w1V, w1H, w2V, w2H, R_CMA, Mu);

            % Switch CMA -> RDE
            if CMAtoRDE && i == N2
                CMAFlag = false;
                RDEFlag = true;
            end

        elseif RDEFlag
            [w1V,w1H,w2V,w2H] = adeq_rde_update( ...
                xV(:,i), xH(:,i), y1(i), y2(i), ...
                w1V, w1H, w2V, w2H, R_RDE, Mu);
        end

        % Reinitialisation for SingleSpike
        if i == N1 && SingleSpike
            w2H = conj(w1V(end:-1:1, 1));
            w2V = -conj(w1H(end:-1:1, 1));
        end
    end

    % Collect output and remove extra samples
    y = [y1, y2];
    y = y(1+NOut:end, :);
end
