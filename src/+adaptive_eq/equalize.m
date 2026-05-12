function y = equalize(x, SpS, NTaps, Mu, SingleSpike, N1, NOut, SignOnly)
%equalize  Adaptive butterfly equalization (CMA).
%
%   y = equalize(x, SpS, NTaps, Mu, SingleSpike, N1, NOut, SignOnly)
%
%   Inputs
%     x             - input signal [samples x 2]
%     SpS           - samples per symbol
%     NTaps         - number of FIR taps
%     Mu            - step size
%     SingleSpike   - true/false for single-spike initialisation
%     N1            - iteration to reinitialise y-pol weights
%     NOut          - samples to discard after equalisation
%     SignOnly      - if true, use sign(error) and complex-sign of y
%                     (sign(real(y)) + j*sign(imag(y))) in the update
%                     instead of the full multiplications

    % CMA radius
    R_CMA = sqrt(2);

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

        % CMA error and conjugate-output factor (optionally sign-reduced)
        if SignOnly
            e1  = sign(R_CMA - abs(y1(i))^2);
            e2  = sign(R_CMA - abs(y2(i))^2);
            yc1 = sign(real(y1(i))) - 1j*sign(imag(y1(i)));
            yc2 = sign(real(y2(i))) - 1j*sign(imag(y2(i)));
        else
            e1  = R_CMA - abs(y1(i))^2;
            e2  = R_CMA - abs(y2(i))^2;
            yc1 = conj(y1(i));
            yc2 = conj(y2(i));
        end

        % CMA update
        w1V = w1V + Mu*xV(:,i)*e1*yc1;
        w1H = w1H + Mu*xH(:,i)*e1*yc1;
        w2V = w2V + Mu*xV(:,i)*e2*yc2;
        w2H = w2H + Mu*xH(:,i)*e2*yc2;

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
