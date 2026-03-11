function [v, ThetaPU] = bps(z, N, NPol, M, B, BlockLen, StepSize, ...
                                Pilots, PilotThreshold)
%bps  Blind Phase Search (BPS) carrier phase recovery.
%
%   [v, ThetaPU] = bps(z, N, NPol, M, B, BlockLen, StepSize,
%                           Pilots, PilotThreshold)
%
%   Inputs
%     z         - input signal [Nsym x NPol]
%     N         - one-sided BPS filter half-length; window = 2*N+1
%     NPol      - number of polarisations
%     M         - QAM order (4, 16, 64, ...)
%     B         - number of blind test phases (must be even)
%     BlockLen  - block length in symbols
%     StepSize  - phase update interval in symbols (1..BlockLen)
%                 Estimation and unwrapping execute every StepSize symbols,
%                 aligned to the start of each block.  Phase is held between
%                 updates.
%                   StepSize = 1        -> symbol-by-symbol (full bandwidth)
%                   StepSize = BlockLen -> one update per block
%     Pilots    - pilot symbols, one per block [NBlocks x NPol]
%                 The first symbol of block b is correlated against Pilots(b,:)
%                 for cycle-slip detection and correction.
%     PilotThreshold - threshold for pilot-based cycle-slip correction in radians
%
%   Outputs
%     v       - phase-corrected signal [Nsym x NPol]
%     ThetaPU - phase estimate [Nsym x NPol]
%               Step positions: estimated, unwrapped, pilot-corrected.
%               Between steps: held from the most recent step.

    %% ----------------------------------------------------------------
    %  Dimensions and test phase vector
    %% ----------------------------------------------------------------
    Nsym    = size(z, 1);
    L       = 2 * N + 1;
    NBlocks = ceil(Nsym / BlockLen);

    % Test phases uniformly covering (-pi/4, pi/4]
    b_vec     = (-B/2 : B/2-1);                % [-B/2 .. B/2-1]
    ThetaTest = (pi/2) * b_vec / B;            % [1 x B]

    %% ================================================================
    %  Step 1 – Pilot correlation  (vectorised per block)
    %
    %  Coherent combining across pols before taking angle(), matching the
    %  fixed-point version.  Each block's pilot slice is a matrix multiply.
    %% ================================================================
    PhiRef = zeros(NBlocks, NPol);

    for b = 1:min(NBlocks, size(Pilots, 1))
        blockStart = (b - 1) * BlockLen + 1;
        if blockStart <= Nsym
            PhiRef(b, :) = angle(conj(Pilots(b, :)) .* z(blockStart, :));
        end
    end

    %% ================================================================
    %  Step 2 – BPS phase estimation  (fully vectorised)
    %
    %  For each polarisation:
    %    (a) Build the sliding window matrix via convmtx: [L x Nsym]
    %    (b) Replicate across B test phases:             [L x B x Nsym]
    %    (c) Rotate each window sample by exp(-j*ThetaTest(b))
    %    (d) Apply QAM slicer element-wise
    %    (e) Sum squared errors over the window (dim 1): [1 x B x Nsym]
    %    (f) argmin over B -> Thetas(:, pol)
    %
    %  Thetas is computed for every symbol. The unwrapper below only
    %  uses step-position values; computing all keeps this fully vectorised.
    %
    %  Memory note: the [L x B x Nsym] arrays are the dominant cost.
    %  For large Nsym, B, L this may be significant — the sequential
    %  fixed-point implementation avoids this by computing one symbol at
    %  a time, at the cost of loop overhead.
    %% ================================================================
    Thetas = zeros(Nsym, NPol);

    % Rotation matrix: [L x B] (same for all symbols and pols)
    RotMat = repmat(exp(-1j * ThetaTest), L, 1);   % [L x B]

    for pol = 1:NPol
        % Build sliding window matrix [L x Nsym]
        zPad    = [zeros(floor(L/2), 1); z(:, pol); zeros(floor(L/2), 1)];
        zBlocks = convmtx(zPad.', L);
        zBlocks = flipud(zBlocks(:, L:end-L+1));     % [L x Nsym]

        % Rotate: [L x B x Nsym]  (broadcast over Nsym)
        % zRot(:, b, i) = zBlocks(:, i) * exp(-j*ThetaTest(b))
        zRot = bsxfun(@times, reshape(zBlocks, L, 1, Nsym), ...
                               reshape(RotMat,  L, B, 1));

        % QPSK decision and squared-error metric [L x B x Nsym]
        zDec = modem.slicer(zRot);
        m    = sum(abs(zRot - zDec).^2, 1);          % [1 x B x Nsym]

        % argmin over B -> [1 x 1 x Nsym]
        [~, im]       = min(m, [], 2);
        Thetas(:, pol) = ThetaTest(im(:)).';
    end

    %% ================================================================
    %  Step 3 – Step-based phase unwrapping + pilot correction
    %
    %  Same structure as carrier_recovery.viterbiViterbi: sequential loop for the
    %  stateful unwrapper, but at step positions only.
    %  Between steps, ThetaPrev is replicated and NOT updated, so the
    %  next step's unwrapper anchors to the last real estimate.
    %
    %  All-pol operations at each step are vectorised across pols.
    %% ================================================================
    ThetaPU   = zeros(Nsym, NPol);
    ThetaPrev = zeros(1, NPol);

    for i = 1:Nsym
        posInBlock = mod(i - 1, BlockLen);
        isStep     = (mod(posInBlock, StepSize) == 0);

        if isStep
            % Unwrap across all pols simultaneously
            n        = floor(0.5 - (Thetas(i, :) - ThetaPrev) / (pi/2));
            theta_uw = Thetas(i, :) + n * (pi/2);

            % Pilot-aided cycle-slip correction (per pol independently)
            BlockIdx  = ceil(i / BlockLen);
            % Wrap to (-pi, pi] for shortest-path error, threshold at ±pi/2
            PhaseDiff = mod(theta_uw - PhiRef(BlockIdx, :) + pi, 2*pi) - pi;
            n_slip    = round(PhaseDiff / PilotThreshold);
            theta_uw  = theta_uw - n_slip * (pi/2);

            ThetaPU(i, :) = theta_uw;
            ThetaPrev     = theta_uw;   % advance anchor to this step

        else
            % Hold: replicate last computed phase
            ThetaPU(i, :) = ThetaPrev;
        end
    end

    %% ================================================================
    %  Step 4 – Phase correction  (fully vectorised)
    %% ================================================================
    v = z .* exp(-1j * ThetaPU);
end