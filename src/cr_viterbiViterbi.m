function [v, ThetaPU] = cr_viterbiViterbi(x, NPol, VVFilter, BlockLen, ...
                                            StepSize, Pilots, UsePilots, PilotThreshold)
%CR_VITERBIVITERBI  Viterbi-Viterbi carrier phase estimation & correction.
%
%   [v, ThetaPU] = cr_viterbiViterbi(x, NPol, VVFilter, BlockLen,
%                                     StepSize, Pilots, UsePilots)
%
%   Inputs
%     x         - input signal [N x NPol]
%     NPol      - number of polarisations
%     VVFilter  - VV filter coefficients [(2*NTaps+1) x 1]
%     BlockLen  - block length in symbols
%                 Pilots are taken from the first P symbols of each block.
%     StepSize  - phase update interval in symbols (1..BlockLen)
%                 Estimation and unwrapping execute every StepSize symbols,
%                 aligned to the start of each block.  Phase is held between
%                 updates.
%                   StepSize = 1        -> symbol-by-symbol (full bandwidth)
%                   StepSize = BlockLen -> one update per block
%     Pilots    - pilot symbols at block start [P x 1]
%     UsePilots - logical: enable pilot-aided cycle-slip correction
%     PilotThreshold - threshold for pilot-based cycle-slip correction in radians
%
%   Outputs
%     v       - phase-corrected signal [N x NPol]
%     ThetaPU - phase estimate [N x NPol]
%               Step positions: estimated, unwrapped, pilot-corrected.
%               Between steps: held from the most recent step.

    %% ----------------------------------------------------------------
    %  Dimensions
    %% ----------------------------------------------------------------
    N      = size(x, 1);
    L_filt = length(VVFilter);
    P      = length(Pilots);
    NBlocks = ceil(N / BlockLen);

    %% ================================================================
    %  Step 1 – Pilot correlation  (vectorised per block)
    %
    %  PhiRef(b, pol) = angle( sum_p conj(Pilots(p)) * x(blockStart+p-1, pol) )
    %
    %  Each block's pilot window is extracted as a matrix slice so the
    %  inner product is a single matrix multiply rather than a loop.
    %% ================================================================
    PhiRef = zeros(NBlocks, NPol);

    if UsePilots
        for b = 1:NBlocks
            blockStart = (b - 1) * BlockLen + 1;
            idxEnd     = min(blockStart + P - 1, N);
            block      = x(blockStart:idxEnd, :);           % [P x NPol]
            % Coherent combining across pols then take angle — matches fxp
            corr = sum(conj(Pilots(1:size(block,1))) .* block, 1); % [1 x NPol]
            PhiRef(b, :) = angle(corr);
        end
    end

    %% ================================================================
    %  Step 2 – VV phase estimation  (fully vectorised)
    %
    %  For each polarisation:
    %    (a) Zero-pad and build the sliding-window matrix via convmtx.
    %    (b) Raise each column to the 4th power.
    %    (c) Filter: ThetaML4 = angle( VVFilter.' * xBlocks4 )
    %    (d) ThetaML = ThetaML4 / 4 - pi/4
    %
    %  ThetaML(:, pol) is computed for every symbol regardless of StepSize;
    %  the unwrapper below only uses values at step positions.  Computing all
    %  values keeps this section fully vectorised — the saving from skipping
    %  non-step symbols is smaller than the cost of the sequential loop that
    %  would be required to select them here.
    %% ================================================================
    ThetaML = zeros(N, NPol);

    for pol = 1:NPol
        xPad    = [zeros(floor(L_filt/2), 1); x(:, pol); zeros(floor(L_filt/2), 1)];
        xBlocks = convmtx(xPad.', L_filt);
        xBlocks = flipud(xBlocks(:, L_filt:end-L_filt+1));   % [L_filt x N]

        xBlocks4          = xBlocks .^ 4;
        ThetaML4_pol      = angle(VVFilter.' * xBlocks4);    % [1 x N]
        ThetaML(:, pol)   = (ThetaML4_pol / 4 - pi/4).';
    end

    %% ================================================================
    %  Step 3 – Step-based phase unwrapping + pilot correction
    %
    %  The unwrapper is inherently sequential (each output depends on the
    %  previous) so a loop is unavoidable.  However:
    %   - At step positions: unwrap, pilot-correct, update ThetaPrev.
    %   - Between steps: replicate ThetaPrev (no estimation, no unwrap).
    %  ThetaPrev is updated ONLY at step positions so that the next step's
    %  unwrapper anchors to the last real estimate, not to a held value.
    %
    %  Step positions within each block use the same 0-indexed boundary
    %  condition as the fixed-point version:
    %    posInBlock = mod(i-1, BlockLen)
    %    isStep     = mod(posInBlock, StepSize) == 0
    %  The first symbol of every block is always a step.
    %% ================================================================
    ThetaPU   = zeros(N, NPol);
    ThetaPrev = zeros(1, NPol);

    for i = 1:N
        posInBlock = mod(i - 1, BlockLen);
        isStep     = (mod(posInBlock, StepSize) == 0);

        if isStep
            % Unwrap all pols simultaneously (vectorised across pols)
            n        = floor(0.5 + (ThetaPrev - ThetaML(i, :)) / (pi/2));
            theta_uw = ThetaML(i, :) + n * (pi/2);

            % Pilot-aided cycle-slip correction (per pol independently)
            if UsePilots
                BlockIdx    = ceil(i / BlockLen);
                % Wrap to (-pi, pi] for shortest-path error, threshold at ±pi/2
                PhaseDiff   = mod(theta_uw - PhiRef(BlockIdx, :) + pi, 2*pi) - pi;
                n_slip      = round(PhaseDiff / PilotThreshold);
                theta_uw    = theta_uw - n_slip * (pi/2);
            end

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
    v = x .* exp(-1j * ThetaPU);
end