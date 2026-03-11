function [v, ThetaPU] = pilots_only(x, NPol, BlockLen, Pilots)
%PILOTS_ONLY  Pilot-only carrier phase recovery.
%
%   [v, ThetaPU] = pilots_only(x, NPol, BlockLen, Pilots)
%
%   Estimates the carrier phase from the single pilot symbol at the start
%   of each block and applies that estimate as a constant correction over
%   the entire block.  No blind estimation (BPS / VV) is performed.
%
%   The phase estimate for block b is:
%
%       ThetaEst(b, pol) = angle( conj(Pilots(b,pol)) * x(blockStart, pol) )
%
%   This is exact (up to noise) because the CPON PRBS pilots have known
%   amplitude, so conj(P)*r = |P|^2 * exp(j*phi) and angle() recovers
%   phi directly in (-pi, pi].  No pi/2 unwrapping is needed.
%
%   Inputs
%     x        - input signal [Nsym x NPol]
%     NPol     - number of polarisations
%     BlockLen - block length in symbols
%     Pilots   - pilot symbols, one per block [NBlocks x NPol]
%                Pilots(b, pol) is the known pilot at the start of block b.
%
%   Outputs
%     v       - phase-corrected signal [Nsym x NPol]
%     ThetaPU - phase estimate [Nsym x NPol]
%               Constant (held) over each block at the value derived from
%               that block's pilot.

    Nsym    = size(x, 1);
    NBlocks = ceil(Nsym / BlockLen);

    %% ================================================================
    %  Step 1 – Pilot phase estimation  (vectorised per block)
    %
    %  ThetaEst(b, pol) = angle( conj(Pilots(b,:)) .* x(blockStart,:) )
    %  where blockStart = (b-1)*BlockLen + 1.
    %  This is the direct phase-noise readout; no ambiguity needs resolving.
    %% ================================================================
    ThetaEst = zeros(NBlocks, NPol);

    for b = 1:NBlocks
        blockStart = (b - 1) * BlockLen + 1;
        if blockStart <= Nsym
            ThetaEst(b, :) = angle(conj(Pilots(b, :)) .* x(blockStart, :));
        end
    end

    %% ================================================================
    %  Step 2 – Hold phase estimate over block  (vectorised)
    %
    %  Build a symbol-indexed phase array by repeating each block's
    %  estimate across its BlockLen symbols.
    %% ================================================================
    ThetaPU = zeros(Nsym, NPol);

    for b = 1:NBlocks
        iStart = (b - 1) * BlockLen + 1;
        iEnd   = min(b * BlockLen, Nsym);
        ThetaPU(iStart:iEnd, :) = repmat(ThetaEst(b, :), iEnd - iStart + 1, 1);
    end

    %% ================================================================
    %  Step 3 – Phase correction  (fully vectorised)
    %% ================================================================
    v = x .* exp(-1j * ThetaPU);
end
