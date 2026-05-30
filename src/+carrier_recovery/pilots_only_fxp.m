function [v, ThetaPU] = pilots_only_fxp(x, NPol, BlockLen, Pilots, CordicIts, T) %#codegen
%PILOTS_ONLY_FXP  Fixed-point pilot-only carrier phase recovery.
%
%   [v, ThetaPU] = pilots_only_fxp(x, NPol, BlockLen, Pilots, CordicIts, T)
%
%   Fixed-point equivalent of carrier_recovery.pilots_only.
%   One phase estimate is obtained per block from pilot correlation at the
%   block start and then held constant over the full block.
%
%   Phase estimate per block, formed INDEPENDENTLY per polarisation:
%       theta_blk(pol) = angle(conj(Pilot(pol)) * x(blockStart, pol))
%
%   After a butterfly equaliser (CMA/RDE) the two polarisations carry
%   different, slowly-drifting carrier phases (the equaliser's per-pol phase
%   ambiguity adapts independently), so collapsing the pols into a single
%   shared phase (angle(sum_pol(...))) tracks neither and floors the BER.
%
%   Both the per-block angle estimate and the final phase correction use
%   CORDIC (cordic.vectoring / cordic.rotate); CordicIts sets the angular
%   resolution (~atan(2^-CordicIts)) and, in the bit-width sweep, equals the
%   swept fractional-bit precision.

    if nargin < 6 || isempty(T)
        T = carrier_recovery.fxp_types('fixed16');
    end

    Nsym    = size(x, 1);
    NBlocks = ceil(Nsym / BlockLen);

    x_fi      = cast(x, 'like', T.x);
    Pilots_fi = cast(Pilots, 'like', T.x);

    ThetaBlk = zeros(NBlocks, NPol, 'like', T.theta);

    % One pilot-based phase estimate per block, per polarisation.
    for b = 1:NBlocks
        blockStart = (b - 1) * BlockLen + 1;
        if blockStart <= Nsym
            for pol = 1:NPol
                rx = x_fi(blockStart, pol);

                pilot_re =  cast(real(Pilots_fi(b, pol)), 'like', T.acc);
                pilot_im = -cast(imag(Pilots_fi(b, pol)), 'like', T.acc);
                rx_re    =  cast(real(rx), 'like', T.acc);
                rx_im    =  cast(imag(rx), 'like', T.acc);

                corr_re = cast(pilot_re * rx_re - pilot_im * rx_im, 'like', T.acc);
                corr_im = cast(pilot_re * rx_im + pilot_im * rx_re, 'like', T.acc);

                ThetaBlk(b, pol) = cordic.vectoring(corr_re, corr_im, CordicIts, T);
            end
        end
    end

    % Hold each polarisation's phase estimate over its block.
    ThetaPU = zeros(Nsym, NPol, 'like', T.theta);
    for b = 1:NBlocks
        iStart = (b - 1) * BlockLen + 1;
        iEnd   = min(b * BlockLen, Nsym);
        for pol = 1:NPol
            ThetaPU(iStart:iEnd, pol) = ThetaBlk(b, pol);
        end
    end

    % Final phase correction via CORDIC rotation.
    v = complex(zeros(Nsym, NPol, 'like', T.x));
    for i = 1:Nsym
        for pol = 1:NPol
            [vr, vi] = cordic.rotate(real(x_fi(i, pol)), imag(x_fi(i, pol)), ...
                                     -ThetaPU(i, pol), CordicIts, T);
            v(i, pol) = complex(vr, vi);
        end
    end
end
