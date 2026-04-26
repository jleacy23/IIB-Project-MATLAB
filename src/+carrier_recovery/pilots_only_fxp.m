function [v, ThetaPU] = pilots_only_fxp(x, NPol, BlockLen, Pilots, ~, T) %#codegen
%PILOTS_ONLY_FXP  Fixed-point pilot-only carrier phase recovery.
%
%   [v, ThetaPU] = pilots_only_fxp(x, NPol, BlockLen, Pilots, CordicIts, T)
%
%   Fixed-point equivalent of carrier_recovery.pilots_only.
%   One phase estimate is obtained per block from pilot correlation at the
%   block start and then held constant over the full block.
%
%   Phase estimate per block (shared across polarisations):
%       theta_blk = angle(sum_pol(conj(Pilot) .* x(blockStart)))
%
%   The final phase correction is applied with complex multiplication.

    if nargin < 6 || isempty(T)
        T = carrier_recovery.fxp_types('fixed16');
    end

    ZERO_ACC = cast(0, 'like', T.acc);

    Nsym    = size(x, 1);
    NBlocks = ceil(Nsym / BlockLen);

    x_fi      = cast(x, 'like', T.x);
    Pilots_fi = cast(Pilots, 'like', T.x);

    ThetaBlk = zeros(NBlocks, 1, 'like', T.theta);

    % One pilot-based phase estimate per block.
    for b = 1:NBlocks
        blockStart = (b - 1) * BlockLen + 1;
        if blockStart <= Nsym
            corr_re = ZERO_ACC;
            corr_im = ZERO_ACC;
            for pol = 1:NPol
                rx = x_fi(blockStart, pol);

                pilot_re =  cast(real(Pilots_fi(b, pol)), 'like', T.acc);
                pilot_im = -cast(imag(Pilots_fi(b, pol)), 'like', T.acc);
                rx_re    =  cast(real(rx), 'like', T.acc);
                rx_im    =  cast(imag(rx), 'like', T.acc);

                corr_re = corr_re + (pilot_re * rx_re - pilot_im * rx_im);
                corr_im = corr_im + (pilot_re * rx_im + pilot_im * rx_re);
            end

            ThetaBlk(b) = cast(atan2(double(corr_im), double(corr_re)), 'like', T.theta);
        end
    end

    % Hold phase estimate over each block for all polarisations.
    ThetaPU = zeros(Nsym, NPol, 'like', T.theta);
    for b = 1:NBlocks
        iStart = (b - 1) * BlockLen + 1;
        iEnd   = min(b * BlockLen, Nsym);
        for pol = 1:NPol
            ThetaPU(iStart:iEnd, pol) = ThetaBlk(b);
        end
    end

    % Final phase correction.
    v = complex(zeros(Nsym, NPol, 'like', T.x));
    for i = 1:Nsym
        for pol = 1:NPol
            s_d = complex(double(real(x_fi(i, pol))), double(imag(x_fi(i, pol))));
            v(i, pol) = cast(s_d * exp(1j * double(-ThetaPU(i, pol))), 'like', T.x);
        end
    end
end
