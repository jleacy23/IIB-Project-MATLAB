function [v, ThetaPU] = pilots_only_fxp(x, NPol, BlockLen, Pilots, CordicIts, T) %#codegen
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
%   The final phase correction is applied with cordicrotate and the input
%   angle is reduced to the principal range for robust CORDIC behaviour.

    if nargin < 6 || isempty(T)
        T = carrier_recovery.fxp_types('fixed16');
    end
    if nargin < 5 || isempty(CordicIts)
        CordicIts = 16;
    end

    CORDIC_ITS = coder.const(CordicIts);
    ZERO_ACC   = cast(0, 'like', T.acc);
    PI_VAL     = cast(pi, 'like', T.theta);
    PI_OVER2   = cast(pi/2, 'like', T.theta);

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

            corr_cplx = complex(cast(corr_re, 'like', T.theta), ...
                                cast(corr_im, 'like', T.theta));
            ThetaBlk(b) = cast(cordicangle(corr_cplx, CORDIC_ITS), 'like', T.theta);
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

    % Final phase correction via CORDIC.
    v = complex(zeros(Nsym, NPol, 'like', T.x));
    for i = 1:Nsym
        for pol = 1:NPol
            theta_d = mod(double(-ThetaPU(i, pol)) + pi, 2*pi) - pi;
            s_in    = x_fi(i, pol);

            if theta_d > pi/2
                theta_d = theta_d - pi;
                s_in    = -s_in;
            elseif theta_d < -pi/2
                theta_d = theta_d + pi;
                s_in    = -s_in;
            end

            theta_safe = cast(theta_d, 'like', T.theta);
            if theta_safe > PI_OVER2
                theta_safe = theta_safe - PI_VAL;
                s_in       = -s_in;
            elseif theta_safe < -PI_OVER2
                theta_safe = theta_safe + PI_VAL;
                s_in       = -s_in;
            end

            v(i, pol) = cast(cordicrotate(theta_safe, s_in, CORDIC_ITS), 'like', T.x);
        end
    end
end
