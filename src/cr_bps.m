function [v, ThetaPU] = cr_bps(z, N, NPol, M, B, BlockLen, Pilots, UsePilots, BlockBased)
    L=2*N+1;
    p = pi/2;

    % Reference phase from pilots
    NBlocks = ceil(size(z,1) / BlockLen);
    PhiRef = zeros(NBlocks, NPol);
    P = length(Pilots);

    for b = 1:NBlocks
        blockStart = (b-1)*BlockLen + 1;
        blockEnd = (b-1)*BlockLen + P;
        block = z(blockStart:blockEnd, :);
        corr = 0;
        for pol = 1:NPol
            corr = corr + sum(conj(Pilots) .* block(:, pol));
        end
        PhiRef(b, :) = angle(corr);
    end
    % Creating a vector of test carrier phase angles:
    b =-B/2:1:B/2-1;
    ThetaTest= p * (b/B);
    % Creating a matrix of test carrier phase angles:
    ThetaTestMatrix = repmat(exp(-1j*ThetaTest),L,1);
    if NPol== 2
        ThetaTestMatrix=cat(3,ThetaTestMatrix,ThetaTestMatrix);
    end
    % Inputblocks:
    % V-pol. orientation:
    zB_V =[zeros(floor(L/2),1);
    z(:,1);
    zeros(floor(L/2),1)];
    zB_V =convmtx(zB_V.',L);
    zB_V =flipud(zB_V(:,L:end-L+1));
    if NPol== 2
        % H-pol. orientation:
        zB_H= [zeros(floor(L/2),1);
        z(:,2);
        zeros(floor(L/2),1)];
        zB_H= convmtx(zB_H.',L);
        zB_H = flipud(zB_H(:,L:end-L+1));
        zBlocks = cat(3,zB_V,zB_H);
        clearvars zBlocks_V zBlocks_H;
    else
        zBlocks= zB_V; clearvars zBlocks_V;
    end
    % Vector of phase estimates and initial phase for the PU:
    ThetaPU = zeros(size(zBlocks,2),NPol);
    ThetaPrev = zeros(1,NPol);
    Thetas = zeros(size(zBlocks,2),NPol);
    % Phasenoise estimates:
    for i =1:size(zBlocks,2)
        % Applying the testphase angles to the symbols:
        zRot= repmat(zBlocks(:,i,:),1,B,1).*ThetaTestMatrix;
        % Decision of the rotated symbols:
        zRot_Decided = qam_slicer(zRot,M);
        % Intermidiate signal to be minimized:
        m =sum(abs(zRot-zRot_Decided).^2,1);
        % Estimating the phase noise as the angle that minimizes 'm':
        [~,im] = min(m,[],2);
        Thetas(i,:) = reshape(ThetaTest(im),1,NPol);
    end

    % Applying the phase unwrapper to the estimated phase angles:
    for i = 1:size(zBlocks,2)
        ThetaPU(i,:) = Thetas(i,:) + floor(0.5 - (Thetas(i,:) - ThetaPrev)./(p)).*(p);
        if UsePilots
            BlockIdx = ceil(i / BlockLen);
            n = 0;
            size(ThetaPU);
            size(PhiRef);
            size(zBlocks,2);
            PhaseDiff = ThetaPU(i, :) - PhiRef(BlockIdx, :);
            if max(PhaseDiff) > pi /2
                n = 1;
            elseif min(PhaseDiff) < -pi/2
                n = -1;
            end
            ThetaPU(i, :) = ThetaPU(i, :) - n * pi/2;
        end
        if BlockBased && mod(i, BlockLen) ~= 1
            ThetaPU(i,:) = ThetaPU(i-1,:);
        end
        % Updating the previous phase variable:
        ThetaPrev = ThetaPU(i,:);
    end
    v = z.*exp(-1i*ThetaPU);
end