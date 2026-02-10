function [w1V, w1H, w2V, w2H] = adeq_rde_update(xV, xH, y1, y2, w1V, w1H, w2V, w2H, R, Mu)
%ADEQ_RDE_UPDATE  RDE weight update for butterfly adaptive equalizer.
%
%   [w1V,w1H,w2V,w2H] = adeq_rde_update(xV,xH,y1,y2,w1V,w1H,w2V,w2H,R,Mu)

    [~, r1] = min(abs(R - abs(y1)));
    [~, r2] = min(abs(R - abs(y2)));

    w1V = w1V + Mu*xV*(R(r1)^2 - abs(y1).^2)*conj(y1);
    w1H = w1H + Mu*xH*(R(r1)^2 - abs(y1).^2)*conj(y1);
    w2V = w2V + Mu*xV*(R(r2)^2 - abs(y2).^2)*conj(y2);
    w2H = w2H + Mu*xH*(R(r2)^2 - abs(y2).^2)*conj(y2);
end
