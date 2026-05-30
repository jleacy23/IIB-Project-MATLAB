function [yr, yi] = rotate(xr, xi, theta, nIters, T, Vproto) %#codegen
%ROTATE  Fixed-point CORDIC phase rotation: (yr + j*yi) = (xr + j*xi)*e^{j*theta}.
%
%   [yr, yi] = cordic.rotate(xr, xi, theta, nIters, T)
%   [yr, yi] = cordic.rotate(xr, xi, theta, nIters, T, Vproto)
%
%   Rotates the vector (xr, xi) by +theta using nIters of the CORDIC
%   rotation recurrence.  To de-rotate a symbol by an estimated phase
%   ThetaPU, call with theta = -ThetaPU.  The CORDIC processing gain is
%   compensated by pre-scaling the input by 1/K(nIters), so the output
%   magnitude matches the input (apart from fixed-point quantisation).
%
%   Inputs
%     xr, xi  - real scalars (fi or double); vector to rotate.
%     theta   - rotation angle [rad] (fi or double scalar).
%     nIters  - number of CORDIC iterations (angular resolution ~atan(2^-nIters)).
%     T       - type table; the angle accumulator uses T.theta.
%     Vproto  - (optional) fi prototype for the rotated vector datapath.
%               Defaults to T.x.  Pass T.acc to synthesise / rotate vectors in
%               the wide accumulator domain (e.g. building a unit-circle FFT
%               input) rather than the narrow signal type.
%
%   Outputs
%     yr, yi  - rotated real / imaginary parts, type Vproto (default T.x).
%
%   Codegen notes
%     - Angle is reduced to [-pi/2, pi/2] (CORDIC convergence range) by an
%       optional +-pi fold; a pi fold is realised by negating the result.
%     - The processing gain K = prod_i sqrt(1 + 2^-2i) depends on nIters and
%       is accumulated in a scalar loop (nIters is a coder.Constant at the
%       MEX entry points), so no precomputed table is needed.

    if nargin < 6 || isempty(Vproto)
        Vproto = T.x;
    end

    ZERO_TH = cast(0, 'like', T.theta);

    %% ----------------------------------------------------------------
    %  Angle reduction to [-pi/2, pi/2].  Fold by +-pi where needed and
    %  flag a vector negation (rotation by pi == negation).
    %% ----------------------------------------------------------------
    phi = mod(double(theta) + pi, 2*pi) - pi;   % wrap to (-pi, pi]
    neg = false;
    if phi > pi/2
        phi = phi - pi;  neg = true;
    elseif phi < -pi/2
        phi = phi + pi;  neg = true;
    end

    %% ----------------------------------------------------------------
    %  CORDIC processing gain and input pre-scaling (gain compensation).
    %% ----------------------------------------------------------------
    Kinv = 1.0;
    for i = 0:nIters-1
        Kinv = Kinv / sqrt(1 + 2^(-2*i));
    end

    x = cast(double(xr) * Kinv, 'like', Vproto);
    y = cast(double(xi) * Kinv, 'like', Vproto);
    z = cast(phi, 'like', T.theta);

    %% ----------------------------------------------------------------
    %  CORDIC rotation recurrence: drive z -> 0, applying the same rotation
    %  to (x, y) so the vector ends up rotated by the original phi.
    %% ----------------------------------------------------------------
    for i = 0:nIters-1
        p2x = cast(2^(-i),      'like', Vproto);
        ang = cast(atan(2^(-i)), 'like', T.theta);
        xs  = cast(x * p2x, 'like', Vproto);
        ys  = cast(y * p2x, 'like', Vproto);
        if z >= ZERO_TH
            x = x - ys;
            y = y + xs;
            z = z - ang;
        else
            x = x + ys;
            y = y - xs;
            z = z + ang;
        end
    end

    if neg
        x = -x;
        y = -y;
    end

    yr = x;
    yi = y;
end
