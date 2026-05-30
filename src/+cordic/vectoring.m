function theta = vectoring(xr, xi, nIters, T) %#codegen
%VECTORING  Fixed-point CORDIC angle estimate, theta ~ atan2(xi, xr).
%
%   theta = cordic.vectoring(xr, xi, nIters, T)
%
%   Computes the phase angle of the vector (xr, xi) using nIters of the
%   CORDIC vectoring recurrence.  Vectoring mode rotates the vector onto the
%   positive real axis while accumulating the rotation, so the returned
%   angle is scale-invariant (only the ratio xi/xr matters).  The CORDIC
%   processing gain therefore needs no compensation here.
%
%   Inputs
%     xr, xi  - real scalars (fi or double); real / imaginary parts of the
%               vector whose angle is wanted.  Cast internally to T.theta,
%               so the caller must ensure their magnitude fits T.theta.
%     nIters  - number of CORDIC iterations.  Angular resolution is
%               ~atan(2^-nIters); in the bit-width sweep this is set equal to
%               the swept fractional-bit precision.
%     T       - type table (carrier_recovery/freq_recovery.fxp_types).  All
%               arithmetic uses T.theta with its SpecifyPrecision fimath.
%
%   Output
%     theta   - phase angle in (-pi, pi], type T.theta.
%
%   Codegen notes
%     - No lookup-table array: atan(2^-i) is evaluated inline each iteration
%       (its argument is a compile-time-derived double), keeping the function
%       free of variable-size data for MATLAB Coder.
%     - nIters is a coder.Constant at the MEX entry points, so the loop is
%       a fixed, unrollable length.

    ZERO    = cast(0,    'like', T.theta);
    HALF_PI = cast(pi/2, 'like', T.theta);

    x = cast(xr, 'like', T.theta);
    y = cast(xi, 'like', T.theta);
    z = ZERO;

    %% ----------------------------------------------------------------
    %  Argument reduction into the right half-plane (x >= 0) so the
    %  remaining angle is within the CORDIC convergence range
    %  (|angle| <= sum_i atan(2^-i) ~ 99.9 deg).  A +-90 deg pre-rotation
    %  is exact (swap/negate of the components).
    %% ----------------------------------------------------------------
    if x < ZERO
        if y >= ZERO
            t = x;  x = y;   y = -t;  z =  HALF_PI;   % undo a -90 deg rotation
        else
            t = x;  x = -y;  y =  t;  z = -HALF_PI;   % undo a +90 deg rotation
        end
    end

    %% ----------------------------------------------------------------
    %  CORDIC vectoring recurrence: drive y -> 0, accumulate angle in z.
    %% ----------------------------------------------------------------
    for i = 0:nIters-1
        p2  = cast(2^(-i),      'like', T.theta);
        ang = cast(atan(2^(-i)), 'like', T.theta);
        xs  = cast(x * p2, 'like', T.theta);
        ys  = cast(y * p2, 'like', T.theta);
        if y >= ZERO
            x = x + ys;
            y = y - xs;
            z = z + ang;
        else
            x = x - ys;
            y = y + xs;
            z = z - ang;
        end
    end

    theta = z;
end
