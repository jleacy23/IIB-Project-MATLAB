function T = equalize_fxp_types(dt) %#codegen
%EQUALIZE_FXP_TYPES  Fixed-point type table for equalize_fxp (adaptive EQ).
%
%   T = equalize_fxp_types(dt)
%
%   Returns a struct of fi prototype objects that define every fixed-point
%   type used inside equalize_fxp.
%
%   Supported configurations
%     'double'              - all types are double (floating-point baseline / reference)
%     'single'              - all types are single
%     'fixed16'             - uniform 16-bit fixed-point, suitable for FPGA / ASIC
%     'fixed32'             - uniform 32-bit fixed-point, higher precision
%     struct('WL',wl,'FL',fl) - custom: uniform word length wl, fraction length fl
%
%   Fields
%     T.x      - input signal (complex QAM samples, tap delay line)
%     T.w      - filter / tap coefficients (complex, updated by CMA)
%     T.y      - equalizer output samples
%     T.acc    - accumulator for the butterfly inner product
%     T.err    - error signal  (R_CMA - |y|^2)
%     T.grad   - unscaled gradient  x * err * conj(y)  before mu scaling
%     T.R_CMA  - CMA target radius
%
%   Weight update strategy
%     The gradient x*err*conj(y) is computed in T.grad fixed-point precision,
%     then converted to double and scaled by mu (a plain double scalar) before
%     being cast back to T.w.  This avoids mu being rounded to zero when it is
%     too small to be represented in the weight bit-width.  T.w must therefore
%     have enough fractional bits to represent mu*grad after the scaling.
%
%   Numerical design notes
%
%   Signal range
%     Unit-average-power QAM samples have |z| ≈ 1; with a CMA target of
%     R_CMA = sqrt(2) the equalised samples sit on |y| ≈ 1.  Tap weights
%     and accumulator stay within a few units, so a few integer bits suffice.
%
%   Gradient range
%     grad = x * err * conj(y) has |grad| < 2 once convergence is approached.
%     T.grad should have enough fractional bits to represent this faithfully
%     before the mu scaling step promotes small updates out of fixed-point.
%
%   Update precision
%     fixed16 (FL = 8, LSB ≈ 3.9e-3) is comfortable for CMA convergence.
%     fixed32 (FL = 16, LSB ≈ 1.5e-5) gives near-floating-point behaviour.
%
%   SpecifyPrecision fimath
%     All fi arithmetic uses SpecifyPrecision so every product and sum is
%     truncated to a known WL/FL with no implicit bit growth — required for
%     deterministic codegen behaviour and to mimic a uniform fixed-point
%     datapath (FPGA / ASIC).

    if isstruct(dt)
        wl = dt.WL;
        fl = dt.FL;
        F = fimath( ...
            'RoundingMethod',       'Floor', ...
            'OverflowAction',       'Wrap',  ...
            'ProductMode',          'SpecifyPrecision', ...
            'ProductWordLength',     wl, ...
            'ProductFractionLength', fl, ...
            'SumMode',              'SpecifyPrecision', ...
            'SumWordLength',         wl, ...
            'SumFractionLength',     fl);
        T.x     = fi([], 1, wl, fl, F);
        T.w     = fi([], 1, wl, fl, F);
        T.y     = fi([], 1, wl, fl, F);
        T.acc   = fi([], 1, wl, fl, F);
        T.err   = fi([], 1, wl, fl, F);
        T.grad  = fi([], 1, wl, fl, F);
        T.R_CMA = fi([], 1, wl, fl, F);
        return;
    end

    switch dt
        % ==============================================================
        case 'double'
            T.x     = double([]);
            T.w     = double([]);
            T.y     = double([]);
            T.acc   = double([]);
            T.err   = double([]);
            T.grad  = double([]);
            T.R_CMA = double([]);

        % ==============================================================
        case 'single'
            T.x     = single([]);
            T.w     = single([]);
            T.y     = single([]);
            T.acc   = single([]);
            T.err   = single([]);
            T.grad  = single([]);
            T.R_CMA = single([]);

        % ==============================================================
        case 'fixed16'
            %  Uniform 16-bit / FL=8 throughout.
            %  Range ±4, LSB = 2^{-8} ≈ 3.9e-3.
            %  mu = 1e-3 ≈ 8 LSBs.  Adequate for CMA convergence.
            F = fimath( ...
                'RoundingMethod',       'Floor', ...
                'OverflowAction',       'Wrap',  ...
                'ProductMode',          'SpecifyPrecision', ...
                'ProductWordLength',     32, ...
                'ProductFractionLength', 8, ...
                'SumMode',              'SpecifyPrecision', ...
                'SumWordLength',         32, ...
                'SumFractionLength',     8);

            T.x     = fi([], 1, 32, 8,  F);
            T.w     = fi([], 1, 32, 24, F);   % extra FL lets mu*grad be non-zero
            T.y     = fi([], 1, 32, 8,  F);
            T.acc   = fi([], 1, 32, 8,  F);
            T.err   = fi([], 1, 32, 8,  F);
            T.grad  = fi([], 1, 32, 8,  F);   % unscaled gradient type (same as acc)
            T.R_CMA = fi([], 1, 32, 8,  F);

        % ==============================================================
        case 'fixed32'
            %  Uniform 32-bit / FL=16 throughout.
            %  Range ±8, LSB = 2^{-16} ≈ 3.7e-9.
            F = fimath( ...
                'RoundingMethod',       'Floor', ...
                'OverflowAction',       'Wrap',  ...
                'ProductMode',          'SpecifyPrecision', ...
                'ProductWordLength',     32, ...
                'ProductFractionLength', 16, ...
                'SumMode',              'SpecifyPrecision', ...
                'SumWordLength',         32, ...
                'SumFractionLength',     16);

            T.x     = fi([], 1, 32, 16, F);
            T.w     = fi([], 1, 32, 28, F);   % extra FL lets mu*grad be non-zero
            T.y     = fi([], 1, 32, 16, F);
            T.acc   = fi([], 1, 32, 16, F);
            T.err   = fi([], 1, 32, 16, F);
            T.grad  = fi([], 1, 32, 16, F);   % unscaled gradient type (same as acc)
            T.R_CMA = fi([], 1, 32, 16, F);

        otherwise
            error('equalize_fxp_types:BadType', ...
                'Unknown type configuration ''%s''.', dt);
    end
end
