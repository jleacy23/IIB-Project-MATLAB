function T = recovery_fxp_types(dt) %#codegen
%RECOVERY_FXP_TYPES  Data-type table for clk_recovery.recovery_fxp (Gardner DPLL).
%
%   T = recovery_fxp_types(dt)
%
%   Returns a struct of fi prototype objects that define every fixed-point
%   type used inside clk_recovery.recovery_fxp.
%
%   Supported configurations:
%     'double'              - all types are double (floating-point baseline)
%     'single'              - all types are single
%     'fixed16'             - 16-bit fixed-point, uniform WL/FL
%     'fixed32'             - 32-bit fixed-point, uniform WL/FL
%     struct('WL',wl,'FL',fl) - custom: uniform word length wl, fraction length fl
%
%   Fields returned
%     T.x     - signal data path (input, interpolator output, Out buffer)
%     T.coef  - cubic Farrow interpolator coefficients (real, in [-1, 1])
%     T.acc   - interpolator MAC accumulator (complex)
%     T.ek    - per-block Gardner timing error accumulator (real)
%     T.lf    - loop filter state (Wk, LF_I) and gain products (ki*e, kp*e)
%     T.nco   - NCO fractional state (Etamn, mun)  in [0, 1)

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
        T.x    = fi([], 1, wl, fl, F);
        T.coef = fi([], 1, wl, fl, F);
        T.acc  = fi([], 1, wl, fl, F);
        T.ek   = fi([], 1, wl, fl, F);
        T.lf   = fi([], 1, wl, fl, F);
        T.nco  = fi([], 1, wl, fl, F);
        return;
    end

    switch dt
        case 'double'
            T.x    = double([]);
            T.coef = double([]);
            T.acc  = double([]);
            T.ek   = double([]);
            T.lf   = double([]);
            T.nco  = double([]);

        case 'single'
            T.x    = single([]);
            T.coef = single([]);
            T.acc  = single([]);
            T.ek   = single([]);
            T.lf   = single([]);
            T.nco  = single([]);

        case 'fixed16'
            %  Uniform 16-bit / FL=8.  NCO fractional state has FL=14
            %  (it lives in [0, 1)) and the loop filter state has FL=20
            %  (Wk ≈ 1; LF_I drifts only by ki * ek per block).
            F = fimath( ...
                'RoundingMethod',       'Floor', ...
                'OverflowAction',       'Wrap',  ...
                'ProductMode',          'SpecifyPrecision', ...
                'ProductWordLength',     32, ...
                'ProductFractionLength', 8, ...
                'SumMode',              'SpecifyPrecision', ...
                'SumWordLength',         32, ...
                'SumFractionLength',     8);
            T.x    = fi([], 1, 32, 8,  F);
            T.coef = fi([], 1, 32, 16, F);
            T.acc  = fi([], 1, 32, 8,  F);
            T.ek   = fi([], 1, 32, 8,  F);
            T.lf   = fi([], 1, 32, 20, F);
            T.nco  = fi([], 1, 32, 14, F);

        case 'fixed32'
            F = fimath( ...
                'RoundingMethod',       'Floor', ...
                'OverflowAction',       'Wrap',  ...
                'ProductMode',          'SpecifyPrecision', ...
                'ProductWordLength',     32, ...
                'ProductFractionLength', 16, ...
                'SumMode',              'SpecifyPrecision', ...
                'SumWordLength',         32, ...
                'SumFractionLength',     16);
            T.x    = fi([], 1, 32, 16, F);
            T.coef = fi([], 1, 32, 24, F);
            T.acc  = fi([], 1, 32, 16, F);
            T.ek   = fi([], 1, 32, 16, F);
            T.lf   = fi([], 1, 32, 24, F);
            T.nco  = fi([], 1, 32, 24, F);

        otherwise
            error('recovery_fxp_types:BadType', ...
                'Unknown type configuration ''%s''.', dt);
    end
end
