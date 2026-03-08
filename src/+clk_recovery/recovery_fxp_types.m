function T = recovery_fxp_types(dt) %#codegen
%recovery_FXP_TYPES  Data-type table for recovery_fxp.
%
%   T = recovery_fxp_types(dt)
%
%   Returns a struct of fi prototype objects that define every fixed-point
%   type used inside recovery_fxp.
%
%   Supported configurations:
%     'double'   - all types are double (floating-point baseline)
%     'single'   - all types are single
%     'fixed16'  - 16-bit fixed-point
%     'fixed32'  - 32-bit fixed-point
%
%   Fields returned
%     T.x     - input / output signal samples
%     T.acc   - accumulator (interpolator arithmetic, TED, loop filter)
%     T.mu    - fractional interval (mun) and NCO state (Etamn, Wk)
%     T.coeff - interpolator polynomial coefficients (1/6, 1/2, 1/3, …)

    switch dt
        % ==============================================================
        case 'double'
            T.x     = double([]);
            T.acc   = double([]);
            T.mu    = double([]);
            T.coeff = double([]);

        % ==============================================================
        case 'single'
            T.x     = single([]);
            T.acc   = single([]);
            T.mu    = single([]);
            T.coeff = single([]);

        % ==============================================================
        case 'fixed16'
            %  Uniform 16-bit / FL=12.
            %  Range ±8, LSB = 2^{-12} ≈ 2.4e-4.
            %  Higher fraction length than other modules because the NCO
            %  state (Etamn, mun) is in [0,1) and needs fine resolution.
            F = fimath( ...
                'RoundingMethod',       'Floor', ...
                'OverflowAction',       'Wrap',  ...
                'ProductMode',          'SpecifyPrecision', ...
                'ProductWordLength',     16, ...
                'ProductFractionLength', 12, ...
                'SumMode',              'SpecifyPrecision', ...
                'SumWordLength',         16, ...
                'SumFractionLength',     12);

            T.x     = fi([], 1, 16, 12, F);
            T.acc   = fi([], 1, 16, 12, F);
            T.mu    = fi([], 1, 16, 12, F);
            T.coeff = fi([], 1, 16, 12, F);

        % ==============================================================
        case 'fixed32'
            %  Uniform 32-bit / FL=24.
            %  Range ±128, LSB = 2^{-24} ≈ 6.0e-8.
            %  Generous fractional resolution for the NCO and loop filter.
            F = fimath( ...
                'RoundingMethod',       'Floor', ...
                'OverflowAction',       'Wrap',  ...
                'ProductMode',          'SpecifyPrecision', ...
                'ProductWordLength',     32, ...
                'ProductFractionLength', 24, ...
                'SumMode',              'SpecifyPrecision', ...
                'SumWordLength',         32, ...
                'SumFractionLength',     24);

            T.x     = fi([], 1, 32, 24, F);
            T.acc   = fi([], 1, 32, 24, F);
            T.mu    = fi([], 1, 32, 24, F);
            T.coeff = fi([], 1, 32, 24, F);

        otherwise
            error('recovery_fxp_types:BadType', ...
                'Unknown type configuration ''%s''.', dt);
    end
end
