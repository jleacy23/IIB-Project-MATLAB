function T = fxp_types(dt) %#codegen
%FXP_TYPES  Data-type table for freq_recovery fixed-point algorithms.
%
%   T = fxp_types(dt)
%
%   Returns a struct of fi prototype objects that define every fixed-point
%   type used inside tretter_kay_fxp, fitz_fxp, and fft_search_fxp.
%
%   Supported configurations
%     'double'   - all types are double  (floating-point baseline)
%     'single'   - all types are single
%     'fixed16'  - 16-bit word length, 8-bit fraction
%     'fixed32'  - 32-bit word length, 16-bit fraction
%
%   Fields returned
%     T.x      - input / output signal prototype
%     T.theta  - phase angle prototype  (must accommodate ±pi)
%     T.acc    - accumulator for weighted sums and autocorrelation

    switch dt

        % ==============================================================
        case 'double'
            T.x     = double([]);
            T.theta = double([]);
            T.acc   = double([]);

        % ==============================================================
        case 'single'
            T.x     = single([]);
            T.theta = single([]);
            T.acc   = single([]);

        % ==============================================================
        case 'fixed16'
            % 16-bit WL, FL = 8.  Range ±128, LSB = 2^{-8} ≈ 3.9e-3.
            % pi < 128 — angle type is safe.
            %
            % T.x / T.theta: sums and products stay at 16 bits so that
            %   variables declared 'like' these types remain the same type
            %   after arithmetic (needed for codegen loop variables).
            %   Phase wrap is natural with OverflowAction='Wrap'.
            %
            % T.acc: 32-bit accumulator for weighted sums; product also
            %   kept at 32 bits to avoid overflow during multiply-accumulate.
            F16 = fimath( ...
                'RoundingMethod',       'Floor', ...
                'OverflowAction',       'Wrap',  ...
                'ProductMode',          'SpecifyPrecision', ...
                'ProductWordLength',     16, ...
                'ProductFractionLength', 8, ...
                'SumMode',              'SpecifyPrecision', ...
                'SumWordLength',         16, ...
                'SumFractionLength',     8);

            F32acc = fimath( ...
                'RoundingMethod',       'Floor', ...
                'OverflowAction',       'Wrap',  ...
                'ProductMode',          'SpecifyPrecision', ...
                'ProductWordLength',     32, ...
                'ProductFractionLength', 8, ...
                'SumMode',              'SpecifyPrecision', ...
                'SumWordLength',         32, ...
                'SumFractionLength',     8);

            T.x     = fi([], 1, 16, 8, F16);
            T.theta = fi([], 1, 16, 8, F16);
            T.acc   = fi([], 1, 32, 8, F32acc);

        % ==============================================================
        case 'fixed32'
            % 32-bit WL, FL = 16.  Range ±32768, LSB = 2^{-16} ≈ 1.5e-5.
            %
            % T.x / T.theta: sums and products at 32 bits (native width).
            % T.acc: also 32 bits — at this precision T.acc == T.x/T.theta.
            F32 = fimath( ...
                'RoundingMethod',       'Floor', ...
                'OverflowAction',       'Wrap',  ...
                'ProductMode',          'SpecifyPrecision', ...
                'ProductWordLength',     32, ...
                'ProductFractionLength', 16, ...
                'SumMode',              'SpecifyPrecision', ...
                'SumWordLength',         32, ...
                'SumFractionLength',     16);

            T.x     = fi([], 1, 32, 16, F32);
            T.theta = fi([], 1, 32, 16, F32);
            T.acc   = fi([], 1, 32, 16, F32);

        otherwise
            error('freq_recovery.fxp_types:BadType', ...
                'Unknown type configuration ''%s''.', dt);
    end
end
