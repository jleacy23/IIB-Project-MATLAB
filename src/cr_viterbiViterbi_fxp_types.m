function T = cr_viterbiViterbi_fxp_types(dt) %#codegen
%CR_VITERBIVITERBI_FXP_TYPES  Data-type table for cr_viterbiViterbi_fxp.
%
%   T = cr_viterbiViterbi_fxp_types(dt)
%
%   Returns a struct of fi prototype objects (empty values) that define
%   every fixed-point type used inside cr_viterbiViterbi_fxp.
%
%   Supported configurations:
%     'double'   - all types are double (floating-point baseline)
%     'single'   - all types are single
%     'fixed16'  - 16-bit fixed-point, uniform WL/FL
%     'fixed32'  - 32-bit fixed-point, uniform WL/FL
%
%   Fields returned
%     T.x      - input / output signal
%     T.w      - VV filter coefficients
%     T.theta  - phase angle (range ±pi)
%     T.acc    - accumulator for inner products / sums

    switch dt
        % ==============================================================
        case 'double'
            T.x     = double([]);
            T.w     = double([]);
            T.theta = double([]);
            T.acc   = double([]);

        % ==============================================================
        case 'single'
            T.x     = single([]);
            T.w     = single([]);
            T.theta = single([]);
            T.acc   = single([]);

        % ==============================================================
        case 'fixed16'
            %  Uniform 16-bit / FL=8 throughout.
            %  Range ±128, LSB = 2^{-8} ≈ 3.9e-3.
            %  pi ≈ 3.14 fits well within ±128.
            F = fimath( ...
                'RoundingMethod',       'Floor', ...
                'OverflowAction',       'Wrap',  ...
                'ProductMode',          'SpecifyPrecision', ...
                'ProductWordLength',     32, ...
                'ProductFractionLength', 8, ...
                'SumMode',              'SpecifyPrecision', ...
                'SumWordLength',         32, ...
                'SumFractionLength',     8);

            T.x     = fi([], 1, 32, 8, F);
            T.w     = fi([], 1, 32, 8, F);
            T.theta = fi([], 1, 32, 8, F);
            T.acc   = fi([], 1, 32, 8, F);

        % ==============================================================
        case 'fixed32'
            %  Uniform 32-bit / FL=16 throughout.
            %  Range ±32768, LSB = 2^{-16} ≈ 1.5e-5.
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
            T.w     = fi([], 1, 32, 16, F);
            T.theta = fi([], 1, 32, 16, F);
            T.acc   = fi([], 1, 32, 16, F);

        otherwise
            error('cr_viterbiViterbi_fxp_types:BadType', ...
                'Unknown type configuration ''%s''.', dt);
    end
end
