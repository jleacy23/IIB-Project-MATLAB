function T = equalize_fxp_types(dt) %#codegen
%equalize_FXP_TYPES  Data-type table for equalize_fxp.
%
%   T = equalize_fxp_types(dt)
%
%   Returns a struct of fi prototype objects (empty values) that define
%   every fixed-point type used inside equalize_fxp.
%
%   Supported configurations:
%     'double'   - all types are double (floating-point baseline)
%     'single'   - all types are single (useful for mismatch checking)
%     'fixed16'  - 16-bit fixed-point, suitable for FPGA / ASIC
%     'fixed32'  - 32-bit fixed-point, higher precision
%
%   You can add your own cases or adjust word / fraction lengths to
%   explore design trade-offs without modifying the algorithm.
%
%   Fields returned
%     T.x      - input signal
%     T.w      - filter (tap) coefficients
%     T.y      - equalizer output samples
%     T.acc    - accumulator for inner-product computation
%     T.err    - error signal  (R - |y|^2)  or  (R^2 - |y|^2)
%     T.mu     - step-size scalar
%     T.R_CMA  - CMA target radius
%     T.R_RDE  - RDE target radii

    % fimath is defined per-configuration below with SpecifyPrecision
    % for both products and sums.  Every arithmetic result is truncated to
    % the same word-length and fraction-length — no bit-growth, no
    % rescaling, matching a real fixed-point datapath.

    switch dt
        % ==============================================================
        case 'double'
            T.x     = double([]);
            T.w     = double([]);
            T.y     = double([]);
            T.acc   = double([]);
            T.err   = double([]);
            T.mu    = double([]);
            T.R_CMA = double([]);
            T.R_RDE = double([]);

        % ==============================================================
        case 'single'
            T.x     = single([]);
            T.w     = single([]);
            T.y     = single([]);
            T.acc   = single([]);
            T.err   = single([]);
            T.mu    = single([]);
            T.R_CMA = single([]);
            T.R_RDE = single([]);

        % ==============================================================
        case 'fixed16'
            %  Uniform 16-bit / FL=8 throughout.
            %  Range ±4, LSB = 2^{-8} ≈ 3.9e-3.
            %  mu = 1e-3 ≈ 8 LSBs.  Adequate for CMA/RDE convergence.
            F = fimath( ...
                'RoundingMethod',       'Floor', ...
                'OverflowAction',       'Wrap',  ...
                'ProductMode',          'SpecifyPrecision', ...
                'ProductWordLength',     32, ...
                'ProductFractionLength', 8, ...
                'SumMode',              'SpecifyPrecision', ...
                'SumWordLength',         32, ...
                'SumFractionLength',     8);

            T.x     = fi([], 1, 32, 8, F);   % all types identical
            T.w     = fi([], 1, 32, 8, F);
            T.y     = fi([], 1, 32, 8, F);
            T.acc   = fi([], 1, 32, 8, F);
            T.err   = fi([], 1, 32, 8, F);
            T.mu    = fi([], 1, 32, 8, F);
            T.R_CMA = fi([], 1, 32, 8, F);
            T.R_RDE = fi([], 1, 32, 8, F);

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
            T.w     = fi([], 1, 32, 16, F);
            T.y     = fi([], 1, 32, 16, F);
            T.acc   = fi([], 1, 32, 16, F);
            T.err   = fi([], 1, 32, 16, F);
            T.mu    = fi([], 1, 32, 16, F);
            T.R_CMA = fi([], 1, 32, 16, F);
            T.R_RDE = fi([], 1, 32, 16, F);

        otherwise
            error('equalize_fxp_types:BadType', ...
                'Unknown type configuration ''%s''.', dt);
    end
end
