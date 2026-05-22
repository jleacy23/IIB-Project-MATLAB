function T = equalize_td_fxp_types(dt) %#codegen
%EQUALIZE_TD_FXP_TYPES  Data-type table for equalize_td_fxp.
%
%   T = equalize_td_fxp_types(dt)
%
%   Returns a struct of fi prototype objects that define every fixed-point
%   type used inside equalize_td_fxp (time-domain FIR CD compensation).
%
%   Supported configurations:
%     'double'              - all types are double (floating-point baseline)
%     'single'              - all types are single
%     'fixed16'             - 16-bit fixed-point, uniform WL/FL
%     'fixed32'             - 32-bit fixed-point, uniform WL/FL
%     struct('WL',wl,'FL',fl) - custom: uniform word length wl, fraction length fl
%
%   Fields returned
%     T.x    - input / output signal
%     T.hcd  - CD impulse-response (FIR tap) coefficients
%     T.acc  - accumulator (FIR multiply-accumulate)
%
%   The fimath attached to every fi prototype uses SpecifyPrecision for
%   both products and sums so that no bit-growth occurs — matching a
%   uniform fixed-point datapath (FPGA / ASIC).
%
%   This mirrors equalize_fxp_types but omits the FFT twiddle field (T.tw),
%   which the time-domain implementation does not use.

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
        T.x   = fi([], 1, wl, fl, F);
        T.hcd = fi([], 1, wl, fl, F);
        T.acc = fi([], 1, wl, fl, F);
        return;
    end

    switch dt
        % ==============================================================
        case 'double'
            T.x   = double([]);
            T.hcd = double([]);
            T.acc = double([]);

        % ==============================================================
        case 'single'
            T.x   = single([]);
            T.hcd = single([]);
            T.acc = single([]);

        % ==============================================================
        case 'fixed16'
            %  Uniform 16-bit / FL=8.
            %  Range ±128, LSB = 2^{-8} ≈ 3.9e-3.
            F = fimath( ...
                'RoundingMethod',       'Floor', ...
                'OverflowAction',       'Wrap',  ...
                'ProductMode',          'SpecifyPrecision', ...
                'ProductWordLength',     32, ...
                'ProductFractionLength', 8, ...
                'SumMode',              'SpecifyPrecision', ...
                'SumWordLength',         32, ...
                'SumFractionLength',     8);

            T.x   = fi([], 1, 32, 8, F);
            T.hcd = fi([], 1, 32, 8, F);
            T.acc = fi([], 1, 32, 8, F);

        % ==============================================================
        case 'fixed32'
            %  Uniform 32-bit / FL=16.
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

            T.x   = fi([], 1, 32, 16, F);
            T.hcd = fi([], 1, 32, 16, F);
            T.acc = fi([], 1, 32, 16, F);

        otherwise
            error('equalize_td_fxp_types:BadType', ...
                'Unknown type configuration ''%s''.', dt);
    end
end
