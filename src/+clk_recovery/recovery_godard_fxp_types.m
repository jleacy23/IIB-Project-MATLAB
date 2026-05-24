function T = recovery_godard_fxp_types(dt) %#codegen
%recovery_godard_fxp_types  Data-type table for recovery_godard_fxp.
%
%   T = recovery_godard_fxp_types(dt)
%
%   Returns a struct of fi prototype objects that define every fixed-point
%   type used inside recovery_godard_fxp.
%
%   Supported configurations:
%     'double'              - all types are double (floating-point baseline)
%     'single'              - all types are single
%     'fixed16'             - 16-bit fixed-point datapath
%     'fixed32'             - 32-bit fixed-point datapath
%     struct('WL',wl,'FL',fl) - custom: uniform word length wl, fraction length fl
%
%   Fields returned
%     T.x    - input / output signal
%     T.tw   - FFT twiddle factors AND phase-ramp sin/cos values
%     T.acc  - accumulator (FFT butterflies, MG sum, freq-domain multiply)
%
%   The fimath attached to every fi prototype uses SpecifyPrecision for
%   both products and sums so that no bit-growth occurs — matching a
%   uniform fixed-point datapath (FPGA / ASIC).

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
        T.tw  = fi([], 1, wl, fl, F);
        T.acc = fi([], 1, wl, fl, F);
        return;
    end

    switch dt
        % ==============================================================
        case 'double'
            T.x   = double([]);
            T.tw  = double([]);
            T.acc = double([]);

        % ==============================================================
        case 'single'
            T.x   = single([]);
            T.tw  = single([]);
            T.acc = single([]);

        % ==============================================================
        case 'fixed16'
            %  Uniform 16-bit / FL=8 (32-bit product/sum container).
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
            T.tw  = fi([], 1, 32, 8, F);
            T.acc = fi([], 1, 32, 8, F);

        % ==============================================================
        case 'fixed32'
            %  Uniform 32-bit / FL=16.
            %  Range ±32768, LSB = 2^{-16} ≈ 1.5e-5.
            %  Sufficient for FFT sizes up to ~1024 with unit-amplitude
            %  inputs without overflow in the FFT or the MG sum.
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
            T.tw  = fi([], 1, 32, 16, F);
            T.acc = fi([], 1, 32, 16, F);

        otherwise
            error('recovery_godard_fxp_types:BadType', ...
                'Unknown type configuration ''%s''.', dt);
    end
end
