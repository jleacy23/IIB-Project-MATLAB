function T = fft_fxp_types(dt) %#codegen
%FFT_FXP_TYPES  Data-type table for fft_fxp.
%
%   T = fft_fxp_types(dt)
%
%   Returns a struct of fi prototype objects that define every fixed-point
%   type used inside fft_fxp.
%
%   Supported configurations:
%     'double'              - all types are double (floating-point baseline)
%     'single'              - all types are single
%     'fixed16'             - 16-bit fixed-point, uniform WL/FL
%     'fixed32'             - 32-bit fixed-point, uniform WL/FL
%     struct('WL',wl,'FL',fl) - custom: uniform word length wl, fraction length fl
%
%   Fields returned
%     T.x    - input signal
%     T.tw   - twiddle factor
%     T.acc  - accumulator / butterfly output
%
%   The fimath attached to every fi prototype uses SpecifyPrecision for
%   both products and sums so that no bit-growth occurs.  This matches a
%   uniform fixed-point datapath (FPGA / ASIC) where the accumulator word
%   length is held constant across all butterfly stages.

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
            %  Uniform 16-bit / FL=8.
            %  Range ±128, LSB = 2^{-8} ≈ 3.9e-3.
            %  Note: fft_fxp now applies 1/2 inter-stage scaling on the
            %  forward transform, so the magnitude stays bounded (no
            %  wrap-around), but with only 8 fractional bits each stage's
            %  right shift loses LSBs — prefer more fractional bits for
            %  accuracy-critical transforms.
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
            %  Suitable for FFT sizes up to ~1024 with unit-amplitude
            %  input without overflow.
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
            error('fft_fxp_types:BadType', ...
                'Unknown type configuration ''%s''.', dt);
    end
end
