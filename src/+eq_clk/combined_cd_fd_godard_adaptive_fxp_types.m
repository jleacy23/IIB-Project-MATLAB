function T = combined_cd_fd_godard_adaptive_fxp_types(cfg) %#codegen
%COMBINED_CD_FD_GODARD_ADAPTIVE_FXP_TYPES  Composite type table for the
%   fixed-point CD-FD + embedded Godard PI + adaptive-EQ block.
%
%   T = combined_cd_fd_godard_adaptive_fxp_types(cfg)
%
%   Builds a struct with one sub-table per pipeline section:
%
%     T.Static  - frequency-domain CD + matched filter (overlap-save)
%                 (cd_eq.equalize_fxp_types fields: x, tw, hcd, acc)
%     T.Godard  - in-loop Godard timing recovery (metric, loop filter,
%                 phase-ramp twiddle).
%                   T.Godard.metric  - per-block metric S accumulator
%                   T.Godard.ek      - imag(S), drives the PI loop
%                   T.Godard.lf      - LF_I, tauSamp, ki*e, kp*e
%                   T.Godard.tw      - phase-ramp exp(-j*2*pi*k*tau/N)
%     T.AdaptEq - adaptive butterfly equaliser
%                 (adaptive_eq.equalize_fxp_types fields)
%
%   Note: coarse CFO correction is applied in floating point inside the
%   combined block (cast fi -> double, run eq_clk.coarse_cfo_fd, cast
%   back to T.Static.x), so no T.CFO sub-table is required.
%
%   Calling conventions for cfg:
%     1. Scalar config (string 'fixed16'|'fixed32'|'double'|'single', or a
%        struct('WL',wl,'FL',fl)) — same config forwarded to every section.
%     2. Struct with section sub-fields cfg.Static, cfg.Godard, cfg.AdaptEq.
%        Missing sub-fields default to 'fixed32'.

    if isstruct(cfg) && (isfield(cfg, 'Static') ...
                      || isfield(cfg, 'Godard') || isfield(cfg, 'AdaptEq'))
        T.Static  = cd_eq.equalize_fxp_types(getOr(cfg, 'Static',  'fixed32'));
        T.Godard  = godardTypes(getOr(cfg, 'Godard', 'fixed32'));
        T.AdaptEq = adaptive_eq.equalize_fxp_types(getOr(cfg, 'AdaptEq', 'fixed32'));
        return;
    end

    T.Static  = cd_eq.equalize_fxp_types(cfg);
    T.Godard  = godardTypes(cfg);
    T.AdaptEq = adaptive_eq.equalize_fxp_types(cfg);
end


function G = godardTypes(dt)
%GODARDTYPES  Builds the embedded-Godard sub-table.
%   Fields: metric (complex), ek (real), lf (real), tw (complex twiddle).

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
        G.metric = fi([], 1, wl, fl, F);
        G.ek     = fi([], 1, wl, fl, F);
        G.lf     = wideLfType();   % wide accumulator, independent of swept fl
        G.tw     = fi([], 1, wl, fl, F);
        return;
    end

    switch dt
        case 'double'
            G.metric = double([]);
            G.ek     = double([]);
            G.lf     = double([]);
            G.tw     = double([]);

        case 'single'
            G.metric = single([]);
            G.ek     = single([]);
            G.lf     = single([]);
            G.tw     = single([]);

        case 'fixed16'
            F = fimath( ...
                'RoundingMethod',       'Floor', ...
                'OverflowAction',       'Wrap',  ...
                'ProductMode',          'SpecifyPrecision', ...
                'ProductWordLength',     32, ...
                'ProductFractionLength', 8, ...
                'SumMode',              'SpecifyPrecision', ...
                'SumWordLength',         32, ...
                'SumFractionLength',     8);
            G.metric = fi([], 1, 32, 8,  F);
            G.ek     = fi([], 1, 32, 8,  F);
            G.lf     = wideLfType();
            G.tw     = fi([], 1, 32, 14, F);

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
            G.metric = fi([], 1, 32, 16, F);
            G.ek     = fi([], 1, 32, 16, F);
            G.lf     = wideLfType();
            G.tw     = fi([], 1, 32, 24, F);

        otherwise
            error('combined_cd_fd_godard_adaptive_fxp_types:BadType', ...
                'Unknown Godard type configuration ''%s''.', dt);
    end
end


function v = getOr(s, name, dflt)
    if isfield(s, name)
        v = s.(name);
    else
        v = dflt;
    end
end


function lf = wideLfType()
%WIDELFTYPE  Wide fixed-point Godard loop-filter accumulator (LF_I, tauSamp,
%   ki*e).  Its width is a FIXED design constant — NOT the swept data-path
%   precision — sized so the tiny PI gains (ki ~ 1e-6/1e-4) and their
%   products neither underflow nor lose the integral.  The swept "specified"
%   precision is applied downstream to the FD phase-ramp twiddle (G.tw).
    WL = 48; FL = 40;
    F  = fimath( ...
        'RoundingMethod',       'Floor', ...
        'OverflowAction',       'Wrap',  ...
        'ProductMode',          'SpecifyPrecision', ...
        'ProductWordLength',     WL, ...
        'ProductFractionLength', FL, ...
        'SumMode',              'SpecifyPrecision', ...
        'SumWordLength',         WL, ...
        'SumFractionLength',     FL);
    lf = fi([], 1, WL, FL, F);
end
