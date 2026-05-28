function T = combined_cd_fd_gardner_adaptive_fxp_types(cfg) %#codegen
%COMBINED_CD_FD_GARDNER_ADAPTIVE_FXP_TYPES  Composite type table for the
%   fixed-point combined CD-FD + Gardner + adaptive-EQ block.
%
%   T = combined_cd_fd_gardner_adaptive_fxp_types(cfg)
%
%   Builds a struct with one sub-table per pipeline section, so each
%   section can have its own bit widths:
%
%     T.Static  - frequency-domain CD + matched filter
%                 (= cd_eq.equalize_fxp_types fields: x, tw, hcd, acc)
%     T.Clk     - Gardner DPLL clock recovery
%                 (= clk_recovery.recovery_fxp_types fields)
%     T.AdaptEq - adaptive butterfly equaliser
%                 (= adaptive_eq.equalize_fxp_types fields)
%
%   Note: coarse CFO correction is applied in floating point inside the
%   combined block (cast fi -> double, run eq_clk.coarse_cfo_fd, cast
%   back to T.Static.x), so no T.CFO sub-table is required.
%
%   Calling conventions for cfg:
%     1. Scalar config (string 'fixed16'|'fixed32'|'double'|'single', or a
%        struct('WL',wl,'FL',fl)).  The same config is forwarded to every
%        sub-table.
%     2. Struct with section sub-fields cfg.Static, cfg.Clk, cfg.AdaptEq.
%        Each sub-field is forwarded to the corresponding sub-table.
%        Missing sub-fields default to 'fixed32'.

    if isstruct(cfg) && (isfield(cfg, 'Static') ...
                      || isfield(cfg, 'Clk')   || isfield(cfg, 'AdaptEq'))
        % Per-section configs
        T.Static  = cd_eq.equalize_fxp_types(getOr(cfg, 'Static',  'fixed32'));
        T.Clk     = clk_recovery.recovery_fxp_types(getOr(cfg, 'Clk',     'fixed32'));
        T.AdaptEq = adaptive_eq.equalize_fxp_types(getOr(cfg, 'AdaptEq', 'fixed32'));
        return;
    end

    %% Single scalar config — forward to every section
    T.Static  = cd_eq.equalize_fxp_types(cfg);
    T.Clk     = clk_recovery.recovery_fxp_types(cfg);
    T.AdaptEq = adaptive_eq.equalize_fxp_types(cfg);
end


function v = getOr(s, name, dflt)
    if isfield(s, name)
        v = s.(name);
    else
        v = dflt;
    end
end
