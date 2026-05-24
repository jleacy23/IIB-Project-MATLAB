function build_cd_eq_equalize_td_fxp_mex(P, cfg)
%BUILD_CDEQ_EQUALIZE_TD_FXP_MEX  Compile cd_eq.equalize_td_fxp to MEX.
%
%   build_cd_eq_equalize_td_fxp_mex(P, cfg)
%
%   Inputs
%     P   - parameter struct from pipeline_params()
%     cfg - coder.MexCodeConfig object
%
%   Required fields in P
%     P.FxpConfig_CD - fixed-point config: 'fixed16' | 'fixed32' | struct('WL',wl,'FL',fl)

    srcDir = fullfile(fileparts(mfilename('fullpath')), '..', 'src');
    fxp = P.FxpConfig_CD;

    % clean mex file
    mexFile = fullfile(srcDir, '+cd_eq', 'equalize_td_fxp_mex.mexw64');
    if isfile(mexFile)
        delete(mexFile);
        fprintf('  Deleted existing MEX file: %s\n', mexFile);
    end

    T_cd = cd_eq.equalize_td_fxp_types(fxp);

    x_cd = fi(complex(0,0), numerictype(T_cd.x), fimath(T_cd.x));
    In_cd_type = coder.typeof(x_cd, [Inf, P.N_pol], [true, false]);

    args_cd = { ...
        In_cd_type, ...                     % In
        double(P.D), ...                    % D
        double(P.L), ...                    % L
        double(P.CWL), ...                  % CLambda
        double(P.Rs), ...                   % Rs
        double(P.N_pol), ...                % NPol
        double(P.SpS), ...                  % SpSIn
        double(0), ...                      % NTap (runtime-variable)
        T_cd};                              % T

    codegen('-config', cfg, ...
            'cd_eq.equalize_td_fxp', ...
            '-args', args_cd, ...
            '-o', fullfile(srcDir, '+cd_eq', 'equalize_td_fxp_mex'));
    fprintf('  cd_eq.equalize_td_fxp_mex  OK\n');
end
