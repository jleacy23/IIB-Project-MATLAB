function build_adaptive_eq_equalize_fxp_mex(P, cfg)
%BUILD_ADEQ_EQUALIZE_FXP_MEX  Compile adeq_equalize_fxp to MEX.
%
%   build_adaptive_eq_equalize_fxp_mex(P, cfg)
%
%   Inputs
%     P   - parameter struct from pipeline_params()
%     cfg - coder.MexCodeConfig object
%
%   Required fields in P
%     P.FxpConfig_AEQ - fixed-point config: 'fixed16' | 'fixed32' | struct('WL',wl,'FL',fl)

    srcDir = fullfile(fileparts(mfilename('fullpath')), '..', 'src');
    fxp = P.FxpConfig_AEQ;

    % clean mex file
    mexFile = fullfile(srcDir, '+adaptive_eq', 'equalize_fxp_mex.mexw64');
    if isfile(mexFile)
        delete(mexFile);
        fprintf('  Deleted existing MEX file: %s\n', mexFile);
    end

    T_aeq = adaptive_eq.equalize_fxp_types(fxp);

    x_aeq = fi(complex(0,0), numerictype(T_aeq.x), fimath(T_aeq.x));
    In_aeq_type = coder.typeof(x_aeq, [Inf, 2], [true, false]);

    Eq_type = coder.typeof('a', [1, 7], [false, true]);

    args_aeq = { ...
        In_aeq_type, ...                    % x
        double(P.SpS), ...                  % SpS
        Eq_type, ...                        % Eq
        double(P.AEQ_NTaps), ...            % NTaps
        double(P.AEQ_Mu), ...               % Mu
        logical(P.AEQ_SingleSpike), ...     % SingleSpike
        double(P.AEQ_N1), ...               % N1
        double(P.AEQ_N2), ...               % N2
        double(P.AEQ_NOut), ...             % NOut
        T_aeq};                             % T

    codegen('-config', cfg, ...
            'adaptive_eq.equalize_fxp', ...
            '-args', args_aeq, ...
            '-o', fullfile(srcDir, '+adaptive_eq', 'equalize_fxp_mex'));
    fprintf('  adaptive_eq.equalize_fxp_mex  OK\n');
end
