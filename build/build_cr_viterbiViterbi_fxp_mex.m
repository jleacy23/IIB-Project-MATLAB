function build_cr_viterbiViterbi_fxp_mex(P, cfg)
%BUILD_CR_VITERBIVITERBI_FXP_MEX  Compile cr_viterbiViterbi_fxp to MEX.
%
%   build_cr_viterbiViterbi_fxp_mex(P, cfg)
%
%   Inputs
%     P   - parameter struct from pipeline_params()
%     cfg - coder.MexCodeConfig object

    srcDir = fullfile(fileparts(mfilename('fullpath')), '..', 'src');
    fxp = P.FxpConfig;

    % clean mex file
    mexFile = fullfile(srcDir, 'cr_viterbiViterbi_fxp_mex.mexw64');
    if isfile(mexFile)
        delete(mexFile);
        fprintf('  Deleted existing MEX file: %s\n', mexFile);
    end

    T_vv = cr_viterbiViterbi_fxp_types(fxp);

    x_vv = fi(complex(0,0), numerictype(T_vv.x), fimath(T_vv.x));
    In_vv_type = coder.typeof(x_vv, [Inf, 2], [true, false]);

    w_proto = fi(0, numerictype(T_vv.w), fimath(T_vv.w));
    w_type  = coder.typeof(w_proto, [Inf, 1], [true, false]);

    args_vv = { ...
        In_vv_type, ...                     % x
        double(P.N_pol), ...                % NPol
        double(P.VV_NTaps), ...             % NTaps
        w_type, ...                         % VVFilter
        T_vv};                              % T

    codegen('-config', cfg, ...
            'cr_viterbiViterbi_fxp', ...
            '-args', args_vv, ...
            '-o', fullfile(srcDir, 'cr_viterbiViterbi_fxp_mex'));
    fprintf('  cr_viterbiViterbi_fxp_mex  OK\n');
end
