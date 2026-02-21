function build_cr_viterbiViterbi_fxp_mex(P, cfg)
%BUILD_CR_VITERBIVITERBI_FXP_MEX  Compile cr_viterbiViterbi_fxp to MEX.
%
%   build_cr_viterbiViterbi_fxp_mex(P, cfg)
%
%   Inputs
%     P   - parameter struct from pipeline_params()
%     cfg - coder.MexCodeConfig object

    srcDir = fullfile(fileparts(mfilename('fullpath')), '..', 'src');
    fxp = P.FxpConfig_VV;

    % Clean existing MEX file
    mexFile = fullfile(srcDir, 'cr_viterbiViterbi_fxp_mex.mexw64');
    if isfile(mexFile)
        delete(mexFile);
        fprintf('  Deleted existing MEX file: %s\n', mexFile);
    end

    T_vv = cr_viterbiViterbi_fxp_types(fxp);

    % ----------------------------------------------------------------
    % x  –  variable-length complex fi matrix [N x NPol]
    %        Row count is unbounded; column count is fixed to N_pol.
    % ----------------------------------------------------------------
    x_vv       = fi(complex(0, 0), numerictype(T_vv.x), fimath(T_vv.x));
    In_vv_type = coder.typeof(x_vv, [Inf, P.N_pol], [true, false]);

    % ----------------------------------------------------------------
    % VVFilter  –  variable-length real fi column vector [(2*NTaps+1) x 1]
    %              Length is fixed at build time but declared variable so
    %              the same MEX works if NTaps is tuned without recompile.
    % ----------------------------------------------------------------
    w_proto = fi(0, numerictype(T_vv.w), fimath(T_vv.w));
    w_type  = coder.typeof(w_proto, [Inf, 1], [true, false]);

    % ----------------------------------------------------------------
    % Pilots  –  fixed-length complex fi column vector [PilotLen x 1]
    %            Size is known at compile time (same fi type as x).
    % ----------------------------------------------------------------
    pilots_proto = fi(complex(0, 0), numerictype(T_vv.x), fimath(T_vv.x));
    pilots_type  = coder.typeof(pilots_proto, [P.PilotLen, 1], [false, false]);

    % ----------------------------------------------------------------
    % Build argument list
    % ----------------------------------------------------------------
    args_vv = { ...
        In_vv_type, ...            % x          [N x NPol]    fi complex
        double(P.N_pol), ...       % NPol        scalar        double
        double(P.VV_NTaps), ...    % NTaps       scalar        double
        w_type, ...                % VVFilter   [L_filt x 1]  fi real
        pilots_type, ...           % Pilots     [PilotLen x 1] fi complex
        double(P.BlockLen), ...    % L           scalar        double
        logical(false), ...        % UsePilots   scalar        logical
        logical(false), ...        % BlockBased  scalar        logical
        T_vv};                     % T           struct of fi prototypes

    codegen('-config', cfg, ...
            'cr_viterbiViterbi_fxp', ...
            '-args', args_vv, ...
            '-o', fullfile(srcDir, 'cr_viterbiViterbi_fxp_mex'));
    fprintf('  cr_viterbiViterbi_fxp_mex  OK\n');
end