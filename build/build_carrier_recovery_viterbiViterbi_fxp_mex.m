function build_carrier_recovery_viterbiViterbi_fxp_mex(P, cfg)
%BUILD_CR_VITERBIVITERBI_FXP_MEX  Compile cr_viterbiViterbi_fxp to MEX.
%
%   build_carrier_recovery_viterbiViterbi_fxp_mex(P, cfg)
%
%   Inputs
%     P   - parameter struct with the fields listed below
%     cfg - coder.MexCodeConfig object
%
%   Required fields in P
%     P.N_pol        - number of polarisations
%     P.VV_NTaps     - one-sided VV filter half-length
%     P.PilotLen     - number of pilot symbols per block
%     P.BlockLen     - block length in symbols
%     P.StepSize     - phase update interval in symbols (1..BlockLen)
%     P.FxpConfig_VV - fixed-point config: 'fixed16' | 'fixed32' | struct('WL',wl,'FL',fl)
%     P.PilotThreshold - threshold for pilot-based cycle-slip correction in radians
%     P.CordicIts      - number of iterations for CORDIC operations

    srcDir = fullfile(fileparts(mfilename('fullpath')), '..', 'src');
    fxp = P.FxpConfig_VV;

    % Clean existing MEX file
    mexFile = fullfile(srcDir, '+carrier_recovery', 'viterbiViterbi_fxp_mex.mexw64');
    if isfile(mexFile)
        delete(mexFile);
        fprintf('  Deleted existing MEX file: %s\n', mexFile);
    end

    T_vv = carrier_recovery.fxp_types(fxp);

    % ----------------------------------------------------------------
    % x  –  variable-length complex fi matrix [N x NPol]
    %        Row count is unbounded; column count is fixed to N_pol.
    % ----------------------------------------------------------------
    x_vv       = fi(complex(0, 0), numerictype(T_vv.x), fimath(T_vv.x));
    In_vv_type = coder.typeof(x_vv, [Inf, P.N_pol], [true, false]);

    % ----------------------------------------------------------------
    % VVFilter  –  variable-length real fi column vector [(2*NTaps+1) x 1]
    %              Declared variable so the same MEX works if NTaps is
    %              tuned without recompilation.
    % ----------------------------------------------------------------
    w_proto = fi(0, numerictype(T_vv.w), fimath(T_vv.w));
    w_type  = coder.typeof(w_proto, [Inf, 1], [true, false]);

    % ----------------------------------------------------------------
    % Pilots  –  variable-size complex fi matrix [NBlocks x NPol]
    %            Row count varies with signal length; column count fixed.
    % ----------------------------------------------------------------
    pilots_proto = fi(complex(0, 0), numerictype(T_vv.x), fimath(T_vv.x));
    pilots_type  = coder.typeof(pilots_proto, [Inf, P.N_pol], [true, false]);

    cordic_its_type = coder.Constant(P.CordicIts);

    % ----------------------------------------------------------------
    % Build argument list — must match cr_viterbiViterbi_fxp signature:
    %   (x, NPol, NTaps, VVFilter, Pilots, BlockLen, StepSize, PilotThreshold, CordicIts, T)
    % ----------------------------------------------------------------
    args_vv = { ...
        In_vv_type, ...            % x          [N x NPol]       fi complex
        double(P.N_pol), ...       % NPol        scalar           double
        double(P.VV_NTaps), ...    % NTaps       scalar           double
        w_type, ...                % VVFilter   [L_filt x 1]     fi real
        pilots_type, ...           % Pilots     [NBlocks x NPol]  fi complex
        double(P.BlockLen), ...    % BlockLen    scalar           double
        double(P.StepSize), ...    % StepSize    scalar           double
        double(P.PilotThreshold), ... % PilotThreshold scalar       double
        cordic_its_type, ...       % CordicIts   scalar           double
        T_vv};                     % T           struct of fi prototypes

    codegen('-config', cfg, ...
            'carrier_recovery.viterbiViterbi_fxp', ...
            '-args', args_vv, ...
            '-o', fullfile(srcDir, '+carrier_recovery', 'viterbiViterbi_fxp_mex'));
    fprintf('  carrier_recovery.viterbiViterbi_fxp_mex  OK\n');
end