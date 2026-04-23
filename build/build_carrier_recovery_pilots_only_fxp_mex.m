function build_carrier_recovery_pilots_only_fxp_mex(P, cfg)
%BUILD_CARRIER_RECOVERY_PILOTS_ONLY_FXP_MEX  Compile pilots_only_fxp to MEX.
%
%   build_carrier_recovery_pilots_only_fxp_mex(P, cfg)
%
%   Inputs
%     P   - parameter struct with the fields listed below
%     cfg - coder.MexCodeConfig object
%
%   Required fields in P
%     P.N_pol        - number of polarisations
%     P.BlockLen     - block length in symbols
%     P.FxpConfig_PO - fixed-point config: 'fixed16' | 'fixed32' | struct('WL',wl,'FL',fl)
%     P.CordicIts    - number of iterations for CORDIC operations

    srcDir = fullfile(fileparts(mfilename('fullpath')), '..', 'src');
    fxp = P.FxpConfig_PO;

    % Clean existing MEX file
    mexFile = fullfile(srcDir, '+carrier_recovery', 'pilots_only_fxp_mex.mexw64');
    if isfile(mexFile)
        delete(mexFile);
        fprintf('  Deleted existing MEX file: %s\n', mexFile);
    end

    T_po = carrier_recovery.fxp_types(fxp);

    % ----------------------------------------------------------------
    % x  - variable-length complex fi matrix [Nsym x NPol]
    %      Row count is unbounded; column count fixed to N_pol.
    % ----------------------------------------------------------------
    x_proto  = fi(complex(0, 0), numerictype(T_po.x), fimath(T_po.x));
    x_type   = coder.typeof(x_proto, [Inf, P.N_pol], [true, false]);

    % ----------------------------------------------------------------
    % Pilots - variable-size complex fi matrix [NBlocks x NPol]
    %          Row count varies with signal length; column count fixed.
    % ----------------------------------------------------------------
    p_proto  = fi(complex(0, 0), numerictype(T_po.x), fimath(T_po.x));
    p_type   = coder.typeof(p_proto, [Inf, P.N_pol], [true, false]);

    cordic_its_type = coder.Constant(P.CordicIts);

    % ----------------------------------------------------------------
    % Build argument list — must match pilots_only_fxp signature:
    %   (x, NPol, BlockLen, Pilots, CordicIts, T)
    % ----------------------------------------------------------------
    args_po = { ...
        x_type, ...                 % x          [Nsym x NPol]       fi complex
        double(P.N_pol), ...        % NPol        scalar              double
        double(P.BlockLen), ...     % BlockLen    scalar              double
        p_type, ...                 % Pilots     [NBlocks x NPol]     fi complex
        cordic_its_type, ...        % CordicIts   scalar              constant
        T_po};                      % T           struct of fi prototypes

    codegen('-config', cfg, ...
            'carrier_recovery.pilots_only_fxp', ...
            '-args', args_po, ...
            '-o', fullfile(srcDir, '+carrier_recovery', 'pilots_only_fxp_mex'));
    fprintf('  carrier_recovery.pilots_only_fxp_mex  OK\n');
end
