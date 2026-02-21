function build_cr_bps_fxp_mex(P, cfg)
%BUILD_CR_BPS_FXP_MEX  Compile cr_bps_fxp to MEX.
%
%   build_cr_bps_fxp_mex(P, cfg)
%
%   Inputs
%     P   - parameter struct from pipeline_params()
%     cfg - coder.MexCodeConfig object
%
%   Required fields in P (add to pipeline_params if not present)
%     P.BPS_N         - one-sided BPS filter half-length (window = 2*N+1)
%     P.BPS_B         - number of blind test phases (must be even)
%     P.FxpConfig_BPS - fixed-point config string: 'fixed16' | 'fixed32'
%
%   Codegen note on M
%     M is passed as a plain double scalar.  The BPS decision step uses
%     qam_slicer() which implements nearest-neighbour QAM decisions via
%     pure arithmetic, removing any dependency on qamdemod / qammod and
%     their compile-time M requirement.

    srcDir = fullfile(fileparts(mfilename('fullpath')), '..', 'src');
    fxp = P.FxpConfig_BPS;

    % Clean existing MEX file
    mexFile = fullfile(srcDir, 'cr_bps_fxp_mex.mexw64');
    if isfile(mexFile)
        delete(mexFile);
        fprintf('  Deleted existing MEX file: %s\n', mexFile);
    end

    T_bps = cr_bps_fxp_types(fxp);

    % ----------------------------------------------------------------
    % z  –  variable-length complex fi matrix [N x NPol]
    %        Row count is unbounded; column count fixed to N_pol.
    % ----------------------------------------------------------------
    z_proto    = fi(complex(0, 0), numerictype(T_bps.x), fimath(T_bps.x));
    In_bps_type = coder.typeof(z_proto, [Inf, P.N_pol], [true, false]);

    % ----------------------------------------------------------------
    % Pilots  –  fixed-length complex fi column vector [PilotLen x 1]
    %            Size is a design constant; declared fixed so codegen
    %            can unroll the pilot correlation loop.
    % ----------------------------------------------------------------
    pilots_proto = fi(complex(0, 0), numerictype(T_bps.x), fimath(T_bps.x));
    pilots_type  = coder.typeof(pilots_proto, [P.PilotLen, 1], [false, false]);

    % ----------------------------------------------------------------
    % Build argument list
    % ----------------------------------------------------------------
    args_bps = { ...
        In_bps_type, ...               % z          [Nsym x NPol]  fi complex
        double(P.BPS_N), ...           % N           scalar         double
        double(P.N_pol), ...           % NPol        scalar         double
        double(P.M), ...               % M           scalar         double
        double(P.BPS_B), ...           % B           scalar         double
        double(P.BlockLen), ...        % BlockLen    scalar         double
        pilots_type, ...               % Pilots     [PilotLen x 1]  fi complex
        logical(false), ...            % UsePilots   scalar         logical
        logical(false), ...            % BlockBased  scalar         logical
        T_bps};                        % T           struct of fi prototypes

    codegen('-config', cfg, ...
            'cr_bps_fxp', ...
            '-args', args_bps, ...
            '-o', fullfile(srcDir, 'cr_bps_fxp_mex'));
    fprintf('  cr_bps_fxp_mex  OK\n');
end