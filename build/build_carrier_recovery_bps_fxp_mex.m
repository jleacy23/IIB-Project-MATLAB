function build_carrier_recovery_bps_fxp_mex(P, cfg)
%BUILD_carrier_recovery.bps_FXP_MEX  Compile carrier_recovery.bps_fxp to MEX.
%
%   build_carrier_recovery_bps_fxp_mex(P, cfg)
%
%   Inputs
%     P   - parameter struct with the fields listed below
%     cfg - coder.MexCodeConfig object
%
%   Required fields in P
%     P.N_pol         - number of polarisations
%     P.BPS_N         - one-sided BPS filter half-length (window = 2*N+1)
%     P.BPS_B         - number of blind test phases (must be even)
%     P.PilotLen      - number of pilot symbols per block
%     P.BlockLen      - block length in symbols
%     P.StepSize      - phase update interval in symbols (1..BlockLen)
%     P.M             - QAM order
%     P.FxpConfig_BPS - fixed-point config string: 'fixed16' | 'fixed32'
%     P.PilotThreshold - threshold for pilot-based cycle-slip correction in radians
%     P.COrdicIts      - number of iterations for CORDIC operations
%
%   Codegen note on M
%     M is passed as a plain double scalar.  The BPS decision step uses
%     modem.slicer() which implements nearest-neighbour QAM decisions via
%     pure arithmetic, removing any dependency on qamdemod / qammod and
%     their compile-time M requirement.

    srcDir = fullfile(fileparts(mfilename('fullpath')), '..', 'src');
    fxp = P.FxpConfig_BPS;

    % Clean existing MEX file
    mexFile = fullfile(srcDir, '+carrier_recovery', 'bps_fxp_mex.mexw64');
    if isfile(mexFile)
        delete(mexFile);
        fprintf('  Deleted existing MEX file: %s\n', mexFile);
    end

    T_bps = carrier_recovery.bps_fxp_types(fxp);

    % ----------------------------------------------------------------
    % z  –  variable-length complex fi matrix [Nsym x NPol]
    %        Row count is unbounded; column count fixed to N_pol.
    % ----------------------------------------------------------------
    z_proto     = fi(complex(0, 0), numerictype(T_bps.x), fimath(T_bps.x));
    In_bps_type = coder.typeof(z_proto, [Inf, P.N_pol], [true, false]);

    % ----------------------------------------------------------------
    % Pilots  –  fixed-length complex fi column vector [PilotLen x 1]
    %            Size is known at compile time; declared fixed so codegen
    %            can unroll the pilot correlation loop.
    % ----------------------------------------------------------------
    pilots_proto = fi(complex(0, 0), numerictype(T_bps.x), fimath(T_bps.x));
    pilots_type  = coder.typeof(pilots_proto, [P.PilotLen, 1], [false, false]);

    cordic_its_type = coder.Constant(P.CordicIts);

    % ----------------------------------------------------------------
    % Build argument list — must match carrier_recovery.bps_fxp signature:
    %   (z, N, NPol, M, B, BlockLen, StepSize, Pilots, UsePilots, PilotThreshold, T)
    % ----------------------------------------------------------------
    args_bps = { ...
        In_bps_type, ...               % z          [Nsym x NPol]   fi complex
        double(P.BPS_N), ...           % N           scalar          double
        double(P.N_pol), ...           % NPol        scalar          double
        double(P.M), ...               % M           scalar          double
        double(P.BPS_B), ...           % B           scalar          double
        double(P.BlockLen), ...        % BlockLen    scalar          double
        double(P.StepSize), ...        % StepSize    scalar          double
        pilots_type, ...               % Pilots     [PilotLen x 1]  fi complex
        logical(false), ...            % UsePilots   scalar          logical
        double(P.PilotThreshold), ...  % PilotThreshold scalar       double
        cordic_its_type, ...       % CordicIts   scalar          double
        T_bps};                        % T           struct of fi prototypes

    codegen('-config', cfg, ...
            'carrier_recovery.bps_fxp', ...
            '-args', args_bps, ...
            '-o', fullfile(srcDir, '+carrier_recovery', 'bps_fxp_mex'));
    fprintf('  carrier_recovery.bps_fxp_mex  OK\n');
end