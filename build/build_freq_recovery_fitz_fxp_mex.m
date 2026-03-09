function build_freq_recovery_fitz_fxp_mex(P, cfg)
%BUILD_FREQ_RECOVERY_FITZ_FXP_MEX  Compile fitz_fxp to MEX.
%
%   build_freq_recovery_fitz_fxp_mex(P, cfg)
%
%   Inputs
%     P   - parameter struct with the fields listed below
%     cfg - coder.MexCodeConfig object
%
%   Required fields in P
%     P.N_pol          - number of polarisations
%     P.TrainingLen    - number of training symbols at subframe start
%     P.FR_Fitz_N      - autocorrelation lag (double, 1 <= N < TrainingLen)
%     P.FxpConfig_FR   - fixed-point config string: 'fixed16' | 'fixed32'
%     P.CordicIts      - number of CORDIC iterations

    srcDir = fullfile(fileparts(mfilename('fullpath')), '..', 'src');
    fxp = P.FxpConfig_FR;

    % Clean existing MEX file
    mexFile = fullfile(srcDir, '+freq_recovery', 'fitz_fxp_mex.mexw64');
    if isfile(mexFile)
        delete(mexFile);
        fprintf('  Deleted existing MEX file: %s\n', mexFile);
    end

    T_fr = freq_recovery.fxp_types(fxp);

    % ----------------------------------------------------------------
    % x  –  variable-length complex fi matrix [Nsym x NPol]
    %        Row count is unbounded; column count fixed to N_pol.
    % ----------------------------------------------------------------
    x_proto    = fi(complex(0, 0), numerictype(T_fr.x), fimath(T_fr.x));
    In_fr_type = coder.typeof(x_proto, [Inf, P.N_pol], [true, false]);

    % ----------------------------------------------------------------
    % training  –  fixed-length complex fi matrix [TrainingLen x NPol]
    %              Size is known at compile time.
    % ----------------------------------------------------------------
    tr_proto   = fi(complex(0, 0), numerictype(T_fr.x), fimath(T_fr.x));
    tr_type    = coder.typeof(tr_proto, [P.TrainingLen, P.N_pol], [false, false]);

    cordic_its_type = coder.Constant(P.CordicIts);

    % ----------------------------------------------------------------
    % Build argument list — must match freq_recovery.fitz_fxp signature:
    %   (x, training, Rs, N, CordicIts, T)
    % ----------------------------------------------------------------
    args = { ...
        In_fr_type, ...          % x         [Nsym x NPol]       fi complex
        tr_type, ...             % training  [TrainingLen x NPol] fi complex
        double(P.Rs), ...        % Rs         scalar              double  [GBd]
        double(P.FR_Fitz_N), ... % N          scalar              double  (autocorrelation lag)
        cordic_its_type, ...     % CordicIts  scalar              double
        T_fr};                   % T          struct of fi prototypes

    codegen('-config', cfg, ...
            'freq_recovery.fitz_fxp', ...
            '-args', args, ...
            '-o', fullfile(srcDir, '+freq_recovery', 'fitz_fxp_mex'));
    fprintf('  freq_recovery.fitz_fxp_mex  OK\n');
end
