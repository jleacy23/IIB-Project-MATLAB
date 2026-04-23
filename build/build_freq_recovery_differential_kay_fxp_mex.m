function build_freq_recovery_differential_kay_fxp_mex(P, cfg)
%BUILD_FREQ_RECOVERY_DIFFERENTIAL_KAY_FXP_MEX  Compile differential_kay_fxp to MEX.
%
%   build_freq_recovery_differential_kay_fxp_mex(P, cfg)
%
%   Compiles the training-aided entry point (data_aided = true).
%   Internally pulls in differential_fxp and tretter_kay_fxp.
%
%   Inputs
%     P   - parameter struct with the fields listed below
%     cfg - coder.MexCodeConfig object
%
%   Required fields in P
%     P.N_pol        - number of polarisations
%     P.TrainingLen  - number of training symbols at subframe start
%     P.Rs           - symbol rate [GBd]
%     P.FxpConfig_FR - fixed-point config: 'fixed16' | 'fixed32' | struct('WL',wl,'FL',fl)
%     P.CordicIts    - number of CORDIC iterations
%     P.MaxFreq      - phase-scaling factor (double)

    srcDir = fullfile(fileparts(mfilename('fullpath')), '..', 'src');
    fxp = P.FxpConfig_FR;

    % Clean existing MEX file
    mexFile = fullfile(srcDir, '+freq_recovery', 'differential_kay_fxp_mex.mexw64');
    if isfile(mexFile)
        delete(mexFile);
        fprintf('  Deleted existing MEX file: %s\n', mexFile);
    end

    T_fr = freq_recovery.fxp_types(fxp);

    % ----------------------------------------------------------------
    % x  –  variable-length complex fi matrix [Nsym x NPol]
    % ----------------------------------------------------------------
    x_proto    = fi(complex(0, 0), numerictype(T_fr.x), fimath(T_fr.x));
    In_fr_type = coder.typeof(x_proto, [Inf, P.N_pol], [true, false]);

    % ----------------------------------------------------------------
    % training  –  fixed-length complex fi matrix [TrainingLen x NPol]
    % ----------------------------------------------------------------
    tr_proto = fi(complex(0, 0), numerictype(T_fr.x), fimath(T_fr.x));
    tr_type  = coder.typeof(tr_proto, [P.TrainingLen, P.N_pol], [false, false]);

    cordic_its_type = coder.Constant(P.CordicIts);

    % ----------------------------------------------------------------
    % Build argument list — matches differential_kay_fxp signature:
    %   (x, training, Rs, CordicIts, T, data_aided, D, max_freq)
    % Explicitly compile the training-aided path with max_freq.
    % ----------------------------------------------------------------
    args = { ...
        In_fr_type, ...          % x           [Nsym x NPol]       fi complex
        tr_type, ...             % training    [TrainingLen x NPol] fi complex
        double(P.Rs), ...        % Rs           scalar              double  [GBd]
        cordic_its_type, ...     % CordicIts    scalar              constant
        T_fr, ...                % T            struct of fi prototypes
        true, ...                % data_aided   scalar              logical
        0, ...                   % D            scalar              double (unused)
        double(P.MaxFreq)};      % max_freq     scalar              double

    codegen('-config', cfg, ...
            'freq_recovery.differential_kay_fxp', ...
            '-args', args, ...
            '-o', fullfile(srcDir, '+freq_recovery', 'differential_kay_fxp_mex'));
    fprintf('  freq_recovery.differential_kay_fxp_mex  OK\n');
end
