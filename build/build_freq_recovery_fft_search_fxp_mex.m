function build_freq_recovery_fft_search_fxp_mex(P, cfg)
%BUILD_FREQ_RECOVERY_FFT_SEARCH_FXP_MEX  Compile fft_search_fxp to MEX.
%
%   build_freq_recovery_fft_search_fxp_mex(P, cfg)
%
%   Inputs
%     P   - parameter struct with the fields listed below
%     cfg - coder.MexCodeConfig object
%
%   Required fields in P
%     P.N_pol           - number of polarisations
%     P.TrainingLen     - number of training symbols at subframe start
%     P.FR_Nfft         - FFT size (integer power of 2, >= TrainingLen)
%     P.FR_Po2Twiddle   - logical: use power-of-2 twiddle factors in fft_fxp
%     P.FxpConfig_FR    - fixed-point config string: 'fixed16' | 'fixed32'
%     P.CordicIts       - number of CORDIC iterations
%     P.MaxFreq         - phase-scaling factor (double)
%     P.FR_BlindD       - blind data length [symbols] (scalar, for type)
%
%   Note on FFT size
%     fft.fft_fxp requires N to be a power of 2.  P.FR_Nfft must satisfy
%     this and be >= P.TrainingLen.  With TrainingLen = 11, the minimum
%     useful choice is 16; Nfft = 128 gives ~11x zero-padding.

    srcDir = fullfile(fileparts(mfilename('fullpath')), '..', 'src');
    fxp = P.FxpConfig_FR;

    % Clean existing MEX file
    mexFile = fullfile(srcDir, '+freq_recovery', 'fft_search_fxp_mex.mexw64');
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
    % Build argument list — must match freq_recovery.fft_search_fxp signature:
    %   (x, training, Rs, Nfft, po2Twiddle, CordicIts, max_freq, T, data_aided, D)
    % ----------------------------------------------------------------
    args = { ...
        In_fr_type, ...                   % x           [Nsym x NPol]       fi complex
        tr_type, ...                      % training    [TrainingLen x NPol] fi complex
        double(P.Rs), ...                 % Rs           scalar              double  [GBd]
        double(P.FR_Nfft), ...            % Nfft         scalar              double  (power of 2)
        logical(P.FR_Po2Twiddle), ...     % po2Twiddle   scalar              logical
        cordic_its_type, ...              % CordicIts    scalar              double
        double(P.MaxFreq), ...            % max_freq     scalar              double
        T_fr, ...                         % T            struct of fi prototypes
        true, ...                         % data_aided   scalar              logical
        double(P.FR_BlindD)};             % D            scalar              double

    codegen('-config', cfg, ...
            'freq_recovery.fft_search_fxp', ...
            '-args', args, ...
            '-o', fullfile(srcDir, '+freq_recovery', 'fft_search_fxp_mex'));
    fprintf('  freq_recovery.fft_search_fxp_mex  OK\n');
end
