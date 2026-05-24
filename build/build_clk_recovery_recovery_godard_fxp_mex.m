function build_clk_recovery_recovery_godard_fxp_mex(P, cfg)
%BUILD_CLK_RECOVERY_RECOVERY_GODARD_FXP_MEX  Compile recovery_godard_fxp to MEX.
%
%   build_clk_recovery_recovery_godard_fxp_mex(P, cfg)
%
%   Inputs
%     P   - parameter struct from pipeline_params()
%     cfg - coder.MexCodeConfig object
%
%   Required fields in P
%     P.CR_NFFT       - FFT block size (power of 2)
%     P.Rolloff       - pulse-shaping roll-off factor (0 < beta <= 1)
%     P.po2Twiddle    - logical: round FFT twiddles to nearest power of 2
%     P.FxpConfig_CR  - fixed-point config: 'fixed16' | 'fixed32' | struct('WL',wl,'FL',fl)

    srcDir = fullfile(fileparts(mfilename('fullpath')), '..', 'src');
    fxp    = P.FxpConfig_CR;

    % Clean existing MEX file
    mexFile = fullfile(srcDir, '+clk_recovery', 'recovery_godard_fxp_mex.mexw64');
    if isfile(mexFile)
        delete(mexFile);
        fprintf('  Deleted existing MEX file: %s\n', mexFile);
    end

    T_cr = clk_recovery.recovery_godard_fxp_types(fxp);

    % ----------------------------------------------------------------
    % In  –  variable-length complex fi column [Nsamp x 1]
    %         recovery_godard_fxp processes one polarisation at a time.
    % ----------------------------------------------------------------
    x_proto    = fi(complex(0, 0), numerictype(T_cr.x), fimath(T_cr.x));
    In_cr_type = coder.typeof(x_proto, [Inf, 1], [true, false]);

    % ----------------------------------------------------------------
    % Argument list — must match recovery_godard_fxp signature:
    %   (In, NSymb, N, beta, po2Twiddle, T)
    % ----------------------------------------------------------------
    args = { ...
        In_cr_type, ...                    % In           [Nsamp x 1]  fi complex
        double(P.Ns), ...                  % NSymb        scalar       double
        double(P.CR_NFFT), ...             % N            scalar       double  (power of 2)
        double(P.Rolloff), ...             % beta         scalar       double
        logical(P.po2Twiddle), ...         % po2Twiddle   scalar       logical
        T_cr};                             % T            struct of fi prototypes

    codegen('-config', cfg, ...
            'clk_recovery.recovery_godard_fxp', ...
            '-args', args, ...
            '-o', fullfile(srcDir, '+clk_recovery', 'recovery_godard_fxp_mex'));
    fprintf('  clk_recovery.recovery_godard_fxp_mex  OK\n');
end
