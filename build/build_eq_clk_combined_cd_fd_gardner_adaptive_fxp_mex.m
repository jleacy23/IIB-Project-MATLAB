function build_eq_clk_combined_cd_fd_gardner_adaptive_fxp_mex(P, cfg)
%BUILD_EQ_CLK_COMBINED_CD_FD_GARDNER_ADAPTIVE_FXP_MEX
%   Compile eq_clk.combined_cd_fd_gardner_adaptive_fxp to MEX.
%
%   build_eq_clk_combined_cd_fd_gardner_adaptive_fxp_mex(P, cfg)
%
%   Inputs
%     P   - parameter struct from pipeline_params()
%     cfg - coder.MexCodeConfig object
%
%   Required fields in P
%     P.NFFT, P.NOverlap, P.D, P.L, P.CWL, P.Rs, P.Rolloff, P.SpS, P.N_pol
%     P.CR_ki, P.CR_kp, P.CR_NLanes        (Gardner DPLL)
%     P.AEQ_NTaps, P.AEQ_Mu, P.AEQ_SingleSpike, P.AEQ_N1, P.AEQ_NOut,
%     P.AEQ_SignOnly, P.AEQ_UpdateStep, P.AEQ_PLanes
%     P.po2Twiddle
%     P.cfoEnable                          (logical)
%     P.Ns                                 (number of symbols)
%     P.FxpConfig_CombGardner              - composite fxp config: string,
%                                            struct(WL,FL), or struct with
%                                            sub-fields Static/CFO/Clk/AdaptEq
%
%   Optional fields in P (defaults match the float reference)
%     P.AEQ_Mode, P.AEQ_BlockLen, P.AEQ_SubframeBlocks

    srcDir = fullfile(fileparts(mfilename('fullpath')), '..', 'src');
    fxp    = P.FxpConfig_CombGardner;

    % Clean existing MEX file
    mexFile = fullfile(srcDir, '+eq_clk', ...
        'combined_cd_fd_gardner_adaptive_fxp_mex.mexw64');
    if isfile(mexFile)
        delete(mexFile);
        fprintf('  Deleted existing MEX file: %s\n', mexFile);
    end

    T = eq_clk.combined_cd_fd_gardner_adaptive_fxp_types(fxp);

    %% Variable-length complex fi input at T.Static.x precision
    x_proto    = fi(complex(0,0), numerictype(T.Static.x), fimath(T.Static.x));
    In_type    = coder.typeof(x_proto, [Inf, P.N_pol], [true, false]);

    %% AdaptOpts struct prototype.  Each field given its concrete type so
    %  codegen can specialise.  Pilots is variable-row fi at T.AdaptEq.y.
    p_proto = fi(complex(0,0), numerictype(T.AdaptEq.y), fimath(T.AdaptEq.y));
    Pilots_type = coder.typeof(p_proto, [Inf, 2], [true, false]);

    if isfield(P, 'AEQ_Mode'),           ModeAEQ   = P.AEQ_Mode;           else, ModeAEQ   = 0;             end
    if isfield(P, 'AEQ_BlockLen'),       BLenAEQ   = P.AEQ_BlockLen;       else, BLenAEQ   = P.AEQ_PLanes;  end
    if isfield(P, 'AEQ_SubframeBlocks'), SubBlkAEQ = P.AEQ_SubframeBlocks; else, SubBlkAEQ = 0;             end

    AdaptOptsProto = struct( ...
        'NTaps',          double(P.AEQ_NTaps), ...
        'Mu',             double(P.AEQ_Mu), ...
        'SingleSpike',    logical(P.AEQ_SingleSpike), ...
        'N1',             double(P.AEQ_N1), ...
        'NOut',           double(P.AEQ_NOut), ...
        'SignOnly',       logical(P.AEQ_SignOnly), ...
        'UpdateStep',     double(P.AEQ_UpdateStep), ...
        'PLanes',         double(P.AEQ_PLanes), ...
        'Mode',           double(ModeAEQ), ...
        'Pilots',         p_proto, ...                % concrete fi for type extraction
        'BlockLen',       double(BLenAEQ), ...
        'SubframeBlocks', double(SubBlkAEQ));

    AdaptOptsType = coder.typeof(AdaptOptsProto);
    AdaptOptsType.Fields.Pilots = Pilots_type;        % override Pilots field

    args = { ...
        In_type, ...                       % In
        double(P.SpS), ...                 % SpS
        double(P.NFFT), ...                % NFFT
        double(P.NOverlap), ...            % NOverlap
        double(P.D), ...                   % D
        double(P.L), ...                   % L
        double(P.CWL), ...                 % CLambda
        double(P.Rs), ...                  % Rs
        double(P.Rolloff), ...             % Rolloff
        double(P.CR_ki), ...               % ki
        double(P.CR_kp), ...               % kp
        double(P.Ns), ...                  % NSymb
        double(P.CR_NLanes), ...           % NLanes
        AdaptOptsType, ...                 % AdaptOpts
        logical(P.cfoEnable), ...          % cfoEnable
        logical(P.po2Twiddle), ...         % po2Twiddle
        T};                                % T (composite fxp types)

    codegen('-config', cfg, ...
            'eq_clk.combined_cd_fd_gardner_adaptive_fxp', ...
            '-args', args, ...
            '-o', fullfile(srcDir, '+eq_clk', ...
                'combined_cd_fd_gardner_adaptive_fxp_mex'));
    fprintf('  eq_clk.combined_cd_fd_gardner_adaptive_fxp_mex  OK\n');
end
