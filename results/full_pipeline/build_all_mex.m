%% build_all_mex.m
%  Master build script: compiles MEX files for every fixed-point DSP block
%  using the parameters defined in pipeline_params.
%
%  Prerequisites: MATLAB Coder, Fixed-Point Designer
%
%  Usage:  run this script from the project root, build/, or
%          results/full_pipeline/.

%% 0 — Paths
thisDir = fileparts(mfilename('fullpath'));
srcDir  = fullfile(thisDir, '..', '..', 'src');
addpath(srcDir);
addpath(thisDir);                      % so pipeline_params is available

P = pipeline_params();
fxp = P.FxpConfig;                     % e.g. 'fixed32'

fprintf('\n========================================\n');
fprintf('  build_all_mex  |  config = ''%s''\n', fxp);
fprintf('========================================\n');

cfg = coder.config('mex');
cfg.GenerateReport      = true;
cfg.EnableMexProfiling  = false;

%% =================================================================
%  1  — CD Equalizer  (cdeq_equalize_fxp)
% ==================================================================
fprintf('\n--- [1/3] cdeq_equalize_fxp ---\n');

T_cd = cdeq_equalize_fxp_types(fxp);

x_cd = fi(complex(0,0), numerictype(T_cd.x), fimath(T_cd.x));
In_cd_type = coder.typeof(x_cd, [Inf, P.N_pol], [true, false]);

args_cd = { ...
    In_cd_type, ...                     % In
    double(P.D), ...                    % D
    double(P.L), ...                    % L
    double(P.CWL), ...                  % CLambda
    double(P.Rs), ...                   % Rs
    double(P.N_pol), ...                % NPol
    double(P.SpS), ...                  % SpSIn
    double(P.NFFT), ...                 % NFFT
    logical(P.po2Twiddle), ...          % po2Twiddle
    T_cd};                              % T

codegen('-config', cfg, ...
        'cdeq_equalize_fxp', ...
        '-args', args_cd, ...
        '-o', fullfile(srcDir, 'cdeq_equalize_fxp_mex'));
fprintf('  cdeq_equalize_fxp_mex  OK\n');

%% =================================================================
%  2  — Adaptive Equalizer  (adeq_equalize_fxp)
% ==================================================================
fprintf('\n--- [2/3] adeq_equalize_fxp ---\n');

T_aeq = adeq_equalize_fxp_types(fxp);

x_aeq = fi(complex(0,0), numerictype(T_aeq.x), fimath(T_aeq.x));
In_aeq_type = coder.typeof(x_aeq, [Inf, 2], [true, false]);

Eq_type = coder.typeof('a', [1, 7], [false, true]);

args_aeq = { ...
    In_aeq_type, ...                    % x
    double(P.SpS), ...                  % SpS
    Eq_type, ...                        % Eq
    double(P.AEQ_NTaps), ...            % NTaps
    double(P.AEQ_Mu), ...               % Mu
    logical(P.AEQ_SingleSpike), ...     % SingleSpike
    double(P.AEQ_N1), ...               % N1
    double(P.AEQ_N2), ...               % N2
    double(P.AEQ_NOut), ...             % NOut
    T_aeq};                             % T

codegen('-config', cfg, ...
        'adeq_equalize_fxp', ...
        '-args', args_aeq, ...
        '-o', fullfile(srcDir, 'adeq_equalize_fxp_mex'));
fprintf('  adeq_equalize_fxp_mex  OK\n');

%% =================================================================
%  3  — Viterbi-Viterbi Carrier Recovery  (cr_viterbiViterbi_fxp)
% ==================================================================
fprintf('\n--- [3/3] cr_viterbiViterbi_fxp ---\n');

T_vv = cr_viterbiViterbi_fxp_types(fxp);

x_vv = fi(complex(0,0), numerictype(T_vv.x), fimath(T_vv.x));
In_vv_type = coder.typeof(x_vv, [Inf, 2], [true, false]);

w_proto = fi(0, numerictype(T_vv.w), fimath(T_vv.w));
w_type  = coder.typeof(w_proto, [Inf, 1], [true, false]);

Pilots_type = coder.typeof(x_vv, [Inf, 2], [true, false]);

args_vv = { ...
    In_vv_type, ...                     % x
    double(P.N_pol), ...                % NPol
    double(P.VV_NTaps), ...             % NTaps
    w_type, ...                         % VVFilter
    Pilots_type, ...                    % Pilots
    double(P.VV_P), ...                 % P
    double(P.VV_L), ...                 % L
    double(P.VV_CSThreshold), ...       % CSThreshold
    logical(P.VV_UsePilots), ...        % UsePilots
    T_vv};                              % T

codegen('-config', cfg, ...
        'cr_viterbiViterbi_fxp', ...
        '-args', args_vv, ...
        '-o', fullfile(srcDir, 'cr_viterbiViterbi_fxp_mex'));
fprintf('  cr_viterbiViterbi_fxp_mex  OK\n');

%% =================================================================
fprintf('\n========================================\n');
fprintf('  All MEX files built successfully.\n');
fprintf('========================================\n');
