%% build_all_mex.m
%  Master build script: compiles MEX files for every fixed-point DSP block
%  using the parameters defined in pipeline_params.
%
%  Prerequisites: MATLAB Coder, Fixed-Point Designer
%
%  Usage:  run this script from the project root, build/, or
%          results/full_pipeline/.

%% 0 — Paths
thisDir  = fileparts(mfilename('fullpath'));
srcDir   = fullfile(thisDir, '..', '..', 'src');
buildDir = fullfile(thisDir, '..', '..', 'build');
addpath(srcDir);
addpath(thisDir);                      % so pipeline_params is available
addpath(buildDir);                     % individual build scripts

P = pipeline_params();

fprintf('\n========================================\n');
fprintf('  build_all_mex\n');
fprintf('    CD  = ''%s''  |  AEQ = ''%s''  |  VV = ''%s''\n', ...
        P.FxpConfig_CD, P.FxpConfig_AEQ, P.FxpConfig_VV);
fprintf('========================================\n');

cfg = coder.config('mex');
cfg.GenerateReport      = true;
cfg.EnableMexProfiling  = false;

%% =================================================================
%  1  — CD Equalizer  (cd_eq.equalize_fxp)
% ==================================================================
fprintf('\n--- [1/3] cd_eq.equalize_fxp ---\n');
build_cd_eq_equalize_fxp_mex(P, cfg);

%% =================================================================
%  2  — Adaptive Equalizer  (adaptive_eq.equalize_fxp)
% ==================================================================
fprintf('\n--- [2/3] adaptive_eq.equalize_fxp ---\n');
build_adaptive_eq_equalize_fxp_mex(P, cfg);

%% =================================================================
%  3  — Viterbi-Viterbi Carrier Recovery  (carrier_recovery.viterbiViterbi_fxp)
% ==================================================================
fprintf('\n--- [3/3] carrier_recovery.viterbiViterbi_fxp ---\n');
build_carrier_recovery_viterbiViterbi_fxp_mex(P, cfg);

%% =================================================================
fprintf('\n========================================\n');
fprintf('  All MEX files built successfully.\n');
fprintf('========================================\n');
