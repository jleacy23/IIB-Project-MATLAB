%% build_cr_viterbiViterbi_fxp_mex.m
%  Build script: checks codegen readiness and compiles a MEX function from
%  cr_viterbiViterbi_fxp using MATLAB Coder + Fixed-Point Designer.
%
%  Prerequisites: MATLAB Coder, Fixed-Point Designer
%
%  Usage:  run this script from the project root or from src/.

%% 0 — Ensure src/ is on the MATLAB path
srcDir = fullfile(fileparts(mfilename('fullpath')));
addpath(srcDir);

%% 1 — Static code-generation readiness check
fprintf('\n=== Step 1: Code generation readiness (coder.screener) ===\n');
coder.screener('cr_viterbiViterbi_fxp');
fprintf('  coder.screener completed — review any flagged issues above.\n');

%% 2 — Define example / prototype input types
fprintf('\n=== Step 2: Defining input argument types ===\n');

% Fixed-point types table
T = cr_viterbiViterbi_fxp_types('fixed16');

% x — complex fi input signal, variable-length rows, up to 2 columns
x_proto = fi(complex(0, 0), numerictype(T.x), fimath(T.x));
x_type  = coder.typeof(x_proto, [Inf, 2], [true, false]);

% Scalar parameters (double)
NPol_ex  = double(2);
NTaps_ex = double(15);

% VVFilter — real fi vector, variable-length (L_filt = 2*NTaps+1)
w_proto  = fi(0, numerictype(T.w), fimath(T.w));
w_type   = coder.typeof(w_proto, [Inf, 1], [true, false]);

% T — struct of fi prototypes
T_ex = T;

args = {x_type, NPol_ex, NTaps_ex, w_type, T_ex};

fprintf('  Input types defined.\n');

%% 3 — Generate MEX
fprintf('\n=== Step 3: Generating MEX via codegen ===\n');

cfg = coder.config('mex');
cfg.GenerateReport = true;
cfg.EnableMexProfiling = false;

codegen('-config', cfg, ...
        'cr_viterbiViterbi_fxp', ...
        '-args', args, ...
        '-o', fullfile(srcDir, 'cr_viterbiViterbi_fxp_mex'));

fprintf('\n=== MEX build complete ===\n');
fprintf('Output: %s\n', fullfile(srcDir, ['cr_viterbiViterbi_fxp_mex.' mexext]));
fprintf('\nCall it exactly like the MATLAB version:\n');
fprintf('  v = cr_viterbiViterbi_fxp_mex(x, NPol, NTaps, VVFilter, T)\n');
