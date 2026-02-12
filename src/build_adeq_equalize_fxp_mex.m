%% build_adeq_equalize_fxp_mex.m
%  Build script: checks codegen readiness and compiles a MEX function from
%  adeq_equalize_fxp using MATLAB Coder + Fixed-Point Designer.
%
%  Prerequisites: MATLAB Coder, Fixed-Point Designer
%
%  Usage:  run this script from the project root or from src/.

%% 0 — Ensure src/ is on the MATLAB path
srcDir = fullfile(fileparts(mfilename('fullpath')));
addpath(srcDir);

%% 1 — Static code-generation readiness check
fprintf('\n=== Step 1: Code generation readiness (coder.screener) ===\n');
coder.screener('adeq_equalize_fxp');
fprintf('  coder.screener completed — review any flagged issues above.\n');

%% 2 — Define example / prototype input types
fprintf('\n=== Step 2: Defining input argument types ===\n');

% Fixed-point types table (all fi prototypes live here)
T = adeq_equalize_fxp_types('fixed32');

% x  — complex fi input signal, variable-length rows, 2 columns
%      Derive the prototype from T.x so WL, FL, and fimath stay in sync.
x_proto = fi(complex(0, 0), numerictype(T.x), fimath(T.x));
x_type  = coder.typeof(x_proto, [Inf, 2], [true, false]);

% Scalar parameters (all double — MATLAB's natural numeric type)
SpS_ex         = double(2);
NTaps_ex       = double(15);
Mu_ex          = double(1e-3);
N1_ex          = double(2000);
N2_ex          = double(4000);   % pass 0 when CMA-only (no switch)
NOut_ex        = double(500);
SingleSpike_ex = true;

% Eq — char vector, variable-length up to 7 chars ('CMA+RDE')
Eq_type = coder.typeof('a', [1, 7], [false, true]);

% T  — struct of fi prototypes (pass in the table directly)
T_ex = T;

args = {x_type, SpS_ex, Eq_type, NTaps_ex, Mu_ex, ...
        SingleSpike_ex, N1_ex, N2_ex, NOut_ex, T_ex};

fprintf('  Input types defined.\n');

%% 3 — Generate MEX
fprintf('\n=== Step 3: Generating MEX via codegen ===\n');

cfg = coder.config('mex');
cfg.GenerateReport = true;           % produce an HTML report
cfg.EnableMexProfiling = false;

codegen('-config', cfg, ...
        'adeq_equalize_fxp', ...
        '-args', args, ...
        '-o', fullfile(srcDir, 'adeq_equalize_fxp_mex'));

fprintf('\n=== MEX build complete ===\n');
fprintf('Output: %s\n', fullfile(srcDir, ['adeq_equalize_fxp_mex.' mexext]));
fprintf('\nCall it exactly like the MATLAB version:\n');
fprintf('  y = adeq_equalize_fxp_mex(x, SpS, Eq, NTaps, Mu, SingleSpike, N1, N2, NOut, T)\n');
fprintf('  (N2 must be a scalar — use 0 when no CMA→RDE switch is needed)\n');
