%% build_fft_fxp_mex.m
%  Build script: checks codegen readiness and compiles a MEX function from
%  fft_fxp using MATLAB Coder + Fixed-Point Designer.
%
%  Prerequisites: MATLAB Coder, Fixed-Point Designer
%
%  Usage:  run this script from the project root or from build/.
%
%  The compiled MEX accepts a fixed-size complex column vector of length
%  N_FFT (defined below).  Change N_FFT to match your application; it
%  must be a power of 2.

%% 0 — Ensure src/ is on the MATLAB path
buildDir = fileparts(mfilename('fullpath'));
srcDir   = fullfile(buildDir, '..', 'src');
addpath(srcDir);

%% 1 — Static code-generation readiness check
fprintf('\n=== Step 1: Code generation readiness (coder.screener) ===\n');
coder.screener('fft_fxp');
fprintf('  coder.screener completed — review any flagged issues above.\n');

%% 2 — Define example / prototype input types
fprintf('\n=== Step 2: Defining input argument types ===\n');

% FFT size (must be a power of 2) — change as needed
N_FFT = 256;

% Fixed-point type configuration
T = fft_fxp_types('fixed32');

% x — complex fi column vector, FIXED size [N_FFT x 1]
x_proto = fi(complex(0, 0), numerictype(T.x), fimath(T.x));
x_type  = coder.typeof(x_proto, [N_FFT, 1], [false, false]);

% Scalar flags (logical)
inverse_ex    = false;
po2Twiddle_ex = false;

% T — struct of fi prototypes
T_ex = T;

args = {x_type, inverse_ex, po2Twiddle_ex, T_ex};

fprintf('  Input types defined (N_FFT = %d, fixed32).\n', N_FFT);

%% 3 — Generate MEX
fprintf('\n=== Step 3: Generating MEX via codegen ===\n');

cfg = coder.config('mex');
cfg.GenerateReport = true;           % produce an HTML report
cfg.EnableMexProfiling = false;

codegen('-config', cfg, ...
        'fft_fxp', ...
        '-args', args, ...
        '-o', fullfile(srcDir, 'fft_fxp_mex'));

fprintf('\n=== MEX build complete ===\n');
fprintf('Output: %s\n', fullfile(srcDir, ['fft_fxp_mex.' mexext]));
fprintf('\nCall it exactly like the MATLAB version:\n');
fprintf('  X = fft_fxp_mex(x, inverse, po2Twiddle, T)\n');
fprintf('  x must be a complex fi column vector of length %d.\n', N_FFT);
