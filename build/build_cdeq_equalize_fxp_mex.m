%% build_cdeq_equalize_fxp_mex.m
%  Build script: checks codegen readiness and compiles a MEX function from
%  cdeq_equalize_fxp using MATLAB Coder + Fixed-Point Designer.
%
%  Prerequisites: MATLAB Coder, Fixed-Point Designer
%
%  Usage:  run this script from the project root or from build/.

%% 0 — Ensure src/ is on the MATLAB path
buildDir = fileparts(mfilename('fullpath'));
srcDir   = fullfile(buildDir, '..', 'src');
addpath(srcDir);

%% 1 — Static code-generation readiness check
fprintf('\n=== Step 1: Code generation readiness (coder.screener) ===\n');
coder.screener('cdeq_equalize_fxp');
fprintf('  coder.screener completed — review any flagged issues above.\n');

%% 2 — Define example / prototype input types
fprintf('\n=== Step 2: Defining input argument types ===\n');

% Fixed-point type configuration
T = cdeq_equalize_fxp_types('fixed32');

% In — complex fi input signal, variable-length rows, 2 columns
x_proto = fi(complex(0, 0), numerictype(T.x), fimath(T.x));
In_type = coder.typeof(x_proto, [Inf, 2], [true, false]);

% Scalar parameters (all double)
D_ex       = double(17);       % dispersion [ps/(nm*km)]
L_ex       = double(80);       % fibre length [km]
CLambda_ex = double(1550);     % central wavelength [nm]
Rs_ex      = double(32);       % symbol rate [GBd]
NPol_ex    = double(2);        % number of polarizations
SpSIn_ex   = double(2);        % samples per symbol
NFFT_ex    = double(512);      % FFT block size

% Logical flag
po2Twiddle_ex = false;

% T — struct of fi prototypes
T_ex = T;

args = {In_type, D_ex, L_ex, CLambda_ex, Rs_ex, NPol_ex, SpSIn_ex, ...
        NFFT_ex, po2Twiddle_ex, T_ex};

fprintf('  Input types defined (NFFT = %d, fixed32).\n', NFFT_ex);

%% 3 — Generate MEX
fprintf('\n=== Step 3: Generating MEX via codegen ===\n');

cfg = coder.config('mex');
cfg.GenerateReport = true;           % produce an HTML report
cfg.EnableMexProfiling = false;

codegen('-config', cfg, ...
        'cdeq_equalize_fxp', ...
        '-args', args, ...
        '-o', fullfile(srcDir, 'cdeq_equalize_fxp_mex'));

fprintf('\n=== MEX build complete ===\n');
fprintf('Output: %s\n', fullfile(srcDir, ['cdeq_equalize_fxp_mex.' mexext]));
fprintf('\nCall it exactly like the MATLAB version:\n');
fprintf('  Out = cdeq_equalize_fxp_mex(In, D, L, CLambda, Rs, NPol, SpSIn, NFFT, po2Twiddle, T)\n');
