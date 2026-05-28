function build_diag_mex(FL, mode)
%BUILD_DIAG_MEX  Build the fixed-point MEX needed for pipeline diagnosis.
%
%   build_diag_mex                 % FL=16, mode='all'
%   build_diag_mex(FL)             % arbitrary fractional length, mode='all'
%   build_diag_mex(FL, 'verify')   % FAST: only the 2 MEX needed to confirm
%                                  %       the clk fix (eq<FL> + eqC16)
%
%   mode = 'all' (default) produces (in src/+...):
%     +eq_clk/combined_cd_fd_gardner_adaptive_fxp_eq<FL>_mex
%     +freq_recovery/differential_kay_fxp_fc<FL>_mex
%     +carrier_recovery/pilots_only_fxp_fc<FL>_mex
%   plus four mixed-precision equaliser MEX used by diag_eq_localize to pin
%   down which sub-stage (Static / Clk / AdaptEq) breaks at FL=16.  Each sets
%   two sub-stages near-lossless (WL=48/FL=32) and one at WL=32/FL=16:
%     ..._eqHI_mex (all lossless)  ..._eqS16_mex (Static FL16)
%     ..._eqC16_mex (Clk FL16)     ..._eqA16_mex (AdaptEq FL16)
%   (the all-FL16 case is just eq<FL>_mex.)  The HI/lo configs here MUST
%   match those hard-coded in diag_eq_localize.
%
%   mode = 'verify' builds ONLY eq<FL>_mex (all-FL16) and eqC16_mex
%   (Clk=FL16, others lossless) -- the two diag_eq_localize rows that prove
%   the recovery_fxp loop-filter fix.  FR/CR and the other eq variants are
%   left as-is (rebuild with 'all' for a full localisation sweep).
%
%   Run this locally (codegen needs the toolchain), then the diagnostic
%   scripts / the sweep can call the freshly-built MEX with no rebuild.

    if nargin < 1 || isempty(FL),   FL = 16;       end
    if nargin < 2 || isempty(mode), mode = 'all';  end
    verifyOnly = strcmpi(mode, 'verify');

    here     = fileparts(mfilename('fullpath'));
    repoRoot = fileparts(fileparts(here));
    addpath(genpath(fullfile(repoRoot, 'src')));
    addpath(fullfile(repoRoot, 'build'));

    % Reuse the sweep's own parameter extraction so the build config here
    % is identical to the sweep's.  Constructing the TestCase does not run
    % TestClassSetup (that only happens under a test runner).
    P = pipeline_fxp_sweep.extractParams(pipeline_fxp_sweep);

    cfgCoder = coder.config('mex');
    cfgCoder.GenerateReport = false;

    hi = struct('WL', 48, 'FL', 32);              % near-lossless
    lo = struct('WL', P.IntBits + FL, 'FL', FL);  % the stage under test

    fprintf('\n=== build_diag_mex: FL=%d, mode=%s ===\n', FL, lower(mode));
    tAll = tic;

    % --- All-FL16 equaliser (the sweep's own binary; picks up clk fix) ---
    mexEq = pipeline_fxp_sweep.mexEqName(FL);
    fprintf('  %s ... ', mexEq); t=tic;
    pipeline_fxp_sweep.buildEqMex(P, cfgCoder, FL, mexEq, repoRoot);
    fprintf('(%.0fs)\n', toc(t));

    % --- Clk=FL16 variant (localiser row that isolates clock recovery) ---
    fprintf('  combined_cd_fd_gardner_adaptive_fxp_eqC16_mex ... '); t=tic;
    buildEqVariant(P, cfgCoder, struct('Static',hi,'Clk',lo,'AdaptEq',hi), ...
        'combined_cd_fd_gardner_adaptive_fxp_eqC16_mex', repoRoot);
    fprintf('(%.0fs)\n', toc(t));

    if verifyOnly
        fprintf('=== done in %.1fs (verify) ===\n', toc(tAll));
        fprintf('Built: %s, eqC16.  Run diag_eq_localize to read the two rows.\n', mexEq);
        return;
    end

    % --- FR + CR (full mode only) ---------------------------------------
    mexFR = pipeline_fxp_sweep.mexFRName('differential_kay', FL);
    fprintf('  %s ... ', mexFR); t=tic;
    pipeline_fxp_sweep.buildFRMex(P, cfgCoder, 'differential_kay', FL, mexFR, repoRoot);
    fprintf('(%.0fs)\n', toc(t));

    mexCR = pipeline_fxp_sweep.mexCRName(FL);
    fprintf('  %s ... ', mexCR); t=tic;
    pipeline_fxp_sweep.buildCRMex(P, cfgCoder, FL, mexCR, repoRoot);
    fprintf('(%.0fs)\n', toc(t));

    % --- Remaining mixed-precision eq variants (eqC16 already built) -----
    variants = {
        'eqHI',  hi, hi, hi
        'eqS16', lo, hi, hi
        'eqA16', hi, hi, lo };
    for vi = 1:size(variants,1)
        mexBase = ['combined_cd_fd_gardner_adaptive_fxp_' variants{vi,1} '_mex'];
        fprintf('  %s ... ', mexBase); t=tic;
        buildEqVariant(P, cfgCoder, ...
            struct('Static',variants{vi,2}, 'Clk',variants{vi,3}, 'AdaptEq',variants{vi,4}), ...
            mexBase, repoRoot);
        fprintf('(%.0fs)\n', toc(t));
    end

    fprintf('=== done in %.1fs ===\n', toc(tAll));
    fprintf('Built base trio:\n  %s\n  %s\n  %s\n', mexEq, mexFR, mexCR);
    fprintf('Built eq variants: eqHI, eqS16, eqC16, eqA16\n');
end

% ---------------------------------------------------------------------
function buildEqVariant(P, cfgCoder, fxpCfg, mexBase, repoRoot)
%   Build the combined equaliser MEX with an arbitrary per-stage fxp
%   config, then rename the generic output to mexBase.  Mirrors the B
%   struct that pipeline_fxp_sweep.buildEqMex assembles, but with a
%   caller-supplied composite FxpConfig_CombGardner.
    B = struct();
    B.FxpConfig_CombGardner = fxpCfg;
    B.SpS        = P.SpS;
    B.NFFT       = P.NFFT;
    B.NOverlap   = 2 * ceil((P.NCD - 1) / 2);
    B.D          = P.D;
    B.L          = P.L_km;
    B.CWL        = P.CWL;
    B.Rs         = P.Rs;
    B.Rolloff    = P.Rolloff;
    B.N_pol      = P.N_pol;
    B.po2Twiddle = false;
    B.cfoEnable  = false;
    B.Ns         = P.N_sub_target * P.SUBFRAME_SYMS;
    B.CR_ki      = P.ki_gardner_po2_off;
    B.CR_kp      = P.kp_gardner_po2_off;
    B.CR_NLanes  = P.NLanesGard;
    B.AEQ_NTaps          = P.NTapsAEQ;
    B.AEQ_Mu             = P.MuAEQ;
    B.AEQ_SingleSpike    = P.SingleSpike;
    B.AEQ_N1             = P.N1AEQ;
    B.AEQ_NOut           = P.NOutAEQ;
    B.AEQ_SignOnly       = P.SignOnly;
    B.AEQ_UpdateStep     = 1;
    B.AEQ_PLanes         = P.PLanesAEQ;
    B.AEQ_Mode           = 0;
    B.AEQ_BlockLen       = P.PLanesAEQ;
    B.AEQ_SubframeBlocks = 0;

    build_eq_clk_combined_cd_fd_gardner_adaptive_fxp_mex(B, cfgCoder);

    srcDir  = fullfile(repoRoot, 'src', '+eq_clk');
    ext     = ['.' mexext];
    srcFile = fullfile(srcDir, ['combined_cd_fd_gardner_adaptive_fxp_mex' ext]);
    dstFile = fullfile(srcDir, [mexBase ext]);
    if isfile(srcFile)
        if isfile(dstFile), delete(dstFile); end
        movefile(srcFile, dstFile);
    end
end
