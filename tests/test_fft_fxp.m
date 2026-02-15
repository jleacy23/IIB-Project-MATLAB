classdef test_fft_fxp < matlab.unittest.TestCase
    %TEST_FFT_FXP  Verify fft_fxp (MATLAB + MEX) against built-in fft/ifft.
    %
    %   Run all tests:
    %     results = runtests('test_fft_fxp');
    %
    %   The MEX tests require that build_fft_fxp_mex has been run first
    %   with the same N_FFT (256).

    properties (Constant)
        N = 256           % FFT size — must match build_fft_fxp_mex.m
    end

    methods (TestMethodSetup)
        function seedRng(~)
            rng(42);
        end
    end

    % ================================================================
    %  Double-precision tests (floating-point baseline)
    % ================================================================
    methods (Test)

        function testFFT_impulse_double(testCase)
            %  FFT of an impulse [1; 0; …; 0] should be all ones.
            N = testCase.N;
            T = fft_fxp_types('double');
            x = zeros(N, 1);
            x(1) = 1;
            X     = fft_fxp(x, false, false, T);
            X_ref = fft(x);
            testCase.verifyEqual(X, X_ref, 'AbsTol', 1e-12, ...
                'FFT of impulse must equal built-in fft.');
        end

        function testFFT_random_double(testCase)
            N = testCase.N;
            T = fft_fxp_types('double');
            x     = randn(N, 1) + 1j*randn(N, 1);
            X     = fft_fxp(x, false, false, T);
            X_ref = fft(x);
            testCase.verifyEqual(X, X_ref, 'AbsTol', 1e-10, ...
                'FFT of random complex vector must match built-in fft.');
        end

        function testIFFT_random_double(testCase)
            N = testCase.N;
            T = fft_fxp_types('double');
            X     = randn(N, 1) + 1j*randn(N, 1);
            x     = fft_fxp(X, true, false, T);
            x_ref = ifft(X);
            testCase.verifyEqual(x, x_ref, 'AbsTol', 1e-10, ...
                'IFFT of random complex vector must match built-in ifft.');
        end

        function testRoundtrip_double(testCase)
            N = testCase.N;
            T = fft_fxp_types('double');
            x     = randn(N, 1) + 1j*randn(N, 1);
            X     = fft_fxp(x, false, false, T);
            x_rec = fft_fxp(X, true, false, T);
            testCase.verifyEqual(x_rec, x, 'AbsTol', 1e-10, ...
                'FFT -> IFFT roundtrip must recover the input.');
        end

        % ============================================================
        %  Fixed-point tests (32-bit)
        % ============================================================

        function testFFT_fxp32(testCase)
            N = testCase.N;
            T = fft_fxp_types('fixed32');
            x     = randn(N, 1) + 1j*randn(N, 1);
            x_fi  = cast(x, 'like', T.x);
            X     = fft_fxp(x_fi, false, false, T);
            X_ref = fft(double(x_fi));

            nrmse = norm(double(X) - X_ref) / norm(X_ref);
            testCase.verifyLessThan(nrmse, 0.01, ...
                'Fixed-point 32-bit FFT NRMSE must be < 1%.');
        end

        function testIFFT_fxp32(testCase)
            N = testCase.N;
            T = fft_fxp_types('fixed32');
            X     = randn(N, 1) + 1j*randn(N, 1);
            X_fi  = cast(X, 'like', T.x);
            x     = fft_fxp(X_fi, true, false, T);
            x_ref = ifft(double(X_fi));

            nrmse = norm(double(x) - x_ref) / norm(x_ref);
            testCase.verifyLessThan(nrmse, 0.01, ...
                'Fixed-point 32-bit IFFT NRMSE must be < 1%.');
        end

        function testRoundtrip_fxp32(testCase)
            N = testCase.N;
            T = fft_fxp_types('fixed32');
            x     = randn(N, 1) + 1j*randn(N, 1);
            x_fi  = cast(x, 'like', T.x);
            X     = fft_fxp(x_fi, false, false, T);
            x_rec = fft_fxp(X, true, false, T);

            nrmse = norm(double(x_rec) - double(x_fi)) / norm(double(x_fi));
            testCase.verifyLessThan(nrmse, 0.02, ...
                'Fixed-point 32-bit FFT->IFFT roundtrip NRMSE must be < 2%.');
        end

        % ============================================================
        %  Power-of-2 twiddle tests
        % ============================================================

        function testPo2Twiddle_double(testCase)
            N = testCase.N;
            T = fft_fxp_types('double');
            x       = randn(N, 1) + 1j*randn(N, 1);
            X_exact = fft_fxp(x, false, false, T);
            X_po2   = fft_fxp(x, false, true, T);

            relErr = norm(X_po2 - X_exact) / norm(X_exact);
            testCase.verifyLessThan(relErr, 1.0, ...
                'Power-of-2 twiddle FFT must be a reasonable approximation.');
        end

        function testPo2Twiddle_fxp32(testCase)
            N = testCase.N;
            T = fft_fxp_types('fixed32');
            x    = randn(N, 1) + 1j*randn(N, 1);
            x_fi = cast(x, 'like', T.x);
            X_exact = fft_fxp(x_fi, false, false, T);
            X_po2   = fft_fxp(x_fi, false, true, T);

            relErr = norm(double(X_po2) - double(X_exact)) / ...
                     norm(double(X_exact));
            testCase.verifyLessThan(relErr, 1.0, ...
                'Po2 twiddle fxp32 FFT must be a reasonable approximation.');
        end

        % ============================================================
        %  MEX tests (bit-exact comparison with MATLAB fixed-point)
        % ============================================================

        function testFFT_mex(testCase)
            assumeMexAvailable(testCase);
            N = testCase.N;
            T = fft_fxp_types('fixed32');
            x    = randn(N, 1) + 1j*randn(N, 1);
            x_fi = cast(x, 'like', T.x);

            X_ml  = fft_fxp(x_fi, false, false, T);
            X_mex = fft_fxp_mex(x_fi, false, false, T);

            testCase.verifyEqual(double(X_mex), double(X_ml), ...
                'MEX FFT output must be bit-exact with MATLAB fxp.');
        end

        function testIFFT_mex(testCase)
            assumeMexAvailable(testCase);
            N = testCase.N;
            T = fft_fxp_types('fixed32');
            X    = randn(N, 1) + 1j*randn(N, 1);
            X_fi = cast(X, 'like', T.x);

            x_ml  = fft_fxp(X_fi, true, false, T);
            x_mex = fft_fxp_mex(X_fi, true, false, T);

            testCase.verifyEqual(double(x_mex), double(x_ml), ...
                'MEX IFFT output must be bit-exact with MATLAB fxp.');
        end

        function testPo2Twiddle_mex(testCase)
            assumeMexAvailable(testCase);
            N = testCase.N;
            T = fft_fxp_types('fixed32');
            x    = randn(N, 1) + 1j*randn(N, 1);
            x_fi = cast(x, 'like', T.x);

            X_ml  = fft_fxp(x_fi, false, true, T);
            X_mex = fft_fxp_mex(x_fi, false, true, T);

            testCase.verifyEqual(double(X_mex), double(X_ml), ...
                'MEX po2-twiddle FFT must be bit-exact with MATLAB fxp.');
        end

        function testFFT_mex_vs_builtin(testCase)
            %  Verify MEX FFT output against built-in fft (with tolerance).
            assumeMexAvailable(testCase);
            N = testCase.N;
            T = fft_fxp_types('fixed32');
            x    = randn(N, 1) + 1j*randn(N, 1);
            x_fi = cast(x, 'like', T.x);

            X_mex = fft_fxp_mex(x_fi, false, false, T);
            X_ref = fft(double(x_fi));

            nrmse = norm(double(X_mex) - X_ref) / norm(X_ref);
            testCase.verifyLessThan(nrmse, 0.01, ...
                'MEX FFT NRMSE vs built-in fft must be < 1%.');
        end

        function testIFFT_mex_vs_builtin(testCase)
            %  Verify MEX IFFT output against built-in ifft (with tolerance).
            assumeMexAvailable(testCase);
            N = testCase.N;
            T = fft_fxp_types('fixed32');
            X    = randn(N, 1) + 1j*randn(N, 1);
            X_fi = cast(X, 'like', T.x);

            x_mex = fft_fxp_mex(X_fi, true, false, T);
            x_ref = ifft(double(X_fi));

            nrmse = norm(double(x_mex) - x_ref) / norm(x_ref);
            testCase.verifyLessThan(nrmse, 0.01, ...
                'MEX IFFT NRMSE vs built-in ifft must be < 1%.');
        end

    end

    % ================================================================
    %  Helpers
    % ================================================================
    methods (Access = private)
        function assumeMexAvailable(testCase)
            testCase.assumeTrue(exist('fft_fxp_mex', 'file') == 3, ...
                'fft_fxp_mex not found — run build_fft_fxp_mex first.');
        end
    end
end
