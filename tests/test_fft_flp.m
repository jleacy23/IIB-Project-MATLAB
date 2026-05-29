classdef test_fft_flp < matlab.unittest.TestCase
    %TEST_FFT_FLP  Verify fft_flp against MATLAB's built-in fft/ifft.
    %
    %   Run all tests:
    %     results = runtests('test_fft_flp');
    %
    %   Normalisation convention: fft_flp spreads the 1/N normalisation
    %   across the forward transform (1/2 per radix-2 stage) and applies no
    %   scaling on the inverse.  So relative to MATLAB's built-ins:
    %       fft_flp(x, false) == fft(x) / N
    %       fft_flp(X, true ) == ifft(X) * N
    %   The forward<->inverse round trip is therefore still the identity.

    properties (Constant)
        N = 512           % FFT size (power of 2)
    end

    methods (TestMethodSetup)
        function seedRng(~)
            rng(42);
        end
    end

    % ================================================================
    %  Normalisation / correctness vs built-ins
    % ================================================================
    methods (Test)

        function testFFT_impulse(testCase)
            %  FFT of an impulse [1; 0; …; 0] is all ones, here 1/N-scaled.
            N = testCase.N;
            x = zeros(N, 1);
            x(1) = 1;
            X     = fft.fft_flp(x, false);
            X_ref = fft(x) / N;     % forward FFT is 1/N-normalised
            testCase.verifyEqual(X, X_ref, 'AbsTol', 1e-12, ...
                'FFT of impulse must equal built-in fft / N.');
        end

        function testFFT_random(testCase)
            N = testCase.N;
            x     = randn(N, 1) + 1j*randn(N, 1);
            X     = fft.fft_flp(x, false);
            X_ref = fft(x) / N;     % forward FFT is 1/N-normalised
            testCase.verifyEqual(X, X_ref, 'AbsTol', 1e-10, ...
                'FFT of random complex vector must match built-in fft / N.');
        end

        function testFFT_realInput(testCase)
            %  Real-valued input should still match fft / N.
            N = testCase.N;
            x     = randn(N, 1);
            X     = fft.fft_flp(x, false);
            X_ref = fft(x) / N;
            testCase.verifyEqual(X, X_ref, 'AbsTol', 1e-10, ...
                'FFT of real vector must match built-in fft / N.');
        end

        function testIFFT_random(testCase)
            N = testCase.N;
            X     = randn(N, 1) + 1j*randn(N, 1);
            x     = fft.fft_flp(X, true);
            x_ref = ifft(X) * N;    % inverse FFT applies no 1/N scaling
            testCase.verifyEqual(x, x_ref, 'AbsTol', 1e-10, ...
                'IFFT of random complex vector must match built-in ifft * N.');
        end

        function testRoundtrip(testCase)
            N = testCase.N;
            x     = randn(N, 1) + 1j*randn(N, 1);
            X     = fft.fft_flp(x, false);
            x_rec = fft.fft_flp(X, true);
            testCase.verifyEqual(x_rec, x, 'AbsTol', 1e-10, ...
                'FFT -> IFFT roundtrip must recover the input.');
        end

        function testDefaultArgs(testCase)
            %  Omitting inverse/po2Twiddle must default to a forward,
            %  exact-twiddle transform.
            N = testCase.N;
            x = randn(N, 1) + 1j*randn(N, 1);
            testCase.verifyEqual(fft.fft_flp(x),            ...
                                 fft.fft_flp(x, false, false), ...
                'AbsTol', 0, 'Default args must equal forward exact FFT.');
            testCase.verifyEqual(fft.fft_flp(x, []),        ...
                                 fft.fft_flp(x, false),     ...
                'AbsTol', 0, 'Empty inverse arg must default to forward.');
        end

        % ============================================================
        %  Multi-column / N-D input (transform along dim 1)
        % ============================================================

        function testFFT_matrixColumns(testCase)
            %  Each column transformed independently (matches fft).
            N = testCase.N;
            M = 4;
            x     = randn(N, M) + 1j*randn(N, M);
            X     = fft.fft_flp(x, false);
            X_ref = fft(x) / N;
            testCase.verifyEqual(X, X_ref, 'AbsTol', 1e-10, ...
                'Column-wise FFT must match built-in fft / N.');
        end

        function testFFT_ndArray(testCase)
            %  Trailing dimensions are processed independently.
            N = testCase.N;
            x     = randn(N, 2, 3) + 1j*randn(N, 2, 3);
            X     = fft.fft_flp(x, false);
            X_ref = fft(x) / N;
            testCase.verifyEqual(size(X), size(x), ...
                'Output must preserve input size.');
            testCase.verifyEqual(X, X_ref, 'AbsTol', 1e-10, ...
                'N-D FFT must match built-in fft / N along dim 1.');
        end

        function testMatrixRoundtrip(testCase)
            N = testCase.N;
            x     = randn(N, 3) + 1j*randn(N, 3);
            X     = fft.fft_flp(x, false);
            x_rec = fft.fft_flp(X, true);
            testCase.verifyEqual(x_rec, x, 'AbsTol', 1e-10, ...
                'Column-wise FFT -> IFFT roundtrip must recover input.');
        end

        % ============================================================
        %  Power-of-2 twiddle approximation
        % ============================================================

        function testPo2Twiddle_approx(testCase)
            %  Po2 twiddles are a coarse but bounded approximation.
            N = testCase.N;
            x       = randn(N, 1) + 1j*randn(N, 1);
            X_exact = fft.fft_flp(x, false, false);
            X_po2   = fft.fft_flp(x, false, true);

            relErr = norm(X_po2 - X_exact) / norm(X_exact);
            testCase.verifyLessThan(relErr, 1.0, ...
                'Power-of-2 twiddle FFT must be a reasonable approximation.');
        end

        % ============================================================
        %  Small sizes / edge cases
        % ============================================================

        function testFFT_N2(testCase)
            %  Smallest valid transform.
            x     = [3 + 1j; -2 + 4j];
            X     = fft.fft_flp(x, false);
            X_ref = fft(x) / 2;
            testCase.verifyEqual(X, X_ref, 'AbsTol', 1e-12, ...
                'N=2 FFT must match built-in fft / N.');
        end

        % ============================================================
        %  Error handling
        % ============================================================

        function testError_nonPo2Length(testCase)
            x = randn(48, 1);    % not a power of two
            testCase.verifyError(@() fft.fft_flp(x, false), ...
                'fft_flp:lenNotPo2', ...
                'Non-power-of-2 length must raise fft_flp:lenNotPo2.');
        end

        function testError_scalarInput(testCase)
            %  N = 1 is below the minimum (N >= 2 required).
            testCase.verifyError(@() fft.fft_flp(5, false), ...
                'fft_flp:lenNotPo2', ...
                'Scalar (N=1) input must raise fft_flp:lenNotPo2.');
        end

    end
end
