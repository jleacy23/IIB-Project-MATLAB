classdef test_Normalise < matlab.unittest.TestCase
    % Tests for modem.normalise.
    % Generates a DP-QPSK signal, applies AWGN and chromatic dispersion,
    % normalises the received signal into the unit box, and plots the
    % before/after constellations for visual inspection.

    properties (Constant)
        % Simulation parameters
        N_pol   = 2
        Ns      = 4096          % symbols per polarization
        SpS     = 2             % samples per symbol
        Rs      = 32            % symbol rate [GBd]
        L       = 80            % fibre length [km]
        SNR_dB  = 25            % SNR [dB]
        D       = 17            % dispersion coeff [ps/nm/km]
        CWL     = 1550          % central wavelength [nm]
        Pct     = 95            % normalisation percentile
    end

    methods (TestMethodSetup)
        function seedRng(~)
            rng(42);            % reproducible noise draws
        end
    end

    % ----- helper ---------------------------------------------------
    methods (Access = private)
        function [rxSym, txSym] = generateRx(testCase)
            % Build a DP-QPSK signal, apply AWGN + CD, return symbol-rate
            % samples of the impaired (un-normalised) signal.
            Nbits   = 4 * testCase.Ns;
            bits    = modem.randomBits(Nbits);
            txSym   = modem.modulate(bits);
            txSig   = modem.rectPulse(txSym, testCase.SpS);

            rxSig = channel.add_awgn(txSig, testCase.SNR_dB);
            rxSig = channel.add_chromatic_dispersion(rxSig, ...
                testCase.L, testCase.SpS, testCase.Rs, testCase.D, testCase.CWL);

            rxSym = rxSig(1:testCase.SpS:end, :);   % downsample to symbol rate
        end
    end

    methods (Test)
        % ====== Behavioural checks ==================================
        function testOutputShape(testCase)
            rxSym = generateRx(testCase);
            [normSym, scale] = modem.normalise(rxSym, testCase.Pct);
            testCase.verifySize(normSym, size(rxSym), ...
                'Normalised output size must match input size.');
            testCase.verifySize(scale, [1, testCase.N_pol], ...
                'Scale must be a 1 x N_pol row vector.');
        end

        function testWithinUnitBox(testCase)
            rxSym = generateRx(testCase);
            normSym = modem.normalise(rxSym, testCase.Pct);
            testCase.verifyLessThanOrEqual(abs(real(normSym(:))), 1, ...
                'Real part must lie within [-1, 1].');
            testCase.verifyLessThanOrEqual(abs(imag(normSym(:))), 1, ...
                'Imag part must lie within [-1, 1].');
        end

        function testFinite(testCase)
            rxSym = generateRx(testCase);
            normSym = modem.normalise(rxSym, testCase.Pct);
            testCase.verifyTrue(all(isfinite(normSym(:))), ...
                'Normalised signal produced NaN/Inf values.');
        end

        % ====== Visual constellation test ===========================
        function testVisualConstellation(testCase)
            [rxSym, txSym] = generateRx(testCase);
            normSym = modem.normalise(rxSym, testCase.Pct);

            figure('Name','Normalise Visual Test','Position',[100 100 1000 700]);
            stages = {txSym, rxSym, normSym};
            labels = {'Tx (ideal)', 'Rx: AWGN + CD', 'Normalised'};

            for k = 1:numel(stages)
                for p = 1:testCase.N_pol
                    idx = (k-1)*testCase.N_pol + p;
                    subplot(3, testCase.N_pol, idx);
                    plot(real(stages{k}(:,p)), imag(stages{k}(:,p)), '.', ...
                         'MarkerSize', 2);
                    grid on; axis equal;
                    title(sprintf('%s  –  Pol %d', labels{k}, p));
                    xlabel('In-Phase'); ylabel('Quadrature');
                end
            end
            sgtitle('Constellation before/after normalisation (AWGN + CD)');

            testCase.verifyTrue(all(isfinite(normSym(:))), ...
                'Normalised signal produced NaN/Inf values.');
        end
    end
end
