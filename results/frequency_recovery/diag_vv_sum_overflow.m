function diag_vv_sum_overflow()
%DIAG_VV_SUM_OVERFLOW  Size the Viterbi-Viterbi 4th-power block sum.
%
%   The fixed-point Viterbi-Viterbi estimator accumulates, per block,
%
%       sum4(blk) = sum_{i in block} sum_{k=1..L_filt} w(k) * s(i+k-halfL)^4
%
%   in the T.theta type before taking its angle with CORDIC.  This script
%   measures how large that accumulator actually gets when the input signal
%   is constrained to the box [-1,1] x [-1,1] (|Re|,|Im| <= 1), so we can
%   decide whether the swept types (integer part = IntBits, here 16, so 15
%   magnitude bits) need overflow handling (saturation / a wider accumulator)
%   or whether the default Wrap-on-overflow is safe.
%
%   Three data populations are tested, all clipped/scaled into [-1,1]:
%     1. Worst-case aligned : every s^4 forced to the same phase (rounded
%        QPSK at the box corners) — the algebraic upper bound.
%     2. Uniform box        : Re,Im ~ U(-1,1) i.i.d. (incoherent stress).
%     3. Realistic channel  : modem.modulate -> CFO -> AWGN -> phase noise,
%        then scaled so the peak |Re|/|Im| sample is 1.
%
%   Run:  diag_vv_sum_overflow

    here = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(here, '..', '..', 'src')));
    rng(42);

    %% ---- Parameters (mirror bit_width_full) ------------------------
    Rs        = 30.5;
    N_pol     = 2;
    NTaps     = 10;
    BlockLen  = 32;
    LW_Hz     = 1000e3;
    DeltaF_Hz = 2e9;
    SNR_dB    = 0:1:30;         % full sweep, to find the worst VV filter
    IntBits   = 16;            % swept types: WL = IntBits + FL
    AvailIntMagBits = IntBits - 1;     % sign takes one bit -> 15 magnitude bits
    AvailRange      = 2^AvailIntMagBits;
    K_CORDIC  = prod(sqrt(1 + 2.^(-2*(0:31))));   % vectoring gain ~1.6468

    L_filt = 2 * NTaps + 1;

    %% ---- VV filter: pick the one with the largest sum|w| -----------
    BITS_PER_SF = 3586 * 2 * 2;
    [tmp, ~, ~, ~] = modem.modulate(modem.randomBits(BITS_PER_SF));
    symEnergy = mean(abs(tmp(:)).^2);

    sumAbsW = zeros(numel(SNR_dB), 1);
    wAll    = cell(numel(SNR_dB), 1);
    for si = 1:numel(SNR_dB)
        w = carrier_recovery.genVVFilter(LW_Hz, Rs, SNR_dB(si), symEnergy, N_pol, NTaps);
        wAll{si}    = w(:);
        sumAbsW(si) = sum(abs(w));
    end
    [maxSumAbsW, wi] = max(sumAbsW);
    w_worst = wAll{wi};
    fprintf('VV filter taps L_filt = %d, sum(w) = 1 by construction.\n', L_filt);
    fprintf('max over SNR of sum|w| = %.4f (at SNR = %d dB)\n\n', maxSumAbsW, SNR_dB(wi));

    %% ---- Algebraic bound -------------------------------------------
    %   |s|<=sqrt(2) => |s^4| <= 4; |Re(s^4)|,|Im(s^4)| <= 4.
    %   |sum4| <= BlockLen * sum|w| * 4.
    boundAbs = BlockLen * maxSumAbsW * 4;
    fprintf('--- Algebraic upper bound (|Re|,|Im|<=1) ---\n');
    fprintf('  max|s^4|                 = %.3f\n', 4);
    fprintf('  worst-case |sum4|        = %.2f  (BlockLen*sum|w|*4)\n', boundAbs);
    fprintf('  integer bits needed      = %d\n', ceil(log2(boundAbs)));
    fprintf('  + CORDIC gain (x%.4f)  = %.2f -> %d integer bits\n\n', ...
        K_CORDIC, boundAbs*K_CORDIC, ceil(log2(boundAbs*K_CORDIC)));

    %% ---- Empirical populations -------------------------------------
    nBlocks = 20000;
    pad     = halfWindowPad(NTaps);

    % (1) Worst-case aligned QPSK at the corners (all s^4 == -4).
    sCorner = (1 + 1i);                          % |Re|=|Im|=1
    maxA1 = empiricalMax(@() repmat(sCorner, BlockLen + 2*pad, 1), ...
                         w_worst, NTaps, BlockLen, nBlocks);

    % (2) Uniform box Re,Im ~ U(-1,1).
    maxA2 = empiricalMax(@() (2*rand(BlockLen + 2*pad,1)-1) + ...
                              1i*(2*rand(BlockLen + 2*pad,1)-1), ...
                         w_worst, NTaps, BlockLen, nBlocks);

    % (3) Realistic channel, scaled into [-1,1].
    rxScaled = buildRealisticData(Rs, N_pol, DeltaF_Hz, LW_Hz, SNR_dB(wi));
    maxA3 = empiricalMaxFromSignal(rxScaled, w_worst, NTaps, BlockLen);

    fprintf('--- Empirical max |sum4| (real or imag part) ---\n');
    fprintf('  (1) aligned corners      = %8.2f -> %d int bits\n', maxA1, ceil(log2(max(maxA1,1))));
    fprintf('  (2) uniform box U(-1,1)  = %8.2f -> %d int bits\n', maxA2, ceil(log2(max(maxA2,1))));
    fprintf('  (3) realistic, scaled    = %8.2f -> %d int bits\n\n', maxA3, ceil(log2(max(maxA3,1))));

    %% ---- Verdict ----------------------------------------------------
    worstObserved = max([boundAbs*K_CORDIC, maxA1*K_CORDIC, maxA2*K_CORDIC, maxA3*K_CORDIC]);
    fprintf('--- Verdict ---\n');
    fprintf('  available integer range  = +/-%d (%d magnitude bits, IntBits=%d)\n', ...
        AvailRange, AvailIntMagBits, IntBits);
    fprintf('  worst-case need (x gain) = %.1f (%d magnitude bits)\n', ...
        worstObserved, ceil(log2(worstObserved)));
    if worstObserved < AvailRange
        fprintf(['  => sum4 fits with %.1f bits of headroom; default Wrap is\n' ...
                 '     SAFE for IntBits=%d. No saturation/wider accumulator needed.\n'], ...
                 AvailIntMagBits - log2(worstObserved), IntBits);
    else
        fprintf(['  => sum4 can OVERFLOW the available integer range.\n' ...
                 '     Add saturation or widen the T.theta integer part.\n']);
    end
end

% =====================================================================
function pad = halfWindowPad(NTaps)
    pad = NTaps;     % one-sided window half-length
end

function s4 = pow4(s)
    s2 = s .* s;
    s4 = s2 .* s2;
end

function m = empiricalMax(genBlock, w, NTaps, BlockLen, nBlocks)
    L_filt = 2*NTaps + 1;
    pad    = NTaps;
    m = 0;
    for b = 1:nBlocks
        s  = genBlock();                 % length BlockLen + 2*pad
        s4 = pow4(s);
        acc = 0;
        for i = 1:BlockLen
            % window centred on sample (pad + i), taps k = 1..L_filt
            idx = (pad + i) - NTaps + (0:L_filt-1);
            acc = acc + sum(w(:).' .* s4(idx).');
        end
        m = max(m, max(abs(real(acc)), abs(imag(acc))));
    end
end

function m = empiricalMaxFromSignal(rx, w, NTaps, BlockLen)
    L_filt = 2*NTaps + 1;
    m = 0;
    for p = 1:size(rx,2)
        s  = rx(:, p);
        s4 = pow4(s);
        N  = numel(s);
        NB = floor(N / BlockLen);
        for blk = 1:NB
            acc = 0;
            for i = 1:BlockLen
                centre = (blk-1)*BlockLen + i;
                idx = centre - NTaps + (0:L_filt-1);
                idx = idx(idx >= 1 & idx <= N);
                acc = acc + sum(w(1:numel(idx)).' .* s4(idx).');
            end
            m = max(m, max(abs(real(acc)), abs(imag(acc))));
        end
    end
end

function rx = buildRealisticData(Rs, N_pol, DeltaF_Hz, LW_Hz, SNR_dB)
    BITS_PER_SF = 3586 * 2 * 2;
    [symbols, ~, ~, ~] = modem.modulate(modem.randomBits(BITS_PER_SF));
    rx = channel.lo_freq_shift(symbols, DeltaF_Hz / 1e6, Rs, 1);
    rx = channel.add_awgn(rx, SNR_dB);
    rx = channel.add_phase_noise(rx, Rs, LW_Hz);
    % Scale so the peak |Re|/|Im| sample sits at the box edge (=> data in [-1,1]).
    peak = max(max(abs(real(rx(:)))), max(abs(imag(rx(:)))));
    rx   = rx / peak;
end
