function [symbols, pilots, training, nSubframes] = modulate(bits)
%MODULATE  DP-QPSK modulation with CPON subframe, pilot and training insertion.
%
%   [symbols, pilots, training, nSubframes] = modulate(bits)
%
%   Builds complete CPON downstream DSP subframes:
%     - Each subframe is 3712 symbols (116 blocks of 32).
%     - Block 1: TS1..TS11 (training), then 21 data symbols.
%       TS1 also serves as pilot index 1.
%     - Blocks 2-116: 31 data symbols + 1 pilot at position 32.
%     - Pilots are at +/-1 +/-1j (unit amplitude, matching data symbols),
%       generated from PRBS10.
%     - Data symbols are at +/-1 +/-1j (standard QPSK).
%
%   Bit-to-symbol mapping (per the CPON spec, for symbol index i):
%     c(4i)   -> I of X-pol
%     c(4i+2) -> Q of X-pol
%     c(4i+1) -> I of Y-pol
%     c(4i+3) -> Q of Y-pol
%   where amplitude = 2*bit - 1  (0 -> -1, 1 -> +1).
%
%   Input
%     bits  - column vector of data bits
%
%   Outputs
%     symbols    - [Nsym x 2] complex QPSK symbols (X-pol, Y-pol) with
%                  pilots and training already inserted
%     pilots     - [116 x 2] pilot symbols for one subframe (same every
%                  subframe; amplitude 1)
%     training   - [11 x 2] training symbols (amplitude 1)
%     nSubframes - number of complete subframes generated

    %% ------------------------------------------------------------------
    %  Constants
    % -------------------------------------------------------------------
    N_pol          = 2;
    SUBFRAME_SYMS  = 3712;           % symbols per full subframe
    BLOCK_LEN      = 32;             % symbols per block
    N_BLOCKS       = 116;            % blocks per subframe
    N_TRAIN        = 11;             % training symbols at start of block 1
    DATA_PER_SUBFRAME = 3586;        % 3712 - 116 pilots - 10 non-pilot TS

    %% ------------------------------------------------------------------
    %  Training symbols (fixed, from spec Table 6.1)
    % -------------------------------------------------------------------
    training = [ ...
        -1+1j, -1-1j; ...   % TS1  (also pilot index 1)
        +1+1j, -1-1j; ...   % TS2
        -1+1j, +1-1j; ...   % TS3
        +1+1j, -1+1j; ...   % TS4
        -1-1j, -1+1j; ...   % TS5
        +1+1j, +1+1j; ...   % TS6
        -1-1j, -1-1j; ...   % TS7
        -1-1j, -1+1j; ...   % TS8
        +1+1j, +1-1j; ...   % TS9
        +1-1j, +1+1j; ...   % TS10
        +1-1j, +1-1j; ...   % TS11
    ];   % [11 x 2]

    %% ------------------------------------------------------------------
    %  Pilot symbols via PRBS10
    % -------------------------------------------------------------------
    pilots = generatePilots(N_BLOCKS);   % [116 x 2]

    %% ------------------------------------------------------------------
    %  Compute number of subframes & pad data bits
    % -------------------------------------------------------------------
    bits = bits(:);
    bitsPerSymbol  = 2;                      % QPSK: 2 bits per symbol per pol
    bitsPerSubframe = DATA_PER_SUBFRAME * N_pol * bitsPerSymbol;

    nSubframes = ceil(length(bits) / bitsPerSubframe);
    nBitsNeeded = nSubframes * bitsPerSubframe;

    if length(bits) < nBitsNeeded
        bits = [bits; randi([0 1], nBitsNeeded - length(bits), 1)];
    end

    %% ------------------------------------------------------------------
    %  Map data bits -> QPSK data symbols  [NdataTotal x 2]
    %
    %  Bit interleaving per spec:
    %    c(4i)   -> X_I,  c(4i+1) -> Y_I,  c(4i+2) -> X_Q,  c(4i+3) -> Y_Q
    % -------------------------------------------------------------------
    nDataTotal = nSubframes * DATA_PER_SUBFRAME;
    bitsUsed   = bits(1 : nDataTotal * N_pol * bitsPerSymbol);

    % Reshape to [nDataTotal, 4]:  columns = c(4i), c(4i+1), c(4i+2), c(4i+3)
    bMat = reshape(bitsUsed, 4, []).';   % [nDataTotal x 4]

    XI = 2*bMat(:,1) - 1;   % c(4i)   -> I of X
    YI = 2*bMat(:,2) - 1;   % c(4i+1) -> I of Y
    XQ = 2*bMat(:,3) - 1;   % c(4i+2) -> Q of X
    YQ = 2*bMat(:,4) - 1;   % c(4i+3) -> Q of Y

    dataX = complex(XI, XQ);   % [nDataTotal x 1]
    dataY = complex(YI, YQ);

    %% ------------------------------------------------------------------
    %  Assemble subframes
    % -------------------------------------------------------------------
    totalSyms = nSubframes * SUBFRAME_SYMS;
    symbols   = zeros(totalSyms, N_pol);

    dataIdx = 0;   % running index into data symbol vectors

    for sf = 1:nSubframes
        base = (sf - 1) * SUBFRAME_SYMS;   % 0-based offset in output

        % --- Block 1: training (11) + data (21) ---
        % Positions 1..11  = training
        symbols(base + (1:N_TRAIN), :) = training;

        % Positions 12..32 = data
        nData1 = BLOCK_LEN - N_TRAIN;   % 21
        symbols(base + N_TRAIN + (1:nData1), 1) = dataX(dataIdx + (1:nData1));
        symbols(base + N_TRAIN + (1:nData1), 2) = dataY(dataIdx + (1:nData1));
        dataIdx = dataIdx + nData1;

        % --- Blocks 2..116: 1 pilot + 31 data ---
        for blk = 2:N_BLOCKS
            blkBase = base + (blk - 1) * BLOCK_LEN;

            % 1 pilot at position 1 of the block (first symbol)
            symbols(blkBase + 1, :) = pilots(blk, :);

            % 31 data symbols at positions 2..32
            nDataBlk = BLOCK_LEN - 1;   % 31
            symbols(blkBase + (2:BLOCK_LEN), 1) = dataX(dataIdx + (1:nDataBlk));
            symbols(blkBase + (2:BLOCK_LEN), 2) = dataY(dataIdx + (1:nDataBlk));
            dataIdx = dataIdx + nDataBlk;
        end

        % Overwrite position 1 with pilot (TS1 = pilot index 1, already done
        % via training, but ensure consistency)
        symbols(base + 1, :) = pilots(1, :);
    end
end

%% ======================================================================
%  Local: PRBS10 pilot generator
% =======================================================================
function pilots = generatePilots(nPilots)
%GENERATEPILOTS  Produce nPilots QPSK pilot symbols for X-pol and Y-pol
%   using the CPON PRBS10 generator (polynomial x^10+x^8+x^4+x^3+1).
%
%   Each pilot is at +/-1 +/-1j (unit amplitude, matching data symbols).

    % Seeds (10-bit, MSB-first)
    seedX = bitand(uint16(hex2dec('19E')), uint16(1023));
    seedY = bitand(uint16(hex2dec('0D0')), uint16(1023));

    bitsX = prbs10(seedX, nPilots * 2);   % 2 bits per pilot symbol
    bitsY = prbs10(seedY, nPilots * 2);

    % Map bit pairs to +/-3 +/-3j
    pilots = zeros(nPilots, 2);
    for k = 1:nPilots
        ix = 2*(k-1) + 1;
        Ix = 2*bitsX(ix)   - 1;   Qx = 2*bitsX(ix+1) - 1;
        Iy = 2*bitsY(ix)   - 1;   Qy = 2*bitsY(ix+1) - 1;
        pilots(k, 1) = complex(Ix, Qx);
        pilots(k, 2) = complex(Iy, Qy);
    end
end

function bits = prbs10(seed, nBits)
%PRBS10  Generate nBits from a 10-bit LFSR.
%   Polynomial: x^10 + x^8 + x^4 + x^3 + 1
%   Taps at bits 10, 8, 4, 3 (1-indexed from MSB).
%   Output is the MSB of the register at each clock.

    reg = seed;
    bits = zeros(1, nBits);
    for k = 1:nBits
        bits(k) = double(bitand(bitshift(reg, -9), uint16(1)));   % MSB out
        % Feedback: XOR of bits 10, 8, 4, 3 (0-indexed: 9, 7, 3, 2)
        fb = bitxor(bitxor(bitshift(reg, -9), bitshift(reg, -7)), ...
             bitxor(bitshift(reg, -3), bitshift(reg, -2)));
        fb = bitand(fb, uint16(1));
        reg = bitand(bitor(bitshift(reg, 1), fb), uint16(1023));
    end
end
