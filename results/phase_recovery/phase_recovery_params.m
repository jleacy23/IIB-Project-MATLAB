function P = phase_recovery_params()
%PHASE_RECOVERY_PARAMS  Parameters for the block-based phase recovery sweep.
%
%   P = phase_recovery_params()
%
%   Returns a struct containing every parameter needed by run_phase_recovery.
%   Edit values here; the test function reads exclusively from P.
%
%   The primary independent variable is P.BlockSizes — the array of block
%   lengths swept by the test.  BlockBased is always true in this test;
%   the block size determines how often the phase estimate is refreshed.

    % ---- Modulation -------------------------------------------------
    P.M      = 4;           % QAM order
    P.N_pol  = 2;           % number of polarisations
    P.Ns     = 2^13;        % symbols per polarisation per trial
                            % Large enough for reliable BER at each block size

    % ---- Pilots -----------------------------------------------------
    P.PilotLen = 4;         % pilot symbols at start of each block
    P.UsePilots = true;     % pilot-aided cycle-slip correction

    % ---- System -----------------------------------------------------
    P.Rs     = 10;          % symbol rate [GBd]
    P.SNR_dB = 17.5;        % AWGN SNR [dB]
    P.LW     = 2400e3;      % laser linewidth [Hz]  (phase noise)

    % ---- Viterbi-Viterbi specific -----------------------------------
    P.VV_NTaps = 5;         % one-sided VV filter half-length

    % ---- BPS specific -----------------------------------------------
    P.BPS_N = 5;            % one-sided BPS window half-length
    P.BPS_B = 64;           % number of blind test phases (must be even)

    % ---- Sweep: block sizes to evaluate ----------------------------
    % Powers of 2 from 8 to 512.  The lower bound is constrained by
    % PilotLen (block must be longer than the pilot sequence).
    % The upper bound is where block-based phase tracking degrades
    % significantly for the linewidth and symbol rate chosen.
    P.BlockSizes = [8, 16, 32, 64, 128, 256, 512, 1024];

    % ---- Fixed-point ------------------------------------------------
    P.FxpConfig = 'fixed16';    % 'fixed16' | 'fixed32'

    % ---- Averaging over trials --------------------------------------
    % Number of independent channel realisations per (algorithm, block size)
    % pair.  Averaging reduces variance in the BER estimate.
    P.NTrials = 20;

    % ---- Plot -------------------------------------------------------
    P.Plot = true;
end