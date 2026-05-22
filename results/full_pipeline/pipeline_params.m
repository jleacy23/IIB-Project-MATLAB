function P = pipeline_params()
%PIPELINE_PARAMS  Central parameter definition for the full receiver pipeline.
%
%   P = pipeline_params()
%
%   Returns a struct containing every parameter needed by run_pipeline
%   and build_all_mex.  Edit values here; everything else reads from P.

    % ---- Modulation -------------------------------------------------
    P.M       = 4;           % QPSK only (kept for carrier-recovery APIs)
    P.N_pol   = 2;            % number of polarisations
    P.Ns      = 2^17;         % approximate symbols per polarisation

    % ---- Pulse shaping ----------------------------------------------
    P.SpS     = 2;            % samples per symbol
    P.Rolloff = 0.25;         % RRC roll-off factor
    P.Span    = 10;           % RRC filter span [symbols]

    % ---- System / channel -------------------------------------------
    P.Rs      = 32;           % symbol rate  [GBd]
    P.L       = 80;           % fibre length [km]
    P.D       = 17;           % dispersion   [ps/(nm·km)]
    P.CWL     = 1550;         % central wavelength [nm]
    P.SNR_dB  = 25;           % AWGN SNR [dB]
    P.LW      = 100e3;        % laser linewidth [Hz]  (phase noise)
    P.DGDSpec = 0.5;          % PMD coefficient [ps/sqrt(km)]
    P.N_pmd   = 10;           % number of PMD sections
    P.ENOBits = 5;            % ADC effective number of bits

    % ---- CD Equalizer -----------------------------------------------
    P.NFFT    = 2^9;          % FFT block size (power of 2)

    % ---- Adaptive Equalizer -----------------------------------------
    P.AEQ_NTaps      = 15;
    P.AEQ_Mu         = 2e-9;
    P.AEQ_SingleSpike = true;
    P.AEQ_N1         = 2000;       % y-pol re-init iteration
    P.AEQ_NOut       = 5000;        % transient discard
    P.AEQ_SignOnly   = false;      % if true, use sign(err) & complex-sign(y) updates
    P.AEQ_PLanes     = 1;          % parallel lanes (1 = serial fixed-point CMA)
    P.AEQ_Mode       = 0;          % 0 = CMA, 1 = pilot-aided LMS
    P.AEQ_BlockLen   = 32;         % CPON weight-update block (32 = one block)

    % ---- Carrier Recovery ----------------------------
    P.VV_NTaps       = 5;         % half-width of VV averaging filter
    P.BPS_B   = 64;         % number of test phases for BPS
    P.BPS_N   = 5;
    P.UsePilots       = true;       % use pilot symbols for phase estimation
    P.BlockBased       = false;      % apply phase correction on a block of symbols

    % ---- Fixed-point configuration ----------------------------------
    P.FxpConfig_CD   = 'fixed16';  % CD equalizer:  'fixed16' | 'fixed32'
    P.FxpConfig_AEQ  = 'fixed16';  % Adaptive EQ:   'fixed16' | 'fixed32'
    P.FxpConfig_VV   = 'fixed16';  % VV carrier recovery: 'fixed16' | 'fixed32'
    P.FxpConfig_BPS  = 'fixed16';  % BPS carrier recovery: 'fixed16' | 'fixed32'
    P.po2Twiddle     = false;      % power-of-2 twiddle factors in FFT

    % ---- Output / display --------------------------------------------
    P.Plot    = true;          % set false to suppress constellation plots

    % ---- Random seed (reproducibility) ------------------------------
    P.Seed    = 42;
end
