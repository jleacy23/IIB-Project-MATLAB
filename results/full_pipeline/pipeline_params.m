function P = pipeline_params()
%PIPELINE_PARAMS  Central parameter definition for the full receiver pipeline.
%
%   P = pipeline_params()
%
%   Returns a struct containing every parameter needed by run_pipeline
%   and build_all_mex.  Edit values here; everything else reads from P.

    % ---- Modulation -------------------------------------------------
    P.M       = 16;           % QAM order (4, 16, 64, …)
    P.N_pol   = 2;            % number of polarisations
    P.Ns      = 2^14;         % symbols per polarisation

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

    % ---- CD Equalizer -----------------------------------------------
    P.NFFT    = 512;          % FFT block size (power of 2)

    % ---- Adaptive Equalizer -----------------------------------------
    P.AEQ_Eq         = 'CMA+RDE';
    P.AEQ_NTaps      = 15;
    P.AEQ_Mu         = 1e-3;
    P.AEQ_SingleSpike = true;
    P.AEQ_N1         = 2000;       % y-pol re-init iteration
    P.AEQ_N2         = 4000;       % CMA→RDE switch iteration
    P.AEQ_NOut       = 500;        % transient discard

    % ---- Viterbi–Viterbi Carrier Recovery ----------------------------
    P.VV_NTaps       = 15;         % half-width of VV averaging filter
    P.VV_P           = 8;          % pilots per block
    P.VV_L           = 256;        % block length [symbols]
    P.VV_CSThreshold = pi/3;       % cycle-slip threshold [rad]
    P.VV_UsePilots   = false;

    % ---- Fixed-point configuration ----------------------------------
    P.FxpConfig      = 'fixed32';  % 'fixed16' | 'fixed32'
    P.po2Twiddle     = false;      % power-of-2 twiddle factors in FFT

    % ---- Output / display --------------------------------------------
    P.Plot    = true;          % set false to suppress constellation plots

    % ---- Random seed (reproducibility) ------------------------------
    P.Seed    = 42;
end
