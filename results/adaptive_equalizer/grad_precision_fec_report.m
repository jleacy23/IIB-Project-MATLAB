function tbl = grad_precision_fec_report(matFile)
%GRAD_PRECISION_FEC_REPORT  Tabulate mean FEC SNR and transmitter energy.
%
%   tbl = grad_precision_fec_report()
%   tbl = grad_precision_fec_report(matFile)
%
%   Loads the sweep produced by grad_precision_fec (default
%   grad_precision_fec_sweep.mat next to this file) and prints, for
%   every implementation (gradient fraction length x weight-update
%   step x SignOnly x Mu-scaling x NTaps) and each network config
%   (fibre length / splitting ratio), the FEC SNR averaged over the
%   Monte-Carlo runs together with the resulting wall-plug transmitter
%   energy per bit.
%
%   The transmitter energy is obtained exactly as in
%   results/frequency_recovery/energy_snr_results.m: each per-trial FEC
%   SNR is pushed through energy.transmitter_shot under the shot-noise
%   network model (35 dB flat OLT-to-ONU loss, 35 dB OLT SNR floor),
%   and the mean and standard deviation are taken directly over the
%   per-trial energies.  The two network configs are
%
%     L = 80 km , K = 16   (1:16  split)
%     L = 20 km , K = 512  (1:512 split)
%
%   Columns
%     GradFL        - gradient-estimate fraction length (bits)
%     UpdateStep    - samples between weight updates
%     SignOnly      - sign-reduced CMA update (false/true)
%     MuScaling     - step size: fixed, or scaled by UpdateStep
%     NTaps         - number of equaliser FIR taps
%     L_km          - fibre length [km]
%     K             - passive splitter ratio
%     MeanFEC_SNR   - mean FEC-limit SNR over trials [dB] (NaN-omitted)
%     StdFEC_SNR    - std of the FEC-limit SNR over trials [dB]
%     MeanE_tx_fJ   - mean transmitter energy per bit [fJ/bit]
%     StdE_tx_fJ    - std of the transmitter energy per bit [fJ/bit]
%     NValid        - number of trials that crossed the FEC limit
%
%   The returned table is sorted by MeanE_tx_fJ ascending (most
%   energy-efficient implementation first); NaN rows (no trial ever
%   crossed the FEC limit) sort last.

    if nargin < 1 || isempty(matFile)
        matFile = fullfile(fileparts(mfilename('fullpath')), ...
            'grad_precision_fec_sweep.mat');
    end
    if ~isfile(matFile)
        error('grad_precision_fec_report:NoFile', ...
            'Sweep file not found: %s\nRun grad_precision_fec first.', ...
            matFile);
    end

    addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src'));

    S = load(matFile, 'fecSNR', 'params');
    fecSNR = S.fecSNR;  % [NFL x NUS x NSO x NMS x NNT x NL x NTrials]
    P      = S.params;

    GradFL_vec     = P.GradFL_vec;
    UpdateStep_vec = P.UpdateStep_vec;
    SignOnly_dim   = P.SignOnly_dim;    % {'false','true'}
    MuScaling_dim  = P.MuScaling_dim;   % {'fixed','scaled_by_UpdateStep'}
    NTaps_vec      = P.NTaps_vec;
    L_km_vec       = P.L_km_vec;        % [80 20]
    K_vec          = P.K_vec;           % [16 512] aligned with L_km_vec

    NFL = numel(GradFL_vec);
    NUS = numel(UpdateStep_vec);
    NSO = numel(SignOnly_dim);
    NMS = numel(MuScaling_dim);
    NNT = numel(NTaps_vec);
    NL  = numel(L_km_vec);

    % Shot-noise network model (matches energy_snr_results.m).
    netc.B       = 30.5e9;   % signal bandwidth = symbol rate [Hz]
    netc.lambda  = 1550e-9;  % wavelength [m]
    netc.M       = 4;        % QPSK
    netc.eta     = 0.1;      % laser wall-plug efficiency
    netc.loss_dB = 35;       % flat OLT-to-ONU power loss [dB]
    netc.nsr0    = 10^(-35/10);  % initial NSR at the OLT (35 dB)

    nRows = NFL * NUS * NSO * NMS * NNT * NL;
    GradFL      = zeros(nRows, 1);
    UpdateStep  = zeros(nRows, 1);
    SignOnly    = strings(nRows, 1);
    MuScaling   = strings(nRows, 1);
    NTaps       = zeros(nRows, 1);
    L_km        = zeros(nRows, 1);
    K           = zeros(nRows, 1);
    MeanFEC_SNR = zeros(nRows, 1);
    StdFEC_SNR  = zeros(nRows, 1);
    MeanE_tx_fJ = zeros(nRows, 1);
    StdE_tx_fJ  = zeros(nRows, 1);
    NValid      = zeros(nRows, 1);

    r = 0;
    for fl = 1:NFL
        for us = 1:NUS
            for so = 1:NSO
                for ms = 1:NMS
                    for nt = 1:NNT
                        for nl = 1:NL
                            r = r + 1;
                            trials = squeeze( ...
                                fecSNR(fl, us, so, ms, nt, nl, :));
                            Kv = K_vec(nl);

                            [mE, sE, nV] = energyStats(trials, Kv, netc);

                            GradFL(r)      = GradFL_vec(fl);
                            UpdateStep(r)  = UpdateStep_vec(us);
                            SignOnly(r)    = string(SignOnly_dim{so});
                            MuScaling(r)   = string(MuScaling_dim{ms});
                            NTaps(r)       = NTaps_vec(nt);
                            L_km(r)        = L_km_vec(nl);
                            K(r)           = Kv;
                            MeanFEC_SNR(r) = mean(trials, 'omitnan');
                            StdFEC_SNR(r)  = std(trials,  'omitnan');
                            MeanE_tx_fJ(r) = mE;
                            StdE_tx_fJ(r)  = sE;
                            NValid(r)      = nV;
                        end
                    end
                end
            end
        end
    end

    tbl = table(GradFL, UpdateStep, SignOnly, MuScaling, NTaps, ...
        L_km, K, MeanFEC_SNR, StdFEC_SNR, MeanE_tx_fJ, StdE_tx_fJ, ...
        NValid);

    % Most energy-efficient first; never-crossed rows (NaN) last.
    [~, order] = sortrows([MeanE_tx_fJ, isnan(MeanE_tx_fJ)], [2 1]);
    tbl = tbl(order, :);

    fprintf(['\nFEC SNR & transmitter energy per implementation ', ...
        '(FEC BER = %.0e, %d trials)\n'], P.FEC_BER, size(fecSNR, 7));
    disp(tbl);
end


function [meanE_fJ, stdE_fJ, nValid] = energyStats(trials_dB, K, netc)
% ENERGYSTATS  Mean / std wall-plug transmitter energy [fJ/bit] from the
%   per-trial FEC SNRs, computed exactly as in energy_snr_results.m:
%   each finite trial SNR is pushed through energy.transmitter_shot and
%   the statistics are taken directly over the per-trial energies.

    snr_dB  = trials_dB(~isnan(trials_dB));
    SNR_lin = 10 .^ (snr_dB / 10);
    N = numel(SNR_lin);

    E_fJ = nan(N, 1);
    for i = 1:N
        try
            e = energy.transmitter_shot(SNR_lin(i), netc.B, ...
                netc.lambda, netc.loss_dB, netc.nsr0, K, ...
                netc.M, netc.eta);
            E_fJ(i) = e * 1e15;
        catch
            % Transmitter model has no solution at this SNR — skip trial
        end
    end

    valid    = ~isnan(E_fJ);
    nValid   = sum(valid);
    if nValid == 0
        meanE_fJ = NaN;
        stdE_fJ  = NaN;
    else
        meanE_fJ = mean(E_fJ(valid));
        stdE_fJ  = std(E_fJ(valid));
    end
end
