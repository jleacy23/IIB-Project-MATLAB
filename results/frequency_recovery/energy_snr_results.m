function [tbl_shot, tbl_amp] = energy_snr_results()
% ENERGY_SNR_RESULTS  Optimal DSP algorithm system energy vs FEC SNR.
%
%   For each FR+CR algorithm pair and PON splitting ratio K, finds the
%   fixed-point precision combo that minimises mean system energy
%   (transmitter + K × receiver), then reports:
%
%     mean / std system energy [fJ/bit]
%     mean / std FEC SNR threshold [dB]
%     mean received optical power [dBm]
%
%   Two tables are returned and printed:
%     tbl_shot — 10 km passive-split PON  (shot-noise limited)
%     tbl_amp  — 3 × 80 km EDFA-amplified link (ASE + NLI limited)
%
%   K ∈ {1, 10, 100} in both cases.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src'));

%% System parameters
Rs     = 30.5e9;        % symbol rate [Hz]
M      = 4;             % QPSK
lambda = 1550e-9;       % wavelength [m]
eta    = 0.1;           % laser wall-plug efficiency
K_vec  = [1, 10, 100];  % PON splitting / WDM fan-out ratios

%% Link parameters -------------------------------------------------------

% Shot-noise limited: 10 km unamplified passive-split PON
sn.B        = Rs;       % signal bandwidth = symbol rate [Hz]
sn.alpha_dB = 0.2;      % fibre loss [dB/km]
sn.L        = 10;       % link length [km]

% Amplified: 3 × 80 km EDFA-amplified spans (gain = span loss)
an.n_spans  = 3;
an.L_span   = 80;                       % span length [km]
an.alpha    = 0.2 * log(10) / 10;       % fibre loss [nepers/km]
an.beta2    = 20e-24;                   % |GVD| [s²/km]
an.gamma    = 1.3;                      % NLI coefficient [W⁻¹ km⁻¹]
an.NF_dB    = 5;                        % EDFA noise figure [dB]

%% Load grid-sweep results -----------------------------------------------
d       = load(fullfile(fileparts(mfilename('fullpath')), ...
               'bit_width_full_grid_sweep.mat'));
sweep   = d.tbl;
E_rx_fJ = sweep.energy_fr_fJ + sweep.energy_cr_fJ;  % total Rx DSP energy [fJ/bit]

%% Build tables ----------------------------------------------------------
tbl_shot = build_table(sweep, E_rx_fJ, K_vec, 'shot', Rs, M, lambda, eta, sn, an);
tbl_amp  = build_table(sweep, E_rx_fJ, K_vec, 'amp',  Rs, M, lambda, eta, sn, an);

fprintf('\n=== Shot-noise limited (10 km passive-split PON) ===\n');
disp(tbl_shot);
fprintf('\n=== Amplified-noise limited (3 × 80 km EDFA spans) ===\n');
disp(tbl_amp);

end


function out = build_table(sweep, E_rx_fJ, K_vec, regime, Rs, M, lambda, eta, sn, an)

    combos = unique(sweep(:, {'fr_algo', 'cr_algo'}), 'rows');
    NC = height(combos);
    NK = numel(K_vec);
    NR = NC * NK;

    v_fr  = strings(NR, 1);  v_cr  = strings(NR, 1);
    v_K   = nan(NR, 1);
    v_mE  = nan(NR, 1);  v_sE  = nan(NR, 1);
    v_mSN = nan(NR, 1);  v_sSN = nan(NR, 1);
    v_Prx = nan(NR, 1);

    row = 0;
    for ci = 1:NC
        fr   = combos.fr_algo(ci);
        cr   = combos.cr_algo(ci);
        mask = sweep.fr_algo == fr & sweep.cr_algo == cr;
        idx  = find(mask);
        Erx  = E_rx_fJ(mask);

        for ki = 1:NK
            K = K_vec(ki);

            % Mean system energy for each fixed-point precision combo
            n_configs  = numel(idx);
            mean_E_sys = inf(n_configs, 1);
            for ri = 1:n_configs
                [E_sys, ~, ~] = trial_energies( ...
                    sweep.fec_snr_trials{idx(ri)}, Erx(ri), K, regime, ...
                    Rs, M, lambda, eta, sn, an);
                if ~isempty(E_sys)
                    mean_E_sys(ri) = mean(E_sys);
                end
            end

            % Optimal config minimising mean system energy
            [~, best] = min(mean_E_sys);
            [E_sys, snr_v, P_rx_W] = trial_energies( ...
                sweep.fec_snr_trials{idx(best)}, Erx(best), K, regime, ...
                Rs, M, lambda, eta, sn, an);

            row = row + 1;
            v_fr(row)  = fr;    v_cr(row)  = cr;    v_K(row)  = K;
            v_mE(row)  = mean(E_sys);
            v_sE(row)  = std(E_sys);
            v_mSN(row) = mean(snr_v);
            v_sSN(row) = std(snr_v);
            v_Prx(row) = mean(10 * log10(P_rx_W * 1e3));   % mean dBm
        end
    end

    out = table(v_fr, v_cr, v_K, v_mE, v_sE, v_mSN, v_sSN, v_Prx, ...
        'VariableNames', {'fr_algo', 'cr_algo', 'K', ...
        'mean_sys_energy_fJ', 'std_sys_energy_fJ', ...
        'mean_fec_snr_dB', 'std_fec_snr_dB', 'prx_dBm'});
end


function [E_sys, snr_dB_valid, P_rx_W] = trial_energies( ...
        snr_trials, E_rx_fJ_row, K, regime, Rs, M, lambda, eta, sn, an)
% TRIAL_ENERGIES  Per-trial system energy [fJ/bit], FEC SNR [dB], and
%   received optical power [W] for one sweep row.
%
%   NaN trials and trials where the transmitter model has no solution
%   (SNR exceeds the link maximum) are silently dropped.

    good  = ~isnan(snr_trials);
    snr_dB  = snr_trials(good);
    SNR_lin = 10 .^ (snr_dB / 10);
    N = numel(SNR_lin);

    E_tx_J = nan(N, 1);
    P_rx_W = nan(N, 1);

    for i = 1:N
        try
            switch regime
                case 'shot'
                    e         = energy.transmitter_shot( ...
                        SNR_lin(i), sn.B, lambda, K, sn.alpha_dB, sn.L, M, eta);
                    P_tx      = e * sn.B * log2(M) * eta;
                    P_rx_W(i) = P_tx / K * 10^(-sn.alpha_dB * sn.L / 10);

                case 'amp'
                    B_wdm     = Rs * K;
                    e         = energy.transmitter_amplified( ...
                        SNR_lin(i), B_wdm, lambda, an.n_spans, an.L_span, ...
                        an.alpha, an.beta2, an.gamma, an.NF_dB, M, eta);
                    P_tx      = e * B_wdm * log2(M) * eta;
                    P_rx_W(i) = P_tx / K;   % per-channel, gain-compensated
            end
            E_tx_J(i) = e / K;
        catch
            % Transmitter model has no solution at this SNR — skip trial
        end
    end

    valid        = ~isnan(E_tx_J);
    E_sys        = E_tx_J(valid) * 1e15 + E_rx_fJ_row;
    snr_dB_valid = snr_dB(valid);
    P_rx_W       = P_rx_W(valid);
end
