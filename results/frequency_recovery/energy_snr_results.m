function tbl = energy_snr_results(ignore_fr_energy)
% ENERGY_SNR_RESULTS  Optimal DSP algorithm system energy vs FEC SNR.
%
%   For each FR+CR algorithm pair and PON splitting ratio K, finds the
%   fixed-point precision combo that minimises mean system energy
%   (transmitter + receiver), then reports:
%
%     mean / std system energy [fJ/bit]
%     mean / std FEC SNR threshold [dB]
%     mean received optical power [dBm]
%
%   The network model assumes a flat 35 dB OLT-to-ONU loss and a 35 dB
%   SNR at the OLT (initial NSR floor), shot-noise limited.  Two
%   configurations are evaluated, differing only in the splitting ratio:
%
%     K = 16  — 80 km, 1:16  split
%     K = 512 — 20 km, 1:512 split
%
%   A single table is returned and printed, and a bar graph with two
%   bars (one per K) per algorithm combo is produced.
%
%   ignore_fr_energy (default false) — when true, drops the frequency
%       recovery energy from the receiver total, modelling the limit where
%       FR cost is amortised over infinite symbols.

if nargin < 1 || isempty(ignore_fr_energy)
    ignore_fr_energy = false;
end

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src'));

%% System parameters
Rs     = 30.5e9;        % symbol rate [Hz]
M      = 4;             % QPSK
lambda = 1550e-9;       % wavelength [m]
eta    = 0.1;           % laser wall-plug efficiency
K_vec  = [16, 512];     % PON splitting ratios (80 km / 20 km configs)

%% Network model ---------------------------------------------------------
% Flat OLT-to-ONU loss budget with a fixed SNR floor at the OLT.
net.B       = Rs;       % signal bandwidth = symbol rate [Hz]
net.loss_dB = 35;       % flat OLT-to-ONU power loss Gamma [dB]
net.snr0_dB = 35;       % SNR at the OLT [dB] (initial NSR floor)

%% Load grid-sweep results -----------------------------------------------
d       = load(fullfile(fileparts(mfilename('fullpath')), ...
               'bit_width_full_grid_sweep.mat'));
sweep   = d.tbl;
if ignore_fr_energy
    E_rx_fJ = sweep.energy_cr_fJ;                      % FR cost amortised away
else
    E_rx_fJ = sweep.energy_fr_fJ + sweep.energy_cr_fJ; % total Rx DSP energy [fJ/bit]
end

%% Build table -----------------------------------------------------------
tbl = build_table(sweep, E_rx_fJ, K_vec, M, lambda, eta, net);

fprintf('\n=== Shot-noise limited (35 dB flat loss, 35 dB OLT SNR) ===\n');
disp(tbl);

print_min_fec_snr(sweep);

plot_energy_split(tbl);

end


function out = build_table(sweep, E_rx_fJ, K_vec, M, lambda, eta, net)

    combos = unique(sweep(:, {'fr_algo', 'cr_algo'}), 'rows');
    NC = height(combos);
    NK = numel(K_vec);
    NR = NC * NK;

    v_fr  = strings(NR, 1);  v_cr  = strings(NR, 1);
    v_K   = nan(NR, 1);
    v_flFR = nan(NR, 1);  v_flCR = nan(NR, 1);
    v_blindD = nan(NR, 1);
    v_mE  = nan(NR, 1);  v_sE  = nan(NR, 1);
    v_mErx = nan(NR, 1);  v_mEtx = nan(NR, 1);
    v_mEfr = nan(NR, 1);  v_mEcr = nan(NR, 1);
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
                    sweep.fec_snr_trials{idx(ri)}, Erx(ri), K, ...
                    M, lambda, eta, net);
                if ~isempty(E_sys)
                    mean_E_sys(ri) = mean(E_sys);
                end
            end

            % Optimal config minimising mean system energy
            [~, best] = min(mean_E_sys);
            [E_sys, snr_v, P_rx_W] = trial_energies( ...
                sweep.fec_snr_trials{idx(best)}, Erx(best), K, ...
                M, lambda, eta, net);

            row = row + 1;
            v_fr(row)  = fr;    v_cr(row)  = cr;    v_K(row)  = K;
            v_flFR(row) = sweep.fl_fr(idx(best));
            v_flCR(row) = sweep.fl_cr(idx(best));
            v_blindD(row) = sweep.blind_d(idx(best));
            v_mE(row)  = mean(E_sys);
            v_sE(row)  = std(E_sys);
            v_mErx(row) = Erx(best);
            v_mEtx(row) = mean(E_sys) - Erx(best);
            v_mEfr(row) = sweep.energy_fr_fJ(idx(best));
            v_mEcr(row) = sweep.energy_cr_fJ(idx(best));
            v_mSN(row) = mean(snr_v);
            v_sSN(row) = std(snr_v);
            v_Prx(row) = mean(10 * log10(P_rx_W * 1e3));   % mean dBm
        end
    end

    out = table(v_fr, v_cr, v_K, v_flFR, v_flCR, v_blindD, v_mE, v_sE, v_mEtx, v_mEfr, v_mEcr, v_mErx, v_mSN, v_sSN, v_Prx, ...
        'VariableNames', {'fr_algo', 'cr_algo', 'K', 'fl_fr', 'fl_cr', 'blind_d', ...
        'mean_sys_energy_fJ', 'std_sys_energy_fJ', ...
        'mean_tx_energy_fJ', 'mean_fr_energy_fJ', 'mean_cr_energy_fJ', 'mean_rx_energy_fJ', ...
        'mean_fec_snr_dB', 'std_fec_snr_dB', 'prx_dBm'});
end


function [E_sys, snr_dB_valid, P_rx_W] = trial_energies( ...
        snr_trials, E_rx_fJ_row, K, M, lambda, eta, net)
% TRIAL_ENERGIES  Per-trial system energy [fJ/bit], FEC SNR [dB], and
%   received optical power [W] for one sweep row.
%
%   NaN trials and trials where the transmitter model has no solution
%   (required NSR below the OLT floor) are silently dropped.

    good  = ~isnan(snr_trials);
    snr_dB  = snr_trials(good);
    SNR_lin = 10 .^ (snr_dB / 10);
    N = numel(SNR_lin);

    E_tx_J = nan(N, 1);
    P_rx_W = nan(N, 1);

    nsr0 = 10^(-net.snr0_dB / 10);

    for i = 1:N
        try
            e         = energy.transmitter_shot( ...
                SNR_lin(i), net.B, lambda, net.loss_dB, nsr0, K, M, eta);
            % transmitter_shot already returns per-ONU energy.
            P_tx      = e * K * net.B * log2(M) * eta;
            P_rx_W(i) = P_tx / (K * 10^(net.loss_dB / 10));
            E_tx_J(i) = e;
        catch
            % Transmitter model has no solution at this SNR — skip trial
        end
    end

    valid        = ~isnan(E_tx_J);
    E_sys        = E_tx_J(valid) * 1e15 + E_rx_fJ_row;
    snr_dB_valid = snr_dB(valid);
    P_rx_W       = P_rx_W(valid);
end


function plot_energy_split(tbl)
    figure('Name', 'Shot-noise limited (35 dB flat loss)', ...
        'Color', 'w', 'Position', [100, 100, 900, 500]);
    ax = axes;

    combos = unique(tbl(:, {'fr_algo', 'cr_algo'}), 'rows', 'stable');
    NC  = height(combos);
    K_u = unique(tbl.K);
    NK  = numel(K_u);
    gap = 1.5;   % extra units of space between algorithm groups

    % Assign x positions: bars within a group are contiguous, groups are separated
    xpos = zeros(1, NC * NK);
    for ci = 1:NC
        xpos((ci-1)*NK + (1:NK)) = (ci-1)*(NK + gap) + (1:NK);
    end

    data = [tbl.mean_tx_energy_fJ, tbl.mean_fr_energy_fJ, tbl.mean_cr_energy_fJ];
    bar(ax, xpos, data, 'stacked');

    % One tick per algorithm combo, centred on its group of K bars
    tickPos = arrayfun(@(ci) (ci-1)*(NK + gap) + (NK + 1)/2, 1:NC);
    tickLbl = strings(NC, 1);
    for ci = 1:NC
        tickLbl(ci) = sprintf('%s / %s', ...
            abbrevAlgo(combos.fr_algo(ci)), abbrevAlgo(combos.cr_algo(ci)));
    end

    ax.XTick = tickPos;
    ax.XTickLabel = tickLbl;
    ax.XTickLabelRotation = 20;
    ax.FontSize = 10;
    grid(ax, 'on');
    ylabel(ax, 'Energy [fJ/bit]', 'FontSize', 11);
    title(ax, sprintf('System energy split (left to right: K = %d, %d)', ...
        K_u(1), K_u(2)), 'FontSize', 11, 'Interpreter', 'none');
    legend(ax, {'TX', 'FR (freq. recovery)', 'CR (phase recovery)'}, ...
        'Location', 'best', 'FontSize', 9, 'Interpreter', 'none');
end


function print_min_fec_snr(sweep)
% PRINT_MIN_FEC_SNR  Lowest mean FEC SNR threshold achievable by each
%   FR/CR algorithm combo, taken over the fixed-point precision grid
%   (i.e. ignoring energy cost), with the precision combo and trial
%   spread at which the minimum occurs.

    combos = unique(sweep(:, {'fr_algo', 'cr_algo'}), 'rows');
    NC = height(combos);

    v_fr  = strings(NC, 1);  v_cr  = strings(NC, 1);
    v_min = nan(NC, 1);      v_std = nan(NC, 1);
    v_flFR = nan(NC, 1);     v_flCR = nan(NC, 1);  v_bd = nan(NC, 1);

    for ci = 1:NC
        fr   = combos.fr_algo(ci);
        cr   = combos.cr_algo(ci);
        idx  = find(sweep.fr_algo == fr & sweep.cr_algo == cr);

        best_min = inf;  best_std = NaN;  best_r = NaN;
        for r = idx.'
            t = sweep.fec_snr_trials{r};
            t = t(~isnan(t));
            if isempty(t)
                continue
            end
            m = mean(t);
            if m < best_min
                best_min = m;
                best_std = std(t);
                best_r   = r;
            end
        end

        v_fr(ci) = fr;  v_cr(ci) = cr;
        if isfinite(best_min)
            v_min(ci)  = best_min;
            v_std(ci)  = best_std;
            v_flFR(ci) = sweep.fl_fr(best_r);
            v_flCR(ci) = sweep.fl_cr(best_r);
            v_bd(ci)   = sweep.blind_d(best_r);
        end
    end

    out = table(v_fr, v_cr, v_min, v_std, v_flFR, v_flCR, v_bd, ...
        'VariableNames', {'fr_algo', 'cr_algo', 'min_fec_snr_dB', ...
        'std_fec_snr_dB', 'fl_fr', 'fl_cr', 'blind_d'});

    fprintf('\n=== Minimum FEC SNR achieved by each algorithm combo ===\n');
    disp(out);
end


function s = abbrevAlgo(name)
    map = {'fft_search',       'R&B';       ...
           'fft_search_blind', 'R&B blind'; ...
           'differential_kay', 'DK';        ...
           'viterbi_viterbi',  'V&V';       ...
           'pilots_only',      'PO'};
    idx = strcmp(map(:, 1), name);
    if any(idx)
        s = map{idx, 2};
    else
        s = name;
    end
end
