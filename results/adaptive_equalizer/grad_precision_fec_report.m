function tbl = grad_precision_fec_report(matFile)
%GRAD_PRECISION_FEC_REPORT  Tabulate mean FEC SNR per implementation.
%
%   tbl = grad_precision_fec_report()
%   tbl = grad_precision_fec_report(matFile)
%
%   Loads the sweep produced by grad_precision_fec (default
%   grad_precision_fec_sweep.mat next to this file) and prints, for
%   every implementation (gradient fraction length x weight-update
%   step x SignOnly x Mu-scaling x NTaps), the FEC SNR averaged over
%   the Monte-Carlo runs.
%
%   Columns
%     GradFL        - gradient-estimate fraction length (bits)
%     UpdateStep    - samples between weight updates
%     SignOnly      - sign-reduced CMA update (false/true)
%     MuScaling     - step size: fixed, or scaled by UpdateStep
%     NTaps         - number of equaliser FIR taps
%     MeanFEC_SNR   - mean FEC-limit SNR over trials [dB] (NaN-omitted)
%     StdFEC_SNR    - std of the FEC-limit SNR over trials [dB]
%     NValid        - number of trials that crossed the FEC limit
%
%   The returned table is sorted by MeanFEC_SNR ascending (best
%   implementation - lowest required SNR - first); NaN rows (no trial
%   ever crossed the FEC limit) sort last.

    if nargin < 1 || isempty(matFile)
        matFile = fullfile(fileparts(mfilename('fullpath')), ...
            'grad_precision_fec_sweep.mat');
    end
    if ~isfile(matFile)
        error('grad_precision_fec_report:NoFile', ...
            'Sweep file not found: %s\nRun grad_precision_fec first.', ...
            matFile);
    end

    S = load(matFile, 'fecSNR', 'params');
    fecSNR = S.fecSNR;          % [NFL x NUS x 2 x 2 x NNTaps x NTrials]
    P      = S.params;

    GradFL_vec     = P.GradFL_vec;
    UpdateStep_vec = P.UpdateStep_vec;
    SignOnly_dim   = P.SignOnly_dim;    % {'false','true'}
    MuScaling_dim  = P.MuScaling_dim;   % {'fixed','scaled_by_UpdateStep'}
    NTaps_vec      = P.NTaps_vec;

    NFL = numel(GradFL_vec);
    NUS = numel(UpdateStep_vec);
    NSO = numel(SignOnly_dim);
    NMS = numel(MuScaling_dim);
    NNT = numel(NTaps_vec);

    nRows = NFL * NUS * NSO * NMS * NNT;
    GradFL      = zeros(nRows, 1);
    UpdateStep  = zeros(nRows, 1);
    SignOnly    = strings(nRows, 1);
    MuScaling   = strings(nRows, 1);
    NTaps       = zeros(nRows, 1);
    MeanFEC_SNR = zeros(nRows, 1);
    StdFEC_SNR  = zeros(nRows, 1);
    NValid      = zeros(nRows, 1);

    r = 0;
    for fl = 1:NFL
        for us = 1:NUS
            for so = 1:NSO
                for ms = 1:NMS
                    for nt = 1:NNT
                        r = r + 1;
                        trials = squeeze(fecSNR(fl, us, so, ms, nt, :));

                        GradFL(r)      = GradFL_vec(fl);
                        UpdateStep(r)  = UpdateStep_vec(us);
                        SignOnly(r)    = string(SignOnly_dim{so});
                        MuScaling(r)   = string(MuScaling_dim{ms});
                        NTaps(r)       = NTaps_vec(nt);
                        MeanFEC_SNR(r) = mean(trials, 'omitnan');
                        StdFEC_SNR(r)  = std(trials,  'omitnan');
                        NValid(r)      = sum(isfinite(trials));
                    end
                end
            end
        end
    end

    tbl = table(GradFL, UpdateStep, SignOnly, MuScaling, NTaps, ...
        MeanFEC_SNR, StdFEC_SNR, NValid);

    % Best (lowest required SNR) first; never-crossed rows (NaN) last.
    [~, order] = sortrows([MeanFEC_SNR, isnan(MeanFEC_SNR)], [2 1]);
    tbl = tbl(order, :);

    fprintf('\nFEC SNR per implementation (FEC BER = %.0e, %d trials)\n', ...
        P.FEC_BER, size(fecSNR, 6));
    disp(tbl);
end
