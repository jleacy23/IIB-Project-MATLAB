function R = run_pipeline(P)
%RUN_PIPELINE  Full receiver pipeline — floating-point and fixed-point paths.
%
%   R = run_pipeline(P)
%
%   Runs:  bits → symbols → RRC pulse shaping → channel impairments
%          → CD equalisation → matched filtering → adaptive equalisation
%          → Viterbi-Viterbi carrier recovery
%
%   Two parallel paths are executed:
%     1. Floating-point  (double)
%     2. Fixed-point via MEX  (uses *_fxp_mex functions)
%
%   Input
%     P  - parameter struct from pipeline_params()
%
%   Output
%     R  - results struct with fields for every stage (see bottom)

    addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src'));
    rng(P.Seed);

    k     = log2(P.M);
    Nbits = k * P.N_pol * P.Ns;

    %% ================================================================
    %  Transmitter (shared by both paths)
    % =================================================================
    fprintf('TX: %d-QAM, %d pol, %d symbols/pol\n', P.M, P.N_pol, P.Ns);

    bits    = qam_randomBits(Nbits);
    symbols = qam_modulate(bits, P.M, P.N_pol);        % [Ns x Npol]
    txSig   = qam_rrcPulse(symbols, P.SpS, P.Rolloff, P.Span);

    %% ================================================================
    %  Channel
    % =================================================================
    fprintf('Channel: SNR=%.0f dB, D=%.0f ps/(nm·km), L=%.0f km, LW=%.0f kHz\n', ...
            P.SNR_dB, P.D, P.L, P.LW/1e3);

    rxSig = channel_add_awgn(txSig, P.SNR_dB);
    rxSig = channel_add_chromatic_dispersion(rxSig, P.L, P.SpS, P.Rs, P.D, P.CWL);
    rxSig = channel_add_phase_noise(rxSig, P.Rs, P.LW);
    rxSig = channel_add_pmd(rxSig, P.L, P.SpS, P.Rs, P.DGDSpec, P.N_pmd);

    %% ================================================================
    %  Generate VV filter & pilot symbols (shared)
    % =================================================================
    SymbolEnergy = 1;   % unit-power constellation
    VVFilter = cr_genVVFilter(P.LW, P.Rs, P.SNR_dB, SymbolEnergy, ...
                              P.N_pol, P.VV_NTaps);

    % Pilot symbols: first P.VV_P symbols of the transmitted sequence
    Pilots = symbols(1:P.VV_P, :);

    %% ================================================================
    %  PATH A — Floating-point (double)
    % =================================================================
    fprintf('\n--- Floating-point path ---\n');

    % CD Equalisation
    fprintf('  CD EQ ... ');
    tic;
    cdOut_fl = cdeq_equalize(rxSig, P.D, P.L, P.CWL, P.Rs, ...
                             P.N_pol, P.SpS, P.NFFT);
    t_cd_fl = toc;
    fprintf('%.3f s\n', t_cd_fl);

    % Matched filter
    mfOut_fl = qam_matched_filter(cdOut_fl, P.SpS, 'rrc', P.Rolloff, P.Span);

    % Adaptive Equalisation
    fprintf('  Adaptive EQ ... ');
    tic;
    aeqOut_fl = adeq_equalize(mfOut_fl, P.SpS, P.AEQ_Eq, ...
                              P.AEQ_NTaps, P.AEQ_Mu, P.AEQ_SingleSpike, ...
                              P.AEQ_N1, P.AEQ_N2, P.AEQ_NOut);
    t_aeq_fl = toc;
    fprintf('%.3f s\n', t_aeq_fl);

    % Viterbi-Viterbi
    fprintf('  VV carrier recovery ... ');
    tic;
    vvOut_fl = cr_viterbiViterbi(aeqOut_fl, P.N_pol, P.VV_NTaps, ...
                                 VVFilter, Pilots, P.VV_P, P.VV_L, ...
                                 P.VV_CSThreshold, P.VV_UsePilots);
    t_vv_fl = toc;
    fprintf('%.3f s\n', t_vv_fl);

    % Decide & compute BER
    dec_fl  = qam_decideSymbols(vvOut_fl, P.M, P.N_pol);
    bits_fl = qam_symbolsToBits(dec_fl, P.M);
    Ncomp   = min(length(bits_fl), length(bits));
    BER_fl  = sum(bits_fl(1:Ncomp) ~= bits(1:Ncomp)) / Ncomp;
    fprintf('  BER (float) = %.2e  (%d / %d bits)\n', BER_fl, ...
            sum(bits_fl(1:Ncomp) ~= bits(1:Ncomp)), Ncomp);

    %% ================================================================
    %  PATH B — Fixed-point (MEX)
    % =================================================================
    fprintf('\n--- Fixed-point MEX path (%s) ---\n', P.FxpConfig);

    % Load types tables
    T_cd  = cdeq_equalize_fxp_types(P.FxpConfig);
    T_aeq = adeq_equalize_fxp_types(P.FxpConfig);
    T_vv  = cr_viterbiViterbi_fxp_types(P.FxpConfig);

    % Cast channel output to fi for the fixed-point path
    rxSig_fi = cast(rxSig, 'like', T_cd.x);

    % CD Equalisation (MEX)
    fprintf('  CD EQ (MEX) ... ');
    tic;
    cdOut_fxp = cdeq_equalize_fxp_mex(rxSig_fi, ...
                    double(P.D), double(P.L), double(P.CWL), ...
                    double(P.Rs), double(P.N_pol), double(P.SpS), ...
                    double(P.NFFT), logical(P.po2Twiddle), T_cd);
    t_cd_fxp = toc;
    fprintf('%.3f s\n', t_cd_fxp);

    % Matched filter (float — not fixed-point)
    mfOut_fxp = qam_matched_filter(double(cdOut_fxp), P.SpS, 'rrc', ...
                                   P.Rolloff, P.Span);

    % Cast back to fi for adaptive EQ
    mfOut_fxp_fi = cast(mfOut_fxp, 'like', T_aeq.x);

    % Adaptive Equalisation (MEX)
    fprintf('  Adaptive EQ (MEX) ... ');
    tic;
    aeqOut_fxp = adeq_equalize_fxp_mex(mfOut_fxp_fi, ...
                     double(P.SpS), P.AEQ_Eq, ...
                     double(P.AEQ_NTaps), double(P.AEQ_Mu), ...
                     logical(P.AEQ_SingleSpike), ...
                     double(P.AEQ_N1), double(P.AEQ_N2), ...
                     double(P.AEQ_NOut), T_aeq);
    t_aeq_fxp = toc;
    fprintf('%.3f s\n', t_aeq_fxp);

    % Cast AEQ output for VV
    aeqOut_vv_fi = cast(double(aeqOut_fxp), 'like', T_vv.x);

    % VV filter & pilots in fi
    VVFilter_fi = cast(VVFilter, 'like', T_vv.w);
    Pilots_fi   = cast(Pilots,   'like', T_vv.x);

    % Viterbi-Viterbi (MEX)
    fprintf('  VV carrier recovery (MEX) ... ');
    tic;
    vvOut_fxp = cr_viterbiViterbi_fxp_mex(aeqOut_vv_fi, ...
                    double(P.N_pol), double(P.VV_NTaps), ...
                    VVFilter_fi, Pilots_fi, ...
                    double(P.VV_P), double(P.VV_L), ...
                    double(P.VV_CSThreshold), logical(P.VV_UsePilots), ...
                    T_vv);
    t_vv_fxp = toc;
    fprintf('%.3f s\n', t_vv_fxp);

    % Decide & compute BER
    dec_fxp  = qam_decideSymbols(double(vvOut_fxp), P.M, P.N_pol);
    bits_fxp = qam_symbolsToBits(dec_fxp, P.M);
    Ncomp_fxp = min(length(bits_fxp), length(bits));
    BER_fxp   = sum(bits_fxp(1:Ncomp_fxp) ~= bits(1:Ncomp_fxp)) / Ncomp_fxp;
    fprintf('  BER (fxp)   = %.2e  (%d / %d bits)\n', BER_fxp, ...
            sum(bits_fxp(1:Ncomp_fxp) ~= bits(1:Ncomp_fxp)), Ncomp_fxp);

    %% ================================================================
    %  Timing summary
    % =================================================================
    fprintf('\n--- Timing summary ---\n');
    fprintf('  Stage            Float [s]   FXP MEX [s]   Speed-up\n');
    fprintf('  CD EQ            %8.3f    %8.3f      %5.1fx\n', t_cd_fl, t_cd_fxp, t_cd_fl/t_cd_fxp);
    fprintf('  Adaptive EQ      %8.3f    %8.3f      %5.1fx\n', t_aeq_fl, t_aeq_fxp, t_aeq_fl/t_aeq_fxp);
    fprintf('  VV recovery      %8.3f    %8.3f      %5.1fx\n', t_vv_fl, t_vv_fxp, t_vv_fl/t_vv_fxp);
    t_total_fl  = t_cd_fl  + t_aeq_fl  + t_vv_fl;
    t_total_fxp = t_cd_fxp + t_aeq_fxp + t_vv_fxp;
    fprintf('  TOTAL            %8.3f    %8.3f      %5.1fx\n', ...
            t_total_fl, t_total_fxp, t_total_fl/t_total_fxp);

    %% ================================================================
    %  Plots — constellation at each stage (float vs fxp, per pol)
    % =================================================================
    if isfield(P, 'Plot') && P.Plot
        plotTitle = sprintf('%d-QAM  |  SNR %.0f dB  |  %s', P.M, P.SNR_dB, P.FxpConfig);
        ms = 1;  % marker size

        cdOut_fxp_d  = double(cdOut_fxp);
        aeqOut_fxp_d = double(aeqOut_fxp);
        vvOut_fxp_d  = double(vvOut_fxp);

        for p = 1:P.N_pol
            figure('Name', sprintf('Pipeline – Pol %d', p), ...
                   'Position', [50+(p-1)*50, 50, 1600, 900]);

            % Row 1: floating-point --------------------------------------
            subplot(2,5,1);
            plot(real(symbols(:,p)), imag(symbols(:,p)), '.', 'MarkerSize', ms);
            grid on; axis equal; title('TX symbols'); xlabel('I'); ylabel('Q');

            subplot(2,5,2);
            plot(real(rxSig(:,p)), imag(rxSig(:,p)), '.', 'MarkerSize', ms);
            grid on; axis equal; title('After channel'); xlabel('I'); ylabel('Q');

            subplot(2,5,3);
            plot(real(cdOut_fl(:,p)), imag(cdOut_fl(:,p)), '.', 'MarkerSize', ms);
            grid on; axis equal; title('After CD EQ'); xlabel('I'); ylabel('Q');

            subplot(2,5,4);
            plot(real(aeqOut_fl(:,p)), imag(aeqOut_fl(:,p)), '.', 'MarkerSize', ms);
            grid on; axis equal; title('After Adaptive EQ'); xlabel('I'); ylabel('Q');

            subplot(2,5,5);
            plot(real(vvOut_fl(:,p)), imag(vvOut_fl(:,p)), '.', 'MarkerSize', ms);
            grid on; axis equal; title('After VV'); xlabel('I'); ylabel('Q');

            % Row 2: fixed-point (MEX) -----------------------------------
            subplot(2,5,6);
            plot(real(symbols(:,p)), imag(symbols(:,p)), '.', 'MarkerSize', ms);
            grid on; axis equal; title('TX symbols'); xlabel('I'); ylabel('Q');

            subplot(2,5,7);
            plot(real(rxSig(:,p)), imag(rxSig(:,p)), '.', 'MarkerSize', ms);
            grid on; axis equal; title('After channel'); xlabel('I'); ylabel('Q');

            subplot(2,5,8);
            plot(real(cdOut_fxp_d(:,p)), imag(cdOut_fxp_d(:,p)), '.', 'MarkerSize', ms);
            grid on; axis equal; title('After CD EQ (fxp)'); xlabel('I'); ylabel('Q');

            subplot(2,5,9);
            plot(real(aeqOut_fxp_d(:,p)), imag(aeqOut_fxp_d(:,p)), '.', 'MarkerSize', ms);
            grid on; axis equal; title('After Adaptive EQ (fxp)'); xlabel('I'); ylabel('Q');

            subplot(2,5,10);
            plot(real(vvOut_fxp_d(:,p)), imag(vvOut_fxp_d(:,p)), '.', 'MarkerSize', ms);
            grid on; axis equal; title('After VV (fxp)'); xlabel('I'); ylabel('Q');

            sgtitle(sprintf('Pol %d  |  %s  |  BER: float=%.1e  fxp=%.1e', ...
                    p, plotTitle, BER_fl, BER_fxp));
        end
    end

    %% ================================================================
    %  Pack results
    % =================================================================
    R.bits    = bits;
    R.symbols = symbols;
    R.txSig   = txSig;
    R.rxSig   = rxSig;

    % Float path
    R.fl.cdOut   = cdOut_fl;
    R.fl.mfOut   = mfOut_fl;
    R.fl.aeqOut  = aeqOut_fl;
    R.fl.vvOut   = vvOut_fl;
    R.fl.dec     = dec_fl;
    R.fl.bits    = bits_fl;
    R.fl.BER     = BER_fl;
    R.fl.time    = struct('cd', t_cd_fl, 'aeq', t_aeq_fl, 'vv', t_vv_fl);

    % FXP path
    R.fxp.cdOut  = cdOut_fxp;
    R.fxp.mfOut  = mfOut_fxp;
    R.fxp.aeqOut = aeqOut_fxp;
    R.fxp.vvOut  = vvOut_fxp;
    R.fxp.dec    = dec_fxp;
    R.fxp.bits   = bits_fxp;
    R.fxp.BER    = BER_fxp;
    R.fxp.time   = struct('cd', t_cd_fxp, 'aeq', t_aeq_fxp, 'vv', t_vv_fxp);
end
