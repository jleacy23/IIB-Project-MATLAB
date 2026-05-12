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

    Nbits = 4 * P.Ns;           % CPON: 4 data bits per dual-pol symbol; modulate pads to integer subframes

    %% ================================================================
    %  Transmitter (shared by both paths)
    % =================================================================
    bits    = modem.randomBits(Nbits);                               % [Nbits x 1]
    [symbols, ~, ~, ~] = modem.modulate(bits);                       % [Ns_actual x 2]
    Ns_actual = size(symbols, 1);
    fprintf('TX: DP-QPSK, %d pol, %d symbols/pol\n', P.N_pol, Ns_actual);
    txSig   = modem.rrcPulse(symbols, P.SpS, P.Rolloff, P.Span);

    %% ================================================================
    %  Channel
    % =================================================================
    fprintf('Channel: SNR=%.0f dB, D=%.0f ps/(nm·km), L=%.0f km, LW=%.0f kHz\n', ...
            P.SNR_dB, P.D, P.L, P.LW/1e3);

    rxSig = channel.add_awgn(txSig, P.SNR_dB);
    rxSig = channel.add_chromatic_dispersion(rxSig, P.L, P.SpS, P.Rs, P.D, P.CWL);
    rxSig = channel.add_phase_noise(rxSig, P.Rs, P.LW);
    rxSig = channel.add_pmd(rxSig, P.L, P.SpS, P.Rs, P.DGDSpec, P.N_pmd);
    rxSig = channel.adc(rxSig, P.ENOBits);

    %% ================================================================
    %  Generate VV filter (shared)
    % =================================================================
    SymbolEnergy = 1;   % unit-power constellation
    VVFilter = carrier_recovery.genVVFilter(P.LW, P.Rs, P.SNR_dB, SymbolEnergy, ...
                              P.N_pol, P.VV_NTaps);

    % The adaptive EQ discards NOut symbols from the front, so the VV
    % input corresponds to TX symbols (NOut+1 : NOut+Nsym_out).
    symOffset = P.AEQ_NOut;   % TX symbol index offset

    %% ================================================================
    %  PATH A — Floating-point (double)
    % =================================================================
    fprintf('\n--- Floating-point path ---\n');

    % CD Equalisation
    fprintf('  CD EQ ... ');
    tic;
    cdOut_fl = cd_eq.equalize(rxSig, P.D, P.L, P.CWL, P.Rs, ...
                             P.N_pol, P.SpS, P.NFFT);
    t_cd_fl = toc;
    fprintf('%.3f s\n', t_cd_fl);

    % Matched filter
    mfOut_fl = modem.matched_filter(cdOut_fl, P.SpS, 'rrc', P.Rolloff, P.Span);

    % Adaptive Equalisation
    fprintf('  Adaptive EQ ... ');
    tic;
    aeqOut_fl = adaptive_eq.equalize(mfOut_fl, P.SpS, ...
                              P.AEQ_NTaps, P.AEQ_Mu, P.AEQ_SingleSpike, ...
                              P.AEQ_N1, P.AEQ_NOut, P.AEQ_SignOnly);
    t_aeq_fl = toc;
    fprintf('%.3f s\n', t_aeq_fl);

    % Viterbi-Viterbi
    fprintf('  VV carrier recovery ... ');
    tic;
    vvOut_fl = carrier_recovery.viterbiViterbi(aeqOut_fl, P.N_pol, P.VV_NTaps, ...
                                 VVFilter);
    t_vv_fl = toc;
    fprintf('%.3f s\n', t_vv_fl);

    % Decide & compute BER — resolve pi/2 ambiguity + pol swap
    Nsym_fl  = size(vvOut_fl, 1);
    refEnd   = min(symOffset + Nsym_fl, Ns_actual);
    Nuse_fl  = refEnd - symOffset;
    refSym_fl = symbols(symOffset+1 : refEnd, :);
    rotations = [1, 1j, -1, -1j];

    totalErr_fl  = 0;
    totalBits_fl = 0;
    bestRotPerPol_fl = zeros(1, P.N_pol);
    bestSrcPerPol_fl = zeros(1, P.N_pol);
    for p = 1:P.N_pol
        refBitsPol = modem.symbolsToBits(refSym_fl(:,p));
        bestPolBER = Inf;
        for q = 1:P.N_pol          % try both EQ outputs (pol swap)
            for ri = 1:4           % try all rotations
                vvRot   = vvOut_fl(1:Nuse_fl, q) * rotations(ri);
                decRot  = modem.decideSymbols(vvRot);
                bitsRot = modem.symbolsToBits(decRot);
                polBER  = sum(bitsRot ~= refBitsPol) / numel(refBitsPol);
                if polBER < bestPolBER
                    bestPolBER = polBER;
                    bestRotPerPol_fl(p) = ri;
                    bestSrcPerPol_fl(p) = q;
                end
            end
        end
        totalErr_fl  = totalErr_fl  + bestPolBER * numel(refBitsPol);
        totalBits_fl = totalBits_fl + numel(refBitsPol);
    end
    BER_fl = totalErr_fl / totalBits_fl;

    % Apply best rotation per pol for constellation plots
    for p = 1:P.N_pol
        vvOut_fl(:,p) = vvOut_fl(:, bestSrcPerPol_fl(p)) * rotations(bestRotPerPol_fl(p));
    end
    dec_fl  = modem.decideSymbols(vvOut_fl(1:Nuse_fl, :));
    bits_fl = modem.symbolsToBits(dec_fl);
    fprintf('  BER (float) = %.2e\n', BER_fl);

    %% ================================================================
    %  PATH B — Fixed-point (MEX)
    % =================================================================
    fprintf('\n--- Fixed-point MEX path (CD=%s, AEQ=%s, VV=%s) ---\n', ...
            P.FxpConfig_CD, P.FxpConfig_AEQ, P.FxpConfig_VV);

    % Load types tables
    T_cd  = cd_eq.equalize_fxp_types(P.FxpConfig_CD);
    T_aeq = adaptive_eq.equalize_fxp_types(P.FxpConfig_AEQ);
    T_vv  = carrier_recovery.viterbiViterbi_fxp_types(P.FxpConfig_VV);

    % Cast channel output to fi for the fixed-point path
    rxSig_fi = cast(rxSig, 'like', T_cd.x);

    % CD Equalisation (MEX)
    fprintf('  CD EQ (MEX) ... ');
    tic;
    cdOut_fxp = cd_eq.equalize_fxp_mex(rxSig_fi, ...
                    double(P.D), double(P.L), double(P.CWL), ...
                    double(P.Rs), double(P.N_pol), double(P.SpS), ...
                    double(P.NFFT), logical(P.po2Twiddle), T_cd);
    t_cd_fxp = toc;
    fprintf('%.3f s\n', t_cd_fxp);

    % Matched filter (float — not fixed-point)
    mfOut_fxp = modem.matched_filter(double(cdOut_fxp), P.SpS, 'rrc', ...
                                   P.Rolloff, P.Span);

    % Cast back to fi for adaptive EQ
    mfOut_fxp_fi = cast(mfOut_fxp, 'like', T_aeq.x);

    % Adaptive Equalisation (MEX)
    fprintf('  Adaptive EQ (MEX) ... ');
    tic;
    aeqOut_fxp = adaptive_eq.equalize_fxp_mex(mfOut_fxp_fi, ...
                     double(P.SpS), ...
                     double(P.AEQ_NTaps), double(P.AEQ_Mu), ...
                     logical(P.AEQ_SingleSpike), ...
                     double(P.AEQ_N1), ...
                     double(P.AEQ_NOut), logical(P.AEQ_SignOnly), T_aeq);
    t_aeq_fxp = toc;
    fprintf('%.3f s\n', t_aeq_fxp);

    % Cast AEQ output for VV
    aeqOut_vv_fi = cast(double(aeqOut_fxp), 'like', T_vv.x);

    % VV filter in fi
    VVFilter_fi = cast(VVFilter, 'like', T_vv.w);

    % Viterbi-Viterbi (MEX)
    fprintf('  VV carrier recovery (MEX) ... ');
    tic;
    vvOut_fxp = carrier_recovery.viterbiViterbi_fxp_mex(aeqOut_vv_fi, ...
                    double(P.N_pol), double(P.VV_NTaps), ...
                    VVFilter_fi, T_vv);
    t_vv_fxp = toc;
    fprintf('%.3f s\n', t_vv_fxp);

    % Decide & compute BER — resolve pi/2 ambiguity + pol swap
    vvOut_fxp_d = double(vvOut_fxp);
    Nsym_fxp    = size(vvOut_fxp_d, 1);
    refEnd_fxp  = min(symOffset + Nsym_fxp, Ns_actual);
    Nuse_fxp    = refEnd_fxp - symOffset;
    refSym_fxp  = symbols(symOffset+1 : refEnd_fxp, :);

    totalErr_fxp  = 0;
    totalBits_fxp = 0;
    bestRotPerPol_fxp = zeros(1, P.N_pol);
    bestSrcPerPol_fxp = zeros(1, P.N_pol);
    for p = 1:P.N_pol
        refBitsPol = modem.symbolsToBits(refSym_fxp(:,p));
        bestPolBER = Inf;
        for q = 1:P.N_pol          % try both EQ outputs (pol swap)
            for ri = 1:4           % try all rotations
                vvRot   = vvOut_fxp_d(1:Nuse_fxp, q) * rotations(ri);
                decRot  = modem.decideSymbols(vvRot);
                bitsRot = modem.symbolsToBits(decRot);
                polBER  = sum(bitsRot ~= refBitsPol) / numel(refBitsPol);
                if polBER < bestPolBER
                    bestPolBER = polBER;
                    bestRotPerPol_fxp(p) = ri;
                    bestSrcPerPol_fxp(p) = q;
                end
            end
        end
        totalErr_fxp  = totalErr_fxp  + bestPolBER * numel(refBitsPol);
        totalBits_fxp = totalBits_fxp + numel(refBitsPol);
    end
    BER_fxp = totalErr_fxp / totalBits_fxp;

    % Apply best rotation per pol for constellation plots
    for p = 1:P.N_pol
        vvOut_fxp_d(:,p) = vvOut_fxp_d(:, bestSrcPerPol_fxp(p)) * rotations(bestRotPerPol_fxp(p));
    end
    dec_fxp  = modem.decideSymbols(vvOut_fxp_d(1:Nuse_fxp, :));
    bits_fxp = modem.symbolsToBits(dec_fxp);
    fprintf('  BER (fxp)   = %.2e\n', BER_fxp);

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
        plotTitle = sprintf('QPSK  |  SNR %.0f dB  |  CD=%s AEQ=%s VV=%s', ...
                    P.SNR_dB, P.FxpConfig_CD, P.FxpConfig_AEQ, P.FxpConfig_VV);
        ms = 1;  % marker size

        cdOut_fxp_d  = double(cdOut_fxp);
        aeqOut_fxp_d = double(aeqOut_fxp);

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
    R.fl.rotPerPol = bestRotPerPol_fl - 1;   % 0..3 = multiples of pi/2
    R.fl.srcPerPol = bestSrcPerPol_fl;        % which EQ output matched each ref pol
    R.fl.time    = struct('cd', t_cd_fl, 'aeq', t_aeq_fl, 'vv', t_vv_fl);

    % FXP path
    R.fxp.cdOut  = cdOut_fxp;
    R.fxp.mfOut  = mfOut_fxp;
    R.fxp.aeqOut = aeqOut_fxp;
    R.fxp.vvOut  = vvOut_fxp_d;         % already rotated, double
    R.fxp.dec    = dec_fxp;
    R.fxp.bits   = bits_fxp;
    R.fxp.BER    = BER_fxp;
    R.fxp.rotPerPol = bestRotPerPol_fxp - 1;
    R.fxp.srcPerPol = bestSrcPerPol_fxp;
    R.fxp.time   = struct('cd', t_cd_fxp, 'aeq', t_aeq_fxp, 'vv', t_vv_fxp);
end
