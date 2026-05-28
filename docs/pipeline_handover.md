# Pipeline FXP Sweep — Handover

## What this test does

Class: `results/pipeline/pipeline_fxp_sweep.m`
Output: `results/pipeline/pipeline_fxp_sweep.mat`

A single `TestCase` (`test_pipeline_fxp_sweep`) runs the **full receiver pipeline** end-to-end at every combination of:

| Axis        | Default values                                                                |
|-------------|--------------------------------------------------------------------------------|
| `po2`       | `[false, true]` — Gardner combined eq with / without power-of-two twiddles    |
| `fr_algo`   | `{differential_kay, fft_search_blind}` — training-aided + blind FR            |
| `FL_eq`     | `FL_eq_vec` (1-D) — fractional bits for `T.Static = T.Clk = T.AdaptEq`        |
| `FL_frcr`   | `FL_frcr_vec` (1-D) — fractional bits for FR `T` and CR `T`                   |
| `SNR_dB`    | `SNR_dB_vec`                                                                  |
| trial       | `1 .. NTrials`                                                                |

`AlgoPo2` and `AlgoFR` are parallel arrays of length 4 (Cartesian-collapsed), then a full grid is done over `FL_eq_vec x FL_frcr_vec` per combo. Word length is `IntBits + FL` with `IntBits = 16`.

### Pipeline

```
modem.modulate(bits) ─┐
                      │
                   pilots/training at ±1±1j (modem.modulate downscales them)
                      │
  modem.rrcPulse → channel.add_chromatic_dispersion → channel.add_pmd (fixed seed)
                  → channel.lo_freq_shift (3 GHz CFO)
                  → channel.apply_timing_error (40 ppm SFO)
                  → channel.add_awgn (swept SNR)
                  → channel.add_phase_noise (1 MHz linewidth)
                      │
                      ▼
      eq_clk.combined_cd_fd_gardner_adaptive_fxp_eq<FL>_mex
         (cfoEnable = true, po2Twiddle = cfg.po2)
                      │
                      ▼
   drop (SUBFRAME_SYMS − NOutAEQ) symbols → aligned to subframe boundary
   trim trailing to integer number of CPON subframes
                      │
                      ▼
   rescale pilot + training positions ×3 (±1±1j → ±3±3j)
                      │
                      ▼
   freq_recovery.<algo>_fxp_fc<FL>_mex
         (training reference also ×3)
                      │
                      ▼
   carrier_recovery.pilots_only_fxp_fc<FL>_mex
         (pilots reference = 3 × repmat(pilotsRef, nSubUsable, 1))
                      │
                      ▼
   computeBER (32 × π/16 phase rotations + pol-swap, mins per pol)
```

### MEX strategy

`TestClassSetup` builds one MEX per unique FL per function (cached on rerun, naming `*_eq<FL>_mex` and `*_fc<FL>_mex`). At runtime each config calls the right MEX via `str2func`.

## Constants worth knowing

```matlab
CFO_GHz   = 3.0     LW_Hz   = 1e6     L_km     = 80     PMD_seed = 12345
SFO_ppm   = 40      Rs      = 30.5    SpS      = 2      Rolloff  = 0.25
NFFT      = 128     NCD     = 22      NLanesGard = 32
NTapsAEQ  = 1       MuAEQ   = 1e-3    N1AEQ    = 500    NOutAEQ  = 1000
                                       SignOnly = true   SingleSpike = true
FR_Nfft = 512   FR_BlindD = 512   MaxFreq = 0.1   TrainingLen = 11
CordicIts = 16  BlockLen_CR = 32
```

PI gains are split per po2:
```matlab
ki_gardner_po2_off = 1e-7    kp_gardner_po2_off = 1e-6
ki_gardner_po2_on  = 1e-6    kp_gardner_po2_on  = 1e-6
```
Both pairs were taken from the existing `combined_eq_clk_sweep` design-sweep results — **tuned at `CFO = 0`**, not 3 GHz.

## Symptom

> BER is still very high for 16 bits of precision and an SNR of 40 dB.

At 16-bit FL and 40 dB SNR a healthy DSP chain should hit BER ≈ 0. The fact that it doesn't means something fundamental is broken — quantisation noise is not the issue.

## Things ruled out (or confirmed sane)

- `cfoEnable = true` is wired through to the equaliser MEX at runtime (line 438 of `pipeline_fxp_sweep.m`). The `B.cfoEnable = false` in `buildEqMex` only sets the codegen prototype to `logical`; runtime value is honoured.
- BER computation resolves π/4 ambiguity (32 × π/16 rotations span the full circle including the four π/2 QPSK rotations) and X↔Y polarisation swap.
- Pilots / training are scaled ×3 in **both** the rx (`eqOut(rescaleIdx,:) = 3*eqOut(rescaleIdx,:)`) **and** the references (`pilotsAll = 3*repmat(pilotsRef, ...)`, `training_fi = cast(3*training, ...)`), so the conj-multiply correlator in `pilots_only_fxp` and the angle subtractor in `differential_kay_fxp` see matched magnitudes.
- `modem.modulate` outputs pilots at `±1±1j` (unit amplitude) — confirmed by reading the source. The ×3 rescale is restoring the CPON-spec `±3±3j` amplitude only at the pipeline stages that expect it.
- The equaliser sees pilots/training **at unit amplitude** (no upstream rescale), so the CMA target isn't disturbed.

## Suspects to investigate next

### 1. (most likely) Time alignment after the Gardner DPLL is unknown

`combined_cd_fd_gardner_adaptive_fxp` runs the **clock recovery loop between the static eq and the adaptive eq**. The Gardner DPLL re-samples the signal with a cubic Farrow interpolator whose NCO converges to the optimum timing offset. The DPLL output's symbol index is **not** simply `Tx symbol index − NOutAEQ` — it's shifted by:
- The DPLL's initial NCO state and `(ki, kp)` lock transient,
- The overlap-save cyclic extension and the static eq's group delay (close to zero for the FD CD/MF, but non-zero for the DPLL convergence).

Our subframe-drop logic assumes `eqOut[1]` corresponds to Tx symbol `NOutAEQ + nDropEq + 1 = 3713` (the start of subframe 2). If the actual shift is ±k symbols, the pilot positions we compute (`(b−1)·32 + 1` within the aligned stream) **don't land on actual pilot symbols** — they land on data symbols. The `pilots_only` CR then forms the phase estimate from `conj(pilot_ref) × random_data_symbol`, which is noise → CR doesn't converge → BER stays at ~0.5.

**Diagnostic**: log `mean(abs(eqOut(pilotIdx, :)))` before the ×3 rescale. Actual pilot symbols (post-eq) sit near magnitude `√2 ≈ 1.4`. If the indices are off, the magnitudes will look like noise (much smaller, no structure).

**Fix candidates**:
- Search over a small lag window (`±32 symbols`) by correlating `abs(eqOut(:,1)) ./ abs(eqOut(:,2))` with the known pilot pattern, and adjust `nDropEq` accordingly.
- Or skip more — use 2 subframes of warmup and search inside the second for the best subframe alignment.
- Or change the equaliser test to use **pilot-aided mode** (`Mode = 1`, `Pilots` non-empty), which fixes a known alignment by design; but this changes the equaliser behaviour and isn't what the upstream sweep tested.

### 2. CFO/MaxFreq margin

`MaxFreq = 0.1` caps the FR estimation range at `0.1 · Rs = 3.05 GHz`. With the 3 GHz CFO, the coarse FFT-centroid CFO correction inside the eq leaves a continuous-valued sub-bin residual; if for any reason the residual sits above ~0.05·Rs the FR's Jacobsen interpolator saturates and the estimate is rubbish.

**Diagnostic**: capture the `cfoBinsApplied` from the eq MEX (second output, currently discarded as `~`) and the FR's returned `frequency_offset` per trial. The sum should be ~3 GHz; the FR's contribution should be small.

**Fix**: raise `MaxFreq` to 0.15 or 0.2 for this sweep, OR widen the eq's coarse CFO accuracy (more averaging blocks in `coarse_cfo_fd`).

### 3. Loop-filter gains tuned at the wrong CFO

`ki_gardner_*` and `kp_gardner_*` defaults are inherited from `combined_eq_clk_sweep`'s design sweep, which tuned at `CFO = 0`. The Gardner DPLL S-curve flattens with CFO (see report eq. 4.14 around line 415 of `report/full/full.tex`); the loop's effective bandwidth at 3 GHz residual is different.

**Fix**: add a one-shot `(ki, kp)` grid search at `CFO = 3 GHz, SNR = SNR_dB_vec(end)` inside `TestClassSetup` (one search per po2), then store and reuse.

### 4. Pol swap and the CR

`pilots_only_fxp` shares one phase estimate **across polarisations**, so a polarisation swap from the equaliser inverts the pilot↔rx pairing on one pol and the per-block phase estimate becomes the average of two unrelated phases. CR collapses, BER stays high.

CMA with `SingleSpike + N1AEQ = 500` rarely swaps, but it does sometimes — worth checking. A quick diagnostic is the per-trial standard deviation of BER across SNR: a working trial shows BER monotone in SNR, a pol-swapped trial shows ~0.5 BER at all SNRs.

### 5. `T.AdaptEq` interpretation

When `EqPrec` is passed as a struct (which is what `pipeline_fxp_sweep` does — `struct('WL', IntBits+FL_eq, 'FL', FL_eq)`), `adaptive_eq.equalize_fxp_types` reads `WL/FL` as the **gradient precision** `T.grad`, holding the data path at HIGH precision (FL = 24/28). That mismatch with `T.Static.x` (which is uniform at the given FL) means the AEQ silently runs at much better precision than the rest of the equaliser. Not a *bug*, but it explains why FL_eq = 4 sweeps don't blow up the AEQ — and it's worth knowing if you're chasing precision-vs-BER curves.

## How to triage in MATLAB (3-step plan)

1. **Verify alignment**. After the subframe-drop, before the ×3 rescale, print
   ```matlab
   pIdx = pipeline_fxp_sweep.pilotTrainingIndices(nSubUsable, P);
   pilotPos = pIdx(P.TrainingLen+1 : P.TrainingLen + P.N_BLOCKS - 1);  % first block, pilots only
   fprintf('|eqOut(pilot)| = %.3f, |eqOut(non-pilot)| = %.3f\n', ...
       mean(abs(eqOut(pilotPos,1))), ...
       mean(abs(eqOut(setdiff(1:size(eqOut,1), pIdx), 1))));
   ```
   If both numbers are similar, the indices are wrong → suspect (1).

2. **Capture the CFO residual**. In `runOneTrial`, change `[yFi, ~]` to `[yFi, cfoBins]` and `[frOut_fi, ~]` to `[frOut_fi, frHz]`. Print `cfoBins`, `frHz`, and their sum (should be ≈ `CFO_GHz * 1e9`).

3. **Try a smaller CFO**. Set `CFO_GHz = 0.5` and rerun one trial at FL_eq=FL_frcr=12, SNR=30. If BER drops to ~0, the problem is CFO/loop-filter-tuning (suspects 2/3). If BER stays ~0.5, the problem is alignment (suspect 1).

## File pointers

- Test: `results/pipeline/pipeline_fxp_sweep.m`
- Equaliser fxp: `src/+eq_clk/combined_cd_fd_gardner_adaptive_fxp.m`
- FR fxp: `src/+freq_recovery/differential_kay_fxp.m`, `src/+freq_recovery/fft_search_fxp.m`
- CR fxp: `src/+carrier_recovery/pilots_only_fxp.m`
- CPON modulator: `src/+modem/modulate.m` (pilots ARE at ±1±1j, not ±3±3j — be sure)
- Comparable sweep without the FR/CR stages: `results/combined_eq_clk/combined_eq_clk_sweep.m`
- Comparable sweep without the equaliser stage: `results/frequency_recovery/bit_width_full.m`
- Pipeline parameter doc: relevant section of `report/full/full.tex` around line 403.
