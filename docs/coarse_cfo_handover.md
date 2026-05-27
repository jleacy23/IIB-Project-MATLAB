# Coarse FD CFO Correction — Debug Handover

## Goal

Add a coarse carrier-frequency-offset (CFO) correction to the two
overlap-save combined blocks
[`src/+eq_clk/combined_cd_fd_gardner_adaptive.m`](../src/+eq_clk/combined_cd_fd_gardner_adaptive.m)
and
[`src/+eq_clk/combined_cd_fd_godard_adaptive.m`](../src/+eq_clk/combined_cd_fd_godard_adaptive.m),
applied in the **frequency domain on the per-block FFT samples after
the CD + matched-filter mask has been multiplied in**. The downstream
timing loop and CMA were observed to degrade with CFO, and the CMA on
its own cannot absorb several-bin CFOs at the 3 GHz CPON worst case.

The correction is supposed to recenter the spectrum to within one FFT
bin so the Modified Godard excess-band window and the downstream
Gardner DPLL see a (near-)zero-CFO signal.

## Algorithm

Per overlap-save block, after `Rfilt = R .* Hstatic` (CD + RRC mask
applied):

1. Compute the magnitude-squared spectrum, summed across polarisations.
2. Estimate the signed FFT-bin centroid:
   $$c[n] = \frac{\sum_k k\,|R_\mathrm{filt}[k]|^2}{\sum_k |R_\mathrm{filt}[k]|^2}$$
   with $k \in \{0, 1, \dots, N/2-1, -N/2, \dots, -1\}$ (natural FFT order).
3. Update the running IIR-filtered estimate:
   $$\hat{c}[n] = (1-\alpha)\,\hat{c}[n-1] + \alpha\,c[n]$$
4. Apply an integer-bin circular shift to recenter:
   `R = circshift(R, -round(\hat{c}), 1)`.

Sub-bin residual CFO is left for the downstream CMA to absorb (CMA is
modulus-only, so it is invariant to any constant phase per symbol).

Bin spacing at NFFT = 128, SpS = 2, Rs = 30.5 GBd:
`binGHz = SpS * Rs / NFFT = 0.477 GHz`. So a 3 GHz CFO is ≈ 6.3 bins,
correctable in principle to within ±0.24 GHz.

## Files Touched

| File | Change |
|---|---|
| [`src/+eq_clk/coarse_cfo_fd.m`](../src/+eq_clk/coarse_cfo_fd.m) | New helper: centroid → IIR → circshift |
| [`src/+eq_clk/combined_cd_fd_godard_adaptive.m`](../src/+eq_clk/combined_cd_fd_godard_adaptive.m) | Calls helper inside the FFT loop, between CD/MF mask and the Godard τ ramp. New `cfoAlpha` arg (default 0.05). |
| [`src/+eq_clk/combined_cd_fd_gardner_adaptive.m`](../src/+eq_clk/combined_cd_fd_gardner_adaptive.m) | Refactored to inline the overlap-save loop (was using `eq_clk.overlap_save_apply`) so the helper can intercept the spectrum. New `cfoAlpha` arg (default 0.05). |
| [`results/combined_eq_clk/combined_eq_clk_sweep.m`](../results/combined_eq_clk/combined_eq_clk_sweep.m) | Exposes `CfoAlpha = 0.05` constant; threaded through `runBlock`. |
| [`tests/test_EqClkCombined.m`](../tests/test_EqClkCombined.m) | New `test_coarse_cfo_fd_accuracy` driving the helper over a sweep of known CFOs (`CfoTruths_GHz = [-3, -1.5, -0.5, 0, 0.5, 1.5, 3]`). |

## Current Symptom

`runtests('test_EqClkCombined/test_coarse_cfo_fd_accuracy')` fails: the
converged bin estimate does not match the true CFO to within one bin.
The exact failure pattern (sign error, bias, magnitude) needs to be
re-run and captured — start by running the test and inspecting the
`truth / est / err` lines that the test prints per CFO point.

This propagates into
`runtests('combined_eq_clk_sweep/test_cfo_sweep')`, where BER stays
high at non-zero CFO because the spectrum never gets recentered.

## Hypotheses to investigate (in priority order)

### 1. Bin-index convention vs FFT layout

`coarse_cfo_fd` assumes `R` is in **natural** FFT order with signed
indices `k_idx = [0:N/2-1, -N/2:-1].'`. Verify:

- Both combined blocks call `fft.fft_flp(InB, false, po2Twiddle)`,
  which should produce natural-order output (matching MATLAB's `fft`).
- The unit test uses MATLAB's `fft(InB, [], 1)` directly, which is
  also natural order.
- `Hstatic` is built as `ifftshift(HCDshift .* HMFshift)` in both
  combined blocks — this is natural order too, so multiplication
  doesn't shift the layout.

Check `fft.fft_flp`'s output ordering convention explicitly — if it
ever reorders to fftshifted form, the centroid estimate uses the wrong
`k_idx` and the sign will be inverted.

### 2. Sign of the centroid vs sign of `lo_freq_shift`

`channel.lo_freq_shift(X, DeltaF, Rs, SpS)` multiplies by
`exp(+1j * 2*pi * deltaF * T * k)` — i.e. a **positive** DeltaF shifts
the spectrum to **positive** frequencies. The centroid we compute is
positive in the same convention. So `circshift(R, -round(centroid))`
should bring it back to DC. Re-verify this sign convention end-to-end
— off-by-one signs will produce 2× CFO error (look for `2 * truth` in
the error column).

### 3. Spectral asymmetry from the matched filter

The RRC matched filter is *symmetric* around DC, so multiplying by
`HMFshift` should not bias the centroid. The CD response
`exp(j * β2/2 * L * ω^2)` is all-pass (`|H_CD| = 1`) and should not
bias `|R|^2` either. Sanity-check by computing the centroid of the
mask alone (i.e. with zero input + AWGN, the centroid should be ≈ 0).

### 4. RRC support is narrower than the FFT band

The signal occupies only `±(1+β)Rs/2 ≈ ±19 GHz` of the `±30.5 GHz`
FFT band at `SpS = 2`. With high CFO, the signal can run off one edge
of the band and wrap (alias) to the other edge. The centroid will
then point in the *wrong* direction. The CPON worst case `±3 GHz`
shifts the band edge from `+19 GHz` to `+22 GHz`, still inside the
`+30.5 GHz` Nyquist, so this should not wrap — but worth confirming
with a synthetic check (set the truth to ±10 GHz and verify the
estimator diverges as expected).

### 5. IIR transient

`alpha = 0.05` means ~20 blocks to settle. The test runs ≈ `Ns *
SpS / (NFFT - NOverlap)` blocks. With `Ns = 2^14 = 16384`, `SpS = 2`,
`NFFT = 256`, `NOverlap = 64`, the loop processes ≈ 170 blocks —
plenty for the IIR to converge. So a transient explanation is
unlikely, but worth verifying by printing `cfoBins` history.

### 6. `pol2cell` reshape in the test

The test code reshapes `rxSig(s:e, :)` into a `NFFT × 1 × 2` array
before calling `fft(InB, [], 1)`. Verify this preserves the per-pol
spectrum layout that `coarse_cfo_fd` expects (the helper indexes
`R(:, 1, p)`). A reshape ordering mistake would scramble polarisations
and corrupt `mag2`.

## How to reproduce

```matlab
addpath(genpath('src'));
runtests('tests/test_EqClkCombined', 'ProcedureName', 'test_coarse_cfo_fd_accuracy')
```

The test prints one line per CFO truth:

```
  truth = +X.XXX GHz, est = +X.XXX GHz, err = +X.XXX GHz  (|err| < 0.477)
```

The first task for the fresh instance is to read that table, identify
the failure pattern (is `est` always wrong by the same constant? wrong
sign? capped at some bin count? bias only at non-zero CFO?), and
match it to one of the hypotheses above.

## What's known to work

- `combined_eq_clk_sweep/test_design_sweep` (CFO = 0) ran the old
  pre-correction blocks fine. The CD/MF + Godard + CMA chain is
  intact.
- The combined blocks both compile and execute — no MATLAB errors,
  just wrong BER. So the wiring (function signatures, call sites,
  state initialisation) is at least surface-correct.

## What's untested

- Whether the residual sub-bin CFO is small enough for CMA to absorb
  even *if* the estimator converges to the correct integer bin. This
  is downstream of the bug above and only matters once the estimator
  is fixed.
- Whether `cfoAlpha = 0.05` is the right time constant — possibly
  worth a quick sweep once accuracy is fixed (the slower the IIR, the
  less the residual jitter but the longer the lock time).
