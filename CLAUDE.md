# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MATLAB implementation of a coherent optical receiver DSP pipeline for DP-QPSK signals, targeting fixed-point hardware (IIB project). The codebase has dual floating-point and fixed-point paths — floating-point (`*.m`) for algorithm development, fixed-point (`*_fxp.m`) compiled to MEX via MATLAB Coder for hardware-efficient implementations.

## Commands

**Running tests (in MATLAB):**
```matlab
% Run all tests
runtests('tests')

% Run a single test file
runtests('tests/test_CarrierRecovery')

% Run a specific test method
results = runtests('tests/test_FreqRecovery', 'ProcedureName', 'test_fft_search');
```

**Building MEX (fixed-point) binaries:**
```matlab
% From repo root, run any build script in build/
run('build/build_carrier_recovery_bps_fxp_mex.m')
run('build/build_freq_recovery_fft_search_fxp_mex.m')
% etc. — each build script generates a MEX for the corresponding *_fxp.m function
```

**Adding src to path (required before running anything):**
```matlab
addpath(genpath('src'))
```

## Architecture

### DSP Pipeline (transmitter → channel → receiver)

```
modulate()          → DP-QPSK symbols [Nsym×2], pilots, training
channel/*           → chromatic dispersion, phase noise, AWGN, PMD, freq shift
─── Receiver ───────────────────────────────────────────────────────────
freq_recovery/*     → coarse frequency offset correction
clk_recovery/*      → symbol timing recovery + interpolation
cd_eq/*             → frequency-domain chromatic dispersion equalization
adaptive_eq/*       → butterfly CMA/RDE adaptive equalizer
carrier_recovery/*  → fine phase recovery → decideSymbols() → symbolsToBits()
```

### Module layout (`src/`)

Each DSP stage is a MATLAB package (`+name/`). Every stage has:
- A floating-point implementation (algorithm development/verification)
- A `*_fxp.m` fixed-point version (MATLAB Coder compatible, uses `fi` objects)
- A `*_fxp_types.m` or `fxp_types.m` defining bit-width configurations

| Package | Purpose |
|---|---|
| `+modem` | Modulation, demodulation, pulse shaping, CPON framing |
| `+freq_recovery` | Frequency offset estimation: `fft_search`, `differential_kay`, `tretter_kay` |
| `+carrier_recovery` | Phase recovery: `bps` (blind phase search), `viterbiViterbi`, `pilots_only` |
| `+cd_eq` | Frequency-domain chromatic dispersion equalization (overlap-save FFT) |
| `+adaptive_eq` | CMA/RDE butterfly equalizer |
| `+clk_recovery` | Timing estimation and interpolation |
| `+channel` | Channel impairment models (AWGN, CD, phase noise, PMD, ADC) |
| `+fft` | Radix-2 fixed-point FFT/IFFT (`fft_fxp.m`) |
| `+energy` | Operation counting for energy modeling |

### CPON Frame Structure

Subframes are 3712 symbols = 116 × 32-symbol blocks. Each block has:
- 1 pilot symbol (PRBS10, seeds `0x19E` / `0x0D0` for X/Y polarizations, amplitude `±3±3j`)
- 31 data symbols

Training sequence: 11 fixed symbols prepended per subframe (X and Y polarizations differ). See [docs/cpon_framing_structure.md](docs/cpon_framing_structure.md) for full specification.

### Signal Dimensions

- Symbols: `[Nsym × 2]` complex (column 1 = X-pol, column 2 = Y-pol)
- Fixed-point types use MATLAB `fi` objects; word/fractional lengths are configured per-module in `fxp_types.m`

### Fixed-Point / MEX Pattern

Each `*_fxp.m` file is MATLAB Coder–compatible (no dynamic allocation, typed inputs). Build scripts in `build/` call `codegen` to produce MEX binaries. When modifying fixed-point code, verify MATLAB Coder codegen restrictions apply (no `nargin`/`nargout` checks, no cell arrays, explicit type casting).
