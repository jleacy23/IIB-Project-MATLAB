# CPON Downstream PMD Layer — Framing Structure

**Based on:** CPON-SP-PMDv1.0-I01-251216  
**Scope:** Downstream only. Subframe is the primary organisational unit. No shortened subframe.

---

## Table of Contents

1. [Modulation and Symbol Mapping](#1-modulation-and-symbol-mapping)
2. [DSP Subframe Structure](#2-dsp-subframe-structure)
3. [Training Symbol Sequence](#3-training-symbol-sequence)
4. [Pilot Symbol Generation](#4-pilot-symbol-generation)
5. [Pilot Insertion Rule](#5-pilot-insertion-rule)
6. [Complete Pilot Sequence Table](#6-complete-pilot-sequence-table)
7. [Overhead Budget](#7-overhead-budget)

---

## 1. Modulation and Symbol Mapping

### 1.1 Modulation Format

- **Modulation:** DP-QPSK (dual-polarization, non-differential)
- **Symbol rate:** 30.504432 Gbaud (±20 ppm)
- **Bits per symbol:** 4 aggregate (2 per polarization)
- **PMD line rate:** 122.017728 Gbps

### 1.2 QPSK Constellation

| Bit pair (I-bit, Q-bit) | I amplitude | Q amplitude |
|------------------------|-------------|-------------|
| (0, 0)                 | −1          | −1          |
| (0, 1)                 | −1          | +1          |
| (1, 0)                 | +1          | −1          |
| (1, 1)                 | +1          | +1          |

> Training and pilot symbols use amplitude values of ±3, placing them outside the data constellation and making them unambiguously identifiable.

### 1.3 Bit-to-Symbol Distribution

For symbol index `i`, bits from the TC layer bitstream are distributed as follows:

| Bit index | Mapped to            |
|-----------|----------------------|
| c(4i)     | I component of X-pol |
| c(4i+2)   | Q component of X-pol |
| c(4i+1)   | I component of Y-pol |
| c(4i+3)   | Q component of Y-pol |

---

## 2. DSP Subframe Structure

### 2.1 Parameters

| Parameter | Value |
|-----------|-------|
| Total symbols per subframe | 3,712 |
| Number of 32-symbol blocks | 116 |
| Training symbols | 11 (block 1 only) |
| Pilot symbols | 116 (one per block) |
| First training symbol (TS1) | Shared as pilot index 1 |

```
3,712 = 116 × 32  ✓
```

### 2.2 Symbol-Level Layout

The subframe consists of 116 consecutive 32-symbol blocks. Block 1 is special due to the training sequence; all subsequent blocks follow a simple data + pilot pattern.

```
BLOCK 1 (symbols 1–32):
┌────┬────┬────┬────┬────┬────┬────┬────┬────┬────┬────┬─────────────────────┐
│TS1 │TS2 │TS3 │TS4 │TS5 │TS6 │TS7 │TS8 │TS9 │TS10│TS11│  DATA × 21         │
│(PS)│    │    │    │    │    │    │    │    │    │    │                     │
└────┴────┴────┴────┴────┴────┴────┴────┴────┴────┴────┴─────────────────────┘
  [1]  [2]  [3]  [4]  [5]  [6]  [7]  [8]  [9] [10] [11]  [12] ......... [32]
   ↑
   TS1 is simultaneously pilot symbol index 1

BLOCKS 2–116 (each 32 symbols):
┌────────────────────────────────────────────────────────────┬────┐
│                       DATA × 31                            │ PS │
└────────────────────────────────────────────────────────────┴────┘
 [1] .................................................... [31] [32]
```

### 2.3 Pilot Symbol Positions

Pilot symbols occupy the first position of every 32-symbol block. In absolute symbol indices within the subframe:

```
Pilot positions:  1, 33, 65, 97, 128, ..., 3,681, 3,712
General rule:     position = 1 + 32k,  for k = 0, 1, 2, ..., 115
```

### 2.4 Symbol Accounting

| Symbol type | Count | Percentage |
|-------------|-------|------------|
| Pilot symbols (incl. TS1) | 116 | 3.13% |
| Training symbols (TS2–TS11) | 10 | 0.27% |
| Data symbols | 3,586 | 96.60% |
| **Total** | **3,712** | **100%** |

---

## 3. Training Symbol Sequence

### 3.1 Values

The training sequence consists of 11 fixed symbols transmitted at the start of every subframe. X-pol and Y-pol carry different sequences:

| Index | X-pol value | Y-pol value | Notes |
|-------|-------------|-------------|-------|
| 1     | −3 + 3j     | −3 − 3j     | Also pilot symbol index 1 |
| 2     | +3 + 3j     | −3 − 3j     | |
| 3     | −3 + 3j     | +3 − 3j     | |
| 4     | +3 + 3j     | −3 + 3j     | |
| 5     | −3 − 3j     | −3 + 3j     | |
| 6     | +3 + 3j     | +3 + 3j     | |
| 7     | −3 − 3j     | −3 − 3j     | |
| 8     | −3 − 3j     | −3 + 3j     | |
| 9     | +3 + 3j     | +3 − 3j     | |
| 10    | +3 − 3j     | +3 + 3j     | |
| 11    | +3 − 3j     | +3 − 3j     | |

> Where j = √(−1).

### 3.2 Properties

- X-pol and Y-pol carry **different sequences** — the asymmetry allows the receiver to resolve polarization swap ambiguity, which pilots alone cannot do.
- The sequence is **identical at every subframe boundary** — it is fixed, not PRBS-generated.
- **TS1 matches the first output of the pilot PRBS at seed reset**, allowing it to serve simultaneously as training symbol 1 and pilot index 1 without discontinuity.

---

## 4. Pilot Symbol Generation

### 4.1 PRBS Generator

Pilot symbols are generated from a **PRBS10** linear feedback shift register:

| Parameter | Value |
|-----------|-------|
| Generator polynomial | x¹⁰ + x⁸ + x⁴ + x³ + 1 |
| Register length | 10 bits |
| X-polarization seed | 0x19E |
| Y-polarization seed | 0x0D0 |
| Sequence reset | At the start of every subframe |
| Maximum sequence length | 2¹⁰ − 1 = 1,023 bits |

### 4.2 Seed Selection Criteria

The seeds satisfy three constraints:

1. **TS1 consistency** — the first PRBS output matches the TS1 value, so TS1 can legitimately serve as both training symbol and pilot index 1.
2. **DC balance** — seeds produce approximately equal numbers of 0s and 1s, giving the pilot sequence good spectral and correlation properties.
3. **X/Y decorrelation** — different seeds on each polarization produce uncorrelated sequences, which aids polarization demultiplexing in the receiver DSP.

### 4.3 PRBS-to-QPSK Mapping

Consecutive PRBS output bit pairs are mapped to QPSK symbols using the same Table 4 mapping as data, then scaled by 3:

```
PRBS bits:   b0  b1 | b2  b3 | b4  b5 | ...
                 ↓       ↓       ↓
Pilot sym:  sym 1   sym 2   sym 3  ...

For symbol k:
  I = +3 if b(2k)   = 1,  else −3
  Q = +3 if b(2k+1) = 1,  else −3
```

### 4.4 Reset Behaviour

The PRBS resets to its seed at the start of every subframe. Therefore:

- The pilot sequence is **identical in every subframe** for a given polarization.
- The receiver always knows the exact expected value at every pilot position.
- There is no pilot sequence accumulation or drift across subframe boundaries.

---

## 5. Pilot Insertion Rule

### 5.1 Rule

```
Insert a pilot symbol at every position p where:

    p = 1 + 32k,   k = 0, 1, 2, ..., 115

counting from symbol 1 of the subframe.
```

The pilot symbol at position `p` is the `(k+1)`th output of the PRBS for that polarization (after seed reset at the subframe boundary).

### 5.2 Pilot Index Tracking

| Symbol position | Pilot index | Notes |
|----------------|-------------|-------|
| 1              | 1           | = TS1 |
| 33             | 2           | |
| 65             | 3           | |
| 97             | 4           | |
| ...            | ...         | |
| 3,681          | 116         | |

### 5.3 Annotated Subframe Symbol Map

```
Pos   1: [PS/TS1]   ← Pilot index 1, Training symbol 1
Pos   2: [TS2]      ← Training symbol 2
Pos   3: [TS3]
Pos   4: [TS4]
Pos   5: [TS5]
Pos   6: [TS6]
Pos   7: [TS7]
Pos   8: [TS8]
Pos   9: [TS9]
Pos  10: [TS10]
Pos  11: [TS11]
Pos  12–32:  [DATA × 21]
Pos  33: [PS]       ← Pilot index 2
Pos  34–64:  [DATA × 31]
Pos  65: [PS]       ← Pilot index 3
Pos  66–96:  [DATA × 31]
     ...
Pos 3,681: [PS]     ← Pilot index 116
Pos 3,682–3,712:  [DATA × 31]
```

---

## 6. Complete Pilot Sequence Table

All 116 pilot symbol values for a full subframe. Values are complex amplitudes with I, Q ∈ {−3, +3}.

| Index | Pilot X | Pilot Y | Index | Pilot X | Pilot Y | Index | Pilot X | Pilot Y |
|-------|---------|---------|-------|---------|---------|-------|---------|---------|
| 1     | −3+3j   | −3−3j   | 41    | +3−3j   | +3−3j   | 81    | +3+3j   | +3−3j   |
| 2     | +3+3j   | −3−3j   | 42    | +3−3j   | +3−3j   | 82    | −3−3j   | −3−3j   |
| 3     | +3−3j   | +3−3j   | 43    | −3+3j   | −3−3j   | 83    | −3−3j   | +3+3j   |
| 4     | −3+3j   | +3+3j   | 44    | −3+3j   | −3−3j   | 84    | +3+3j   | −3−3j   |
| 5     | +3−3j   | −3−3j   | 45    | −3−3j   | +3+3j   | 85    | +3−3j   | −3−3j   |
| 6     | +3−3j   | +3+3j   | 46    | −3+3j   | −3+3j   | 86    | −3+3j   | −3−3j   |
| 7     | −3−3j   | −3+3j   | 47    | −3−3j   | +3+3j   | 87    | +3+3j   | +3−3j   |
| 8     | +3+3j   | −3+3j   | 48    | +3+3j   | −3+3j   | 88    | +3−3j   | −3+3j   |
| 9     | −3+3j   | −3−3j   | 49    | +3+3j   | +3−3j   | 89    | −3−3j   | −3+3j   |
| 10    | +3+3j   | +3+3j   | 50    | −3+3j   | −3+3j   | 90    | +3−3j   | +3−3j   |
| 11    | +3+3j   | +3+3j   | 51    | +3−3j   | +3+3j   | 91    | +3−3j   | +3+3j   |
| 12    | −3−3j   | −3−3j   | 52    | +3−3j   | −3+3j   | 92    | −3+3j   | +3−3j   |
| 13    | +3+3j   | +3+3j   | 53    | +3−3j   | −3+3j   | 93    | −3−3j   | +3−3j   |
| 14    | +3−3j   | +3+3j   | 54    | −3−3j   | +3+3j   | 94    | +3+3j   | −3+3j   |
| 15    | +3+3j   | +3−3j   | 55    | +3−3j   | −3+3j   | 95    | −3−3j   | +3−3j   |
| 16    | +3−3j   | +3+3j   | 56    | +3+3j   | −3+3j   | 96    | −3−3j   | +3−3j   |
| 17    | +3+3j   | +3+3j   | 57    | −3+3j   | −3−3j   | 97    | +3+3j   | −3+3j   |
| 18    | +3−3j   | −3+3j   | 58    | −3−3j   | +3−3j   | 98    | −3+3j   | +3−3j   |
| 19    | −3+3j   | −3−3j   | 59    | +3−3j   | +3−3j   | 99    | +3−3j   | −3−3j   |
| 20    | −3−3j   | +3−3j   | 60    | +3+3j   | −3+3j   | 100   | −3−3j   | +3+3j   |
| 21    | +3+3j   | +3−3j   | 61    | +3−3j   | +3+3j   | 101   | +3+3j   | −3−3j   |
| 22    | −3+3j   | +3+3j   | 62    | −3−3j   | −3−3j   | 102   | −3+3j   | −3+3j   |
| 23    | −3+3j   | −3+3j   | 63    | +3−3j   | +3+3j   | 103   | −3−3j   | −3+3j   |
| 24    | +3−3j   | +3−3j   | 64    | −3+3j   | −3+3j   | 104   | −3−3j   | +3+3j   |
| 25    | −3+3j   | +3−3j   | 65    | +3−3j   | +3−3j   | 105   | +3+3j   | −3+3j   |
| 26    | −3+3j   | +3+3j   | 66    | +3+3j   | −3+3j   | 106   | −3−3j   | +3−3j   |
| 27    | −3+3j   | −3+3j   | 67    | +3−3j   | −3−3j   | 107   | +3+3j   | +3+3j   |
| 28    | −3+3j   | +3+3j   | 68    | −3+3j   | +3−3j   | 108   | −3+3j   | −3+3j   |
| 29    | −3−3j   | +3+3j   | 69    | +3−3j   | −3+3j   | 109   | −3−3j   | +3+3j   |
| 30    | +3−3j   | +3−3j   | 70    | −3+3j   | −3+3j   | 110   | −3+3j   | −3−3j   |
| 31    | −3−3j   | −3+3j   | 71    | +3+3j   | −3+3j   | 111   | −3−3j   | −3+3j   |
| 32    | +3+3j   | −3−3j   | 72    | −3−3j   | −3−3j   | 112   | −3+3j   | +3−3j   |
| 33    | −3+3j   | +3−3j   | 73    | −3−3j   | −3+3j   | 113   | −3+3j   | −3+3j   |
| 34    | −3+3j   | −3−3j   | 74    | +3−3j   | +3+3j   | 114   | +3+3j   | +3+3j   |
| 35    | −3+3j   | −3−3j   | 75    | −3+3j   | −3−3j   | 115   | +3+3j   | +3−3j   |
| 36    | +3−3j   | +3−3j   | 76    | +3−3j   | −3−3j   | 116   | −3−3j   | +3−3j   |
| 37    | +3−3j   | +3−3j   | 77    | −3+3j   | −3−3j   |       |         |         |
| 38    | −3−3j   | −3−3j   | 78    | −3−3j   | +3+3j   |       |         |         |
| 39    | −3−3j   | +3+3j   | 79    | +3+3j   | −3−3j   |       |         |         |
| 40    | +3−3j   | −3−3j   | 80    | +3+3j   | −3−3j   |       |         |         |

---

## 7. Overhead Budget

### 7.1 Per Subframe

| Symbol type | Count | Percentage |
|-------------|-------|------------|
| Pilot symbols (incl. TS1) | 116 | 3.13% |
| Training symbols (TS2–TS11) | 10 | 0.27% |
| Data symbols | 3,586 | 96.60% |
| **Total** | **3,712** | **100%** |

### 7.2 Line Rate Chain

```
100 Gbps        Net user data rate
     ↓  + TC layer overhead (FEC parity, PSBd framing)
99.5328 Gbps    TC layer line rate
     ↓  + PMD overhead (pilots 3.13%, training 0.27%)
122.017728 Gbps PMD line rate  (= 4 × 30.504432 Gbaud)
```

---

*Based on CPON-SP-PMDv1.0-I01-251216. Downstream only. Subframe is the primary structural unit.*
