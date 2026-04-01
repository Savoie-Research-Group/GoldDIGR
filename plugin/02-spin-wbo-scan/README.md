# Stage 2: Two-Pass Spin/WBO Scan

## Overview

Scans each IRC frame for the optimal spin state and extracts Wiberg bond orders.

## Usage

```bash
python spin_wbo_scan.py \
    --trj finished_irc.trj \
    --charge 0 \
    --mult -1 \
    --uhf-max 6 \
    --workdir xTB-scan \
    --nprocs 4
```

## Two-Pass Architecture

**Pass 1 (tblite):** For each IRC frame, compute single-point energies at UHF = 0, 2, 4, 6
(or 1, 3, 5 for odd-electron systems). Select the lowest-energy spin state per frame.
Uses `tblite` for efficient spin-polarized calculations.

**Pass 2 (standard xTB):** Re-run at the best UHF per frame using standard `xtb` to extract
Mulliken charges and Wiberg bond orders from the `.wbo` output file.

## Outputs

| File | Content |
|------|---------|
| `spin_scan_summary.csv` | Per-frame energies at all UHF values + best UHF |
| `wbo_timeseries.csv` | WBO values for all atom pairs across frames |
| `spin_crossover_frames.txt` | Frames where the ground-state spin multiplicity changes |
| `bem_snapshots/step_NNNN.json` | BEM-lite JSON per frame |
| `frame_NNNNN_bond_electrons.csv` | Full BEM CSV per frame |
| `reactive_summary.csv` | Time series of entries with |Δ| ≥ 0.5 between first/last frames |

## Helper Scripts

- `run_analysis.sh` — Shell wrapper for batch execution
- `compute_sco_wbo_stats.py` — Aggregate spin-crossover and WBO slope statistics
- `summarize_sco_wbo.py` — Summary report across multiple reactions
