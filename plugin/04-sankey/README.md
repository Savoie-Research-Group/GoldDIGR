# Stage 4: Sankey Diagram Generation

## Overview

Builds electron-flow Sankey diagrams from the IRC BEM time series,
decomposing electron redistribution into transport (atom → atom)
and storage (bond filling/emptying) components.

## Usage

```bash
python irc_sankey_zoomsvg_fixed.py \
    --root IRC_Analysis/finished_irc \
    --window 50 \
    --topK 999 \
    --min-frac 0.0 \
    --out IRC_Analysis/sankey
```

## Outputs

| File | Content |
|------|---------|
| `sankey_transport.html` | Interactive transport Sankey (atom → atom flows) |
| `sankey_two_rail.html` | Transport + bond-storage two-rail Sankey |
| `nodes.csv` / `links_*.csv` | Raw data for reproducibility |
| `event_times.csv` | Per-bond-change event timing relative to TS |
| `*.svg` | Exported SVG with zoomed variants (if kaleido available) |

## Method

- Ownership electrons: R_i(k) = U_i(k) + Σ_j V_ij(k) (diagonal + incident bonds)
- Transport flows f satisfy: B·f = ΔR (minimum-norm least-squares solution)
- Integrated over a TS-centered window (default ±50 frames)
- Links below `--min-frac` of total flow are suppressed
