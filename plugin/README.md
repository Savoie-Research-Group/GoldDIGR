# GoldDIGR Computational Plugin

Reproducible workflow scripts for mining organometallic reaction pathways from supplementary information at scale.

## Overview

This plugin processes XYZ structures extracted by GoldDIGR through a five-stage computational pipeline:

| Stage | Directory | Description |
|-------|-----------|-------------|
| 1 | `01-charge-scan-tsopt-irc/` | Charge sampling (−1, 0, +1) → xTB frequency → TSOPT → IRC |
| 2 | `02-spin-wbo-scan/` | Two-pass spin-polarized energy scan + Wiberg bond order extraction |
| 3 | `03-irc-analysis/` | YARP bond-electron matrix analysis of IRC endpoints and trajectory |
| 4 | `04-sankey/` | Electron-transport Sankey diagrams from BEM time series |
| 5 | `05-reaction-classification/` | Rule-based classification into canonical organometallic reaction classes |

## Pipeline Flow

```
XYZ structures (from GoldDIGR SI extraction)
    │
    ▼
[01] Charge scan: charges {-1, 0, +1} × auto multiplicity
     → xTB frequency (gate: |ν_imag| > 20 cm⁻¹)
     → Pysisyphus TSOPT (RS-P-RFO)
     → Pysisyphus IRC (Euler predictor–corrector)
    │
    ▼
[02] Spin/WBO scan along IRC trajectory
     Pass 1: tblite spin-state screening (UHF 0–6)
     Pass 2: standard xTB property extraction at best spin state
    │
    ▼
[03] IRC analysis: YARP BEM construction for each IRC frame
     → Bond-electron matrices, adjacency lists, geometry graphs
    │
    ▼
[04] Sankey diagram generation
     → Transport + storage electron-flow diagrams
    │
    ▼
[05] Reaction classification
     → OA, RE, MI, β-atom elimination, C–H activation, transmetalation
```

## Environment

```bash
conda env create -f environment.yml
conda activate another-yarp
```

Key software: xTB 6.5.1, tblite 0.3.0, pysisyphus 0.7.6, YARP 0.0.1, molSimplify 1.8.0, RDKit 2025.3.3.

## HPC Requirements

- SLURM scheduler (templates in `01-charge-scan-tsopt-irc/templates/`)
- Modules: `conda/2025.02`, `intel-mkl`, `openmpi`

## Citation

> Zhao Li, Ivan Yiwen Yang, Brett M. Savoie. "Reaction Database for Catalysis and Organometallics via Freely Available Supplementary Information." *ArXiv* (2025).
