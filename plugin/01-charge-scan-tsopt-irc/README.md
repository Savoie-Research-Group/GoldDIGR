# Stage 1: Charge Scan → TSOPT → IRC

## Overview

Samples net charges {−1, 0, +1} for each extracted XYZ structure, determines the lowest spin
multiplicity for each charge, then runs xTB frequency analysis → pysisyphus TSOPT → pysisyphus IRC.

## Usage

```bash
# Edit AAAS_folders.txt with one XYZ directory path per line
bash 5_master_workflow.sh
```

## Key Files

| File | Purpose |
|------|---------|
| `5_master_workflow.sh` | Master orchestrator: loops over XYZ files × charges, submits SLURM jobs |
| `calculate_mult.py` | Determines lowest spin multiplicity from electron count and charge |
| `6_restart_untouched.sh` | Restarts jobs that were queued but never ran |
| `check_job_status.sh` | Checks completion status of submitted jobs |
| `cleanup.sh` | Archives completed job directories into zip files |
| `zip_crashed.sh` | Archives crashed job directories separately |
| `glue/runner.sh` | Per-job runner: executes the 6-stage pipeline sequentially |
| `glue/submit_all.sh` | Batch SLURM submission for all (XYZ, charge) combinations |
| `plugin.yaml` | Machine-readable workflow manifest |
| `templates/` | Pysisyphus YAML configs and SLURM job templates |

## Charge/Multiplicity Logic

- For each XYZ, three jobs are submitted at charges −1, 0, and +1.
- `calculate_mult.py` sums atomic numbers, subtracts the charge, and returns
  multiplicity 1 (singlet) for even electron count or 2 (doublet) for odd.
- Output directory naming: `{xyz_stem}_{charge}_{multiplicity}/` (e.g., `19_-1_2/`).

## Computational Settings

- **Frequency**: GFN2-xTB Hessian; gate = |ν_imag| > 20 cm⁻¹
- **TSOPT**: RS-P-RFO via pysisyphus; Hessian recalc every step; Gaussian convergence; max 50 cycles
- **IRC**: Euler predictor–corrector; forward + backward; endpoint optimization max 200 cycles
