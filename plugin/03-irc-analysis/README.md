# Stage 3: IRC Analysis (YARP Bond-Electron Matrices)

## Overview

Runs YARP on IRC endpoint structures and trajectory frames to generate bond-electron
matrices (BEMs), adjacency matrices, and compressed adjacency-list representations.

## Usage

```bash
python yarp_results_builder.py /path/to/job_dir/ --charge 0
```

The job directory should contain `finished_first.xyz`, `finished_last.xyz`,
`finished_irc.trj`, and optionally `ts_final_geometry.xyz`.

## Outputs

Creates `IRC_Analysis/` with per-source subdirectories:

```
IRC_Analysis/
├── finished_first/
│   ├── frame_00000.json          # BEM + adjacency + metadata
│   ├── frame_00000_adjacency.csv
│   └── frame_00000_bond_electrons.csv
├── finished_last/
│   └── ...
├── finished_irc/
│   ├── frame_00000.json
│   ├── ...
│   └── ts_frame.txt              # TS frame index (max energy)
└── ts_final_geometry/
    └── ...
```

## Dependencies

- `yarp` (YARP) — bond-electron matrix generation
- `metal_ligand/yarp_helpers.py` — silent wrappers for YARP yarpecule objects
