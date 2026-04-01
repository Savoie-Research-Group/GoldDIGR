# Stage 5: Reaction Classification

## Overview

Classifies IRC-validated reactions into canonical organometallic reaction classes
based on bond-forming/breaking patterns at metal centers.

## Usage

```bash
# Classify from tar.gz workflow archives
python batch_classify.py /path/to/targz_dir/ -o results.jsonl

# Classify specific archives
python batch_classify.py file1.tar.gz file2.tar.gz -o results.jsonl
```

## Supported Reaction Classes

| Class | Checker | Description |
|-------|---------|-------------|
| OA | `FailOA_only_checker.py` | Oxidative addition |
| RE | `FailRE_only_checker.py` | Reductive elimination |
| MI | `FailMI_only_checker.py` | Migratory insertion |
| BXA | `FailBXA_only_checker.py` | β-Atom elimination |
| BXE | `FailBXE_only_checker.py` | β-Heteroatom elimination |
| CH | `FailCH_only_checker.py` | C–H activation |
| MET | `FailMET_only_checker.py` | σ-Bond metathesis |
| TM | `FailTM_only_checker.py` | Transmetalation |

## Architecture

- `batch_classify.py` — Entry point; extracts BEL graphs from tar.gz archives,
  deduplicates by graph hash, runs all checkers, writes JSONL output.
- `rxnclass_helper.py` — Shared utilities for bond-change parsing, metal detection,
  electron counting, and BEL graph manipulation.
- `bel_parser.py` — Parser for bond-electron line (BEL) notation.
- `metal_ligand/` — Package providing atom labeling, matrix operations, electron
  counting, and coordination geometry utilities.

## Output Format

One JSON line per unique reaction in `results.jsonl`:

```json
{
  "hash": "abc123...",
  "source_tar": "path/to/file.tar.gz",
  "reactant_bel": "...",
  "product_bel": "...",
  "bond_changes_sigma": "C13-O43:1>n;O42-H44:1>n;O43-H44:n>1",
  "positive_classes": ["reductive_elimination"],
  "classifications": { ... }
}
```
