````markdown
# no-three-in-line---CP-SAT

CP-SAT (OR-Tools) code for finding maximal configurations for the **no-three-in-line** problem on an \(n\times n\) grid.

> **No-three-in-line:** place as many grid points as possible so that no three are collinear (lines of any slope).  
> Under the common “2-per-row” formulation, the target size is \(2n\).

<!-- After Zenodo setup, add a DOI badge here.
Example (replace with your Zenodo badge snippet):
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
-->

## What’s in this repo

- `no_three_in_line.py` — main script:
  - **Direct model:** variables on all grid points with collinearity constraints.
  - **Symmetry-reduced model (`--sym`):** reduces variables/constraints via rotational symmetry (as described in the companion paper).
  - Optional **verification** (`--verify`) to check the no-three-in-line property by enumerating (deduped) line intersection sets.
- `Table 1.txt` — table/data included to accompany the companion paper (verbatim export).

## Requirements

- Python 3.9+ recommended
- OR-Tools (CP-SAT)

Install dependency:

```bash
python -m pip install ortools
````

## Quick start

Run the solver on an `n×n` grid:

```bash
python no_three_in_line.py --n 30 --time_limit 600 --seed 1 --workers 8
```

Use the symmetry-reduced model:

```bash
python no_three_in_line.py --n 30 --sym --time_limit 600 --seed 1 --workers 8
```

Verify the resulting configuration (slower, but useful for sanity checks):

```bash
python no_three_in_line.py --n 30 --sym --time_limit 600 --seed 1 --workers 8 --verify
```

### Full CLI

```text
usage: no_three_in_line.py [-h] --n N [--sym] [--time_limit TIME_LIMIT]
                           [--seed SEED] [--workers WORKERS] [--log_search]
                           [--no_dedupe_line_orbits] [--no_dedupe_incidence]
                           [--verify]
```

Notes:

* `--workers` controls parallelism inside OR-Tools CP-SAT.
* The `--no_dedupe_*` flags are primarily for debugging / benchmarking orbit and incidence deduplication.

## Reproducibility

For archival, reproducible citation, use a tagged GitHub release that has been archived to Zenodo (DOI per release).
See `CITATION.cff` for the preferred citation metadata.

## License

MIT (see `LICENSE`).

## Acknowledgements

This codebase was generated from the companion paper using an LLM and then validated and edited for correctness.
(See the paper for the formal model definitions and discussion.)

```
::contentReference[oaicite:0]{index=0}
```
