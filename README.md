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
