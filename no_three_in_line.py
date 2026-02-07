#!/usr/bin/env python3
"""
No-three-in-line (2 per row) CP-SAT models:

(1) Direct model  F_n  (Eq. (CPSAT) in the LaTeX): variables x_{ij} on G_n
    - for every line with >=3 grid points: sum_{(i,j) in line} x_{ij} <= 2
    - for every row j: sum_i x_{ij} == 2

(2) Symmetry-reduced model F_n^sym (Eq. (CPSATred) in the LaTeX): variables y_{ij} on H_n
    - rotational symmetry via rho(i,j)=(j,n-1-i), with modified rule for odd n:
      anti-diagonal empty; diagonal reps expand under rho^2 only.
    - line reps: one line per Gamma-orbit, Gamma=<rho> (even) or <rho^2> (odd)
    - constraints rewritten using coefficients c_ell(i,j) and d_r(i,j)

Requires: ortools (CP-SAT).
"""

from __future__ import annotations

from dataclasses import dataclass
from math import gcd
from typing import Dict, Iterable, List, Optional, Tuple, Set
import argparse
import time

from ortools.sat.python import cp_model


# -----------------------------
# Basic geometry / encoding
# -----------------------------

def p_of_ij(i: int, j: int, n: int) -> int:
    """Flatten (i,j) to integer id."""
    return i * n + j

def ij_of_p(p: int, n: int) -> Tuple[int, int]:
    """Unflatten integer id to (i,j)."""
    return (p // n, p % n)

def in_grid(i: int, j: int, n: int) -> bool:
    return 0 <= i < n and 0 <= j < n

def rho_90(i: int, j: int, n: int) -> Tuple[int, int]:
    """rho(i,j) = (j, n-1-i)."""
    return (j, n - 1 - i)

def rotate_point(p: int, n: int, k: int) -> int:
    """Apply rho^k to point p (k mod 4)."""
    k %= 4
    i, j = ij_of_p(p, n)
    for _ in range(k):
        i, j = rho_90(i, j, n)
    return p_of_ij(i, j, n)

def rotate_line_key(line_key: Tuple[int, ...], n: int, k: int) -> Tuple[int, ...]:
    """Rotate a line (given by sorted tuple of points) and return sorted tuple."""
    pts = [rotate_point(p, n, k) for p in line_key]
    pts.sort()
    return tuple(pts)


# -----------------------------
# Lines: enumerate and deduplicate by intersection sets
# -----------------------------

def primitive_directions(n: int) -> Iterable[Tuple[int, int]]:
    """
    Primitive direction vectors (dx,dy) with gcd(dx,|dy|)=1 and a canonical orientation:
      dx > 0, or dx == 0 and dy > 0
    dy ranges negative to allow negative slopes.
    """
    for dx in range(0, n):
        for dy in range(-(n - 1), n):
            if dx == 0 and dy == 0:
                continue
            # canonical orientation to avoid sign duplicates
            if dx < 0 or (dx == 0 and dy <= 0):
                continue
            g = gcd(dx, abs(dy))
            if g != 1:
                continue
            yield dx, dy

def generate_line_intersection_sets(n: int, min_points: int = 3) -> Dict[Tuple[int, ...], Tuple[int, ...]]:
    """
    Return dict: key -> points, where key is canonical encoding of the intersection set (sorted tuple of point ids).
    Values are also stored as sorted tuples (same as keys).
    This matches the LaTeX description: "dictionary whose keys are canonical encodings of distinct intersection sets".
    """
    lines: Dict[Tuple[int, ...], Tuple[int, ...]] = {}

    for dx, dy in primitive_directions(n):
        for i0 in range(n):
            for j0 in range(n):
                # start point is minimal in this direction if predecessor is out of bounds
                if in_grid(i0 - dx, j0 - dy, n):
                    continue

                pts: List[int] = []
                i, j = i0, j0
                while in_grid(i, j, n):
                    pts.append(p_of_ij(i, j, n))
                    i += dx
                    j += dy

                if len(pts) >= min_points:
                    pts.sort()
                    key = tuple(pts)
                    # Deduplicate by intersection set key (should already be unique with our construction,
                    # but we keep the dict to follow the paper and be robust).
                    if key not in lines:
                        lines[key] = key

    return lines


# -----------------------------
# Symmetry: H_n, orbits O(i,j), Gamma-line orbit reps
# -----------------------------

def fundamental_domain_H(n: int) -> List[Tuple[int, int]]:
    """
    H_n as defined in the LaTeX:
      even n: [n/2]_0 x [n/2]_0
      odd  n: [(n+1)/2]_0 x [(n-1)/2]_0
    """
    if n % 2 == 0:
        a = n // 2
        return [(i, j) for i in range(a) for j in range(a)]
    else:
        a = (n + 1) // 2
        b = (n - 1) // 2
        return [(i, j) for i in range(a) for j in range(b)]

def orbit_O(rep: Tuple[int, int], n: int) -> List[Tuple[int, int]]:
    """
    Orbit O(i,j) as in the LaTeX:
      - if n even OR i != j: full <rho>-orbit size 4
      - if n odd and i == j: <rho^2>-orbit size 2: {(i,i),(n-1-i,n-1-i)}
    (Anti-diagonal handling for odd n is enforced separately by leaving those points unmapped / fixed to 0.)
    """
    i, j = rep
    if n % 2 == 1 and i == j:
        return [(i, i), (n - 1 - i, n - 1 - i)]
    # full 4-cycle
    pts = []
    ii, jj = i, j
    for _ in range(4):
        pts.append((ii, jj))
        ii, jj = rho_90(ii, jj, n)
    # distinct by construction for our use-cases; still, be safe:
    out = []
    seen = set()
    for t in pts:
        if t not in seen:
            seen.add(t)
            out.append(t)
    return out

def anti_diagonal_points(n: int) -> Set[int]:
    """A_n = {(i,n-1-i)} as flattened ids."""
    return {p_of_ij(i, n - 1 - i, n) for i in range(n)}

def gamma_rotations(n: int) -> List[int]:
    """Gamma_n: <rho> (even) -> [0,1,2,3]; <rho^2> (odd) -> [0,2]."""
    return [0, 1, 2, 3] if n % 2 == 0 else [0, 2]

def canonical_under_gamma(line_key: Tuple[int, ...], n: int) -> Tuple[int, ...]:
    """Return canonical representative (lexicographically minimal rotated key) under Gamma_n."""
    rots = gamma_rotations(n)
    keys = [rotate_line_key(line_key, n, k) for k in rots]
    return min(keys)

def line_representatives_under_gamma(
    lines: Dict[Tuple[int, ...], Tuple[int, ...]],
    n: int
) -> List[Tuple[int, ...]]:
    """
    Keep one representative per Gamma_n-orbit of lines: choose the canonical (min) one.
    Returns list of canonical line keys.
    """
    reps: List[Tuple[int, ...]] = []
    seen: Set[Tuple[int, ...]] = set()
    for key in lines.keys():
        canon = canonical_under_gamma(key, n)
        if canon in seen:
            continue
        seen.add(canon)
        reps.append(canon)
    return reps


# -----------------------------
# Model builders
# -----------------------------

@dataclass
class DirectModel:
    model: cp_model.CpModel
    x: Dict[int, cp_model.IntVar]  # point id -> BoolVar
    lines: List[Tuple[int, ...]]   # line intersection sets used

@dataclass
class SymModel:
    model: cp_model.CpModel
    y: Dict[Tuple[int, int], cp_model.IntVar]  # rep in H_n -> BoolVar
    H: List[Tuple[int, int]]
    lines_rep: List[Tuple[int, ...]]           # Gamma representatives used (post-incidence dedupe optionally)
    point_to_rep: Dict[int, Tuple[int, int]]   # point id -> representative in H_n (anti-diagonal points absent for odd n)

def build_direct_model(n: int, dedupe_lines: bool = True) -> DirectModel:
    """
    Build the direct model (Eq. CPSAT): x_{ij} on G_n, all line constraints, and 2-per-row constraints.
    """
    model = cp_model.CpModel()
    x: Dict[int, cp_model.IntVar] = {}

    # Variables x_{ij}
    for i in range(n):
        for j in range(n):
            p = p_of_ij(i, j, n)
            x[p] = model.NewBoolVar(f"x_{i}_{j}")

    # Row constraints: sum_i x_{ij} = 2 for each row j
    for j in range(n):
        model.Add(sum(x[p_of_ij(i, j, n)] for i in range(n)) == 2)

    # Line constraints: sum_{p in line} x_p <= 2
    lines_dict = generate_line_intersection_sets(n, min_points=3)
    line_keys = list(lines_dict.keys()) if dedupe_lines else [v for v in lines_dict.values()]

    for key in line_keys:
        model.Add(sum(x[p] for p in key) <= 2)

    return DirectModel(model=model, x=x, lines=line_keys)

def build_symmetry_reduced_model(
    n: int,
    *,
    dedupe_line_orbits: bool = True,
    dedupe_incidence_vectors: bool = True
) -> SymModel:
    """
    Build the symmetry reduced model (Eq. CPSATred) with y_{ij} on H_n.

    - For odd n: anti-diagonal points are fixed to 0 by NOT mapping them to any y-variable.
    - Line family: one rep per Gamma-orbit of lines (optional).
    - Optionally deduplicate further by incidence vectors on H_n.
    """
    model = cp_model.CpModel()
    H = fundamental_domain_H(n)
    y: Dict[Tuple[int, int], cp_model.IntVar] = {}
    rep_to_varid: Dict[Tuple[int, int], int] = {}

    for idx, rep in enumerate(H):
        i, j = rep
        y[rep] = model.NewBoolVar(f"y_{i}_{j}")
        rep_to_varid[rep] = idx

    # Build point_to_rep mapping via orbits O(i,j)
    point_to_rep: Dict[int, Tuple[int, int]] = {}
    if n % 2 == 1:
        anti = anti_diagonal_points(n)
    else:
        anti = set()

    # Row contribution bookkeeping: for each row r, collect (rep, coeff)
    row_terms: List[List[Tuple[Tuple[int, int], int]]] = [[] for _ in range(n)]
    # We'll compute per rep a small dict row->coeff
    for rep in H:
        counts_by_row: Dict[int, int] = {}
        for (ii, jj) in orbit_O(rep, n):
            p = p_of_ij(ii, jj, n)
            if p in anti:
                # odd n: anti-diagonal is empty (fixed to 0), so skip mapping
                continue
            # map this point to its representative orbit variable
            point_to_rep[p] = rep
            counts_by_row[jj] = counts_by_row.get(jj, 0) + 1
        for r, c in counts_by_row.items():
            row_terms[r].append((rep, c))

    # Row constraints: sum_{(i,j) in H_n} d_r(i,j) y_{ij} = 2
    for r in range(n):
        model.Add(sum(c * y[rep] for (rep, c) in row_terms[r]) == 2)

    # Lines in the full grid
    lines_dict = generate_line_intersection_sets(n, min_points=3)
    all_line_keys = list(lines_dict.keys())

    # Restrict to one representative per Gamma-orbit of lines, if requested
    if dedupe_line_orbits:
        line_keys = line_representatives_under_gamma(lines_dict, n)
    else:
        line_keys = all_line_keys

    # Add line constraints in reduced form using incidence counts per orbit variable
    seen_incidence: Set[Tuple[Tuple[int, int], ...]] = set()
    final_line_keys: List[Tuple[int, ...]] = []

    for key in line_keys:
        coeffs: Dict[Tuple[int, int], int] = {}
        for p in key:
            rep = point_to_rep.get(p)
            if rep is None:
                # (odd n) anti-diagonal points are absent => x=0, so contribute nothing
                continue
            coeffs[rep] = coeffs.get(rep, 0) + 1

        if not coeffs:
            # Entire line lies on anti-diagonal (possible only for odd n); constraint is 0<=2 so skip.
            continue

        if sum(coeffs.values()) <= 2:
            # Tautological: with y ∈ {0,1}, LHS ≤ sum(coeffs) ≤ 2
            continue

        if dedupe_incidence_vectors:
            # Canonical key for this inequality: sorted (varid, coeff)
            inc = tuple(sorted((rep_to_varid[rep], c) for rep, c in coeffs.items()))
            if inc in seen_incidence:
                continue
            seen_incidence.add(inc)

        model.Add(sum(c * y[rep] for rep, c in coeffs.items()) <= 2)
        final_line_keys.append(key)

    return SymModel(model=model, y=y, H=H, lines_rep=final_line_keys, point_to_rep=point_to_rep)


# -----------------------------
# Solve + extract solutions
# -----------------------------

def solve_model(
    model: cp_model.CpModel,
    *,
    time_limit: Optional[float] = None,
    seed: Optional[int] = None,
    workers: Optional[int] = None,
    log_search: bool = False
) -> Tuple[cp_model.CpSolver, cp_model.CpSolverStatus]:
    solver = cp_model.CpSolver()
    if time_limit is not None:
        solver.parameters.max_time_in_seconds = float(time_limit)
    if seed is not None:
        solver.parameters.random_seed = int(seed)
    if workers is not None:
        solver.parameters.num_search_workers = int(workers)
    solver.parameters.log_search_progress = bool(log_search)
    status = solver.Solve(model)
    return solver, status

def extract_direct_solution(n: int, dm: DirectModel, solver: cp_model.CpSolver) -> List[Tuple[int, int]]:
    occ = []
    for p, var in dm.x.items():
        if solver.Value(var):
            occ.append(ij_of_p(p, n))
    occ.sort()
    return occ

def extract_sym_solution(n: int, sm: SymModel, solver: cp_model.CpSolver) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]]]:
    """
    Return (occupied_points_in_Gn, selected_reps_in_Hn).
    """
    reps = [rep for rep, var in sm.y.items() if solver.Value(var)]
    reps.sort()

    occ_pts: List[Tuple[int, int]] = []
    anti = anti_diagonal_points(n) if n % 2 == 1 else set()

    for rep in reps:
        for (ii, jj) in orbit_O(rep, n):
            p = p_of_ij(ii, jj, n)
            if p in anti:
                continue
            occ_pts.append((ii, jj))

    occ_pts.sort()
    return occ_pts, reps

def verify_no_three_in_line(n: int, occupied: Set[int], lines: Iterable[Tuple[int, ...]]) -> bool:
    """
    Verify: for every line intersection set, at most 2 occupied points.
    (This is an expensive check if you pass all lines; intended for sanity checks.)
    """
    for key in lines:
        cnt = 0
        for p in key:
            if p in occupied:
                cnt += 1
                if cnt > 2:
                    return False
    return True


# -----------------------------
# CLI
# -----------------------------

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--sym", action="store_true", help="Use symmetry-reduced model (Eq. CPSATred). Default: direct model (Eq. CPSAT).")
    ap.add_argument("--time_limit", type=float, default=None)
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--log_search", action="store_true")
    ap.add_argument("--no_dedupe_line_orbits", action="store_true", help="Sym model: do not reduce lines by Gamma-orbits.")
    ap.add_argument("--no_dedupe_incidence", action="store_true", help="Sym model: do not dedupe constraints by incidence vectors.")
    ap.add_argument("--verify", action="store_true", help="Verify no-three-in-line by checking all (deduped) line intersection sets.")
    args = ap.parse_args()

    n = args.n
    t0 = time.time()

    if not args.sym:
        dm = build_direct_model(n, dedupe_lines=True)
        build_time = time.time() - t0
        print(f"[direct] n={n} vars={len(dm.x)} line_constraints={len(dm.lines)} build_time={build_time:.3f}s")

        solver, status = solve_model(dm.model, time_limit=args.time_limit, seed=args.seed, workers=args.workers, log_search=args.log_search)
        print(f"[direct] status={solver.StatusName(status)} walltime={solver.WallTime():.3f}s")

        if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            occ = extract_direct_solution(n, dm, solver)
            print(f"[direct] found |occ|={len(occ)} (expected 2n={2*n})")
            print(occ)

            if args.verify:
                occ_set = {p_of_ij(i, j, n) for (i, j) in occ}
                ok = verify_no_three_in_line(n, occ_set, dm.lines)
                print(f"[direct] verify_no_three_in_line={ok}")

    else:
        sm = build_symmetry_reduced_model(
            n,
            dedupe_line_orbits=not args.no_dedupe_line_orbits,
            dedupe_incidence_vectors=not args.no_dedupe_incidence
        )
        build_time = time.time() - t0
        print(f"[sym] n={n} |H|={len(sm.H)} vars={len(sm.y)} line_constraints={len(sm.lines_rep)} build_time={build_time:.3f}s")

        solver, status = solve_model(sm.model, time_limit=args.time_limit, seed=args.seed, workers=args.workers, log_search=args.log_search)
        print(f"[sym] status={solver.StatusName(status)} walltime={solver.WallTime():.3f}s")

        if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            occ_pts, reps = extract_sym_solution(n, sm, solver)
            print(f"[sym] found |occ|={len(occ_pts)} (expected 2n={2*n})  |reps|={len(reps)}")
            print("[sym] reps in H_n:")
            print(reps)
            print("[sym] occupied points in G_n:")
            print(occ_pts)

            if args.verify:
                # verify against the FULL set of line intersection sets (not just Gamma reps),
                # since the property is Gamma-invariant this is mostly a sanity check.
                all_lines = generate_line_intersection_sets(n, min_points=3).keys()
                occ_set = {p_of_ij(i, j, n) for (i, j) in occ_pts}
                ok = verify_no_three_in_line(n, occ_set, all_lines)
                print(f"[sym] verify_no_three_in_line={ok}")


if __name__ == "__main__":
    main()
