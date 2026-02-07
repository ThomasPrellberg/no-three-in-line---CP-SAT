# no-three-in-line---CP-SAT
CP-SAT code for no-three-in-line

This repository contains a python script to find maximal solutions for the no-three-in-line-problem on an nxn grid.

It has been created from the companion paper using ChatGPT and subsequently validated.

usage: no_three_in_line.py [-h] --n N [--sym] [--time_limit TIME_LIMIT] [--seed SEED] [--workers WORKERS] [--log_search] [--no_dedupe_line_orbits] [--no_dedupe_incidence] [--verify]

options:

  -h, --help            show this help message and exit
  
  --n N
  
  --sym                 Use symmetry-reduced model (Eq. CPSATred). Default: direct model (Eq. CPSAT).
  
  --time_limit TIME_LIMIT
  
  --seed SEED
  
  --workers WORKERS
  
  --log_search
  
  --no_dedupe_line_orbits
  
                        Sym model: do not reduce lines by Gamma-orbits.
                        
  --no_dedupe_incidence
  
                        Sym model: do not dedupe constraints by incidence vectors.
                        
  --verify              Verify no-three-in-line by checking all (deduped) line intersection sets.
