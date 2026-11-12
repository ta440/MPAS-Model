# Kinematic transport tests for MPAS-Atmosphere

These tests run the production MPAS-A scalar transport with the dynamics switched off
and the mass fluxes prescribed from an analytic flow, so that results can be compared
with known exact solutions. Everything else (mesh, decomposition, halos, I/O, timers)
is the normal model.

Location in the fork: `testing_and_setup/atmosphere/kinematic_tests/`.

## Cases

| `config_kinematic_case` | Flow | Exact solution | Status |
|---|---|---|---|
| `lauritzen_nondiv` | Lauritzen et al. (2012, GMD 5, 887–901), eqs. 18–20: non-divergent deformation with a background solid-body rotation, T = 12 days | the initial state at t = T | available |
| `lauritzen_div` | Lauritzen et al. (2012), eqs. 21–22: divergent deformation with the same background rotation, T = 12 days | the initial state at T, and density 1 at T | available |
| `dcmip2012_1_2` | DCMIP 2012 test 1-2 (test case document v1.7, eqs. 36–38; Kent et al. 2014, QJRMS 140, 1279–1293): solid-body zonal flow plus five Hadley-like cells in latitude–height that reverse at T/2, T = 1 day, ztop = 12 km | the initial state at T, and q1 at T/2 from backward trajectories (`dcmip2012_1_2_exact_T_half.npz`) | available |
| `dcmip2012_1_1` | DCMIP 2012 test 1-1: 3D deformation | initial state at T | planned |
| `terminator` | Lauritzen et al. (2015) terminator toy chemistry | conservation of Cl + 2 Cl2 | planned |

Each case has a directory with namelist and streams templates. The tracers are the
passive-tracer array (`config_num_passive_tracers`).

- Lauritzen cases: tracer `n` has shape `mod(n-1,5)`: 0 uniform (q = 1), 1 Gaussian
  hills, 2 cosine bells, 3 slotted cylinders, 4 correlated cosine bells (ξ = −0.8 χ² + 0.9,
  with χ the cosine bells; Lauritzen et al. 2012, eqs. 13–15). Shapes 2 and 4 form the
  tracer pair for the mixing diagnostics.
- `dcmip2012_1_2`: tracer `n` has shape `mod(n-1,3)`: 0 uniform, 1 q1 (the layer between
  2 and 5 km, DCMIP eq. 43), 2 ξ = −0.8 q1² + 0.9. Shapes 1 and 2 form the mixing pair.

Tracers `n` and `n + (number of shapes)` are identical copies, which must stay bitwise
identical.

### How the fluxes are prescribed

- `lauritzen_nondiv`: the streamfunction is averaged over the step (3-point
  Gauss–Legendre in time) at the cell vertices, and the edge flux is its difference along
  the edge. This is the exact edge- and step-mean flux, and its discrete divergence is
  zero, so ρ_zz stays 1.
- `lauritzen_div`: the normal velocity is averaged over each edge and over the step with
  3 × 3 Gauss–Legendre points. The density starts at 1 and is advanced on the mesh with
  the stages of the MPAS RK3, as the dynamics would; `ruAvg` is the final-stage mass flux,
  and ρ_zz at the new time level is updated from it with the transport's own arithmetic,
  so a uniform tracer stays uniform. The halo of ρ_zz is exchanged after each stage.
  The exact density at T is 1, so its error at T is reported as well.
- `dcmip2012_1_2`: the zonal flux is ρ_zz times the vertex-streamfunction flux of
  u0 cos φ. The meridional flux factors as (time) × (latitude) × (height): the step mean
  and the layer mean are exact, and the edge mean in latitude uses 3-point Gauss–Legendre.
  `wwAvg` is then derived from the discrete horizontal divergence, so the 3D mass
  divergence of every cell is zero to roundoff (the run stops if it exceeds 1e-12
  relative) and ρ_zz keeps its initial isothermal profile. This `wwAvg` differs from the
  point-sampled analytic ρw at second order in the mesh spacing.

## Requirements

- MPAS-Model built from this branch with `PRECISION=double`, for both cores:
  ```
  make gfortran CORE=init_atmosphere PRECISION=double
  make clean CORE=atmosphere
  make gfortran CORE=atmosphere PRECISION=double
  ```
  This leaves `init_atmosphere_model` and `atmosphere_model` in the top directory.
  Clean with `CORE=atmosphere` (the core about to be built): `make clean
  CORE=init_atmosphere` would delete the `init_atmosphere_model` just built.
- MPI (`mpirun`, or set `MPAS_MPIRUN`, e.g. `MPAS_MPIRUN="srun -n"`).
- Quasi-uniform meshes from the MPAS mesh page
  (https://mpas-dev.github.io/atmosphere/atmosphere_meshes.html). Each tarball holds
  `x1.N.grid.nc`, `x1.N.graph.info` and partition files for common rank counts:

  | Mesh | Spacing | Default `dt` (s), Lauritzen | Default `dt` (s) and levels, `dcmip2012_1_2` |
  |---|---|---|---|
  | `x1.10242` | 240 km | 900 | 300, L30 |
  | `x1.40962` | 120 km | 450 | 150, L60 |
  | `x1.163842` | 60 km | 225 | 75, L120 |

  Unpack each into its own directory, for example `meshes/x1.10242/`. If a partition
  file for your rank count is missing, the run script makes one with `gpmetis` when it
  is on the path.
- Python 3 with `numpy` and `netCDF4` for the checker. The scripts call `python3`; to use
  another interpreter (for example a virtual environment), set `PYTHON=/path/to/python`.

## Running one test

```
./run_kinematic_test.sh -c lauritzen_nondiv -m x1.10242 -d /path/to/meshes/x1.10242
```

| Option | Meaning | Default |
|---|---|---|
| `-c` | case | (required) |
| `-m` | mesh name | (required) |
| `-d` | directory holding the mesh files | (required) |
| `-b` | directory holding the two executables | top of this repository |
| `-n` | MPI ranks | 4 |
| `-t` | time step (s) | by case and mesh, table above |
| `-k` | Runge–Kutta order, `config_time_integration_order` (2 or 3) | 2 |
| `-L` | FCT limiter, `config_monotonic` (`on` or `off`) | `on` |
| `-l` | vertical levels | Lauritzen: 5; `dcmip2012_1_2`: by mesh, table above |
| `-N` | passive tracers | Lauritzen: 5; `dcmip2012_1_2`: 3 |
| `-r` | run duration | Lauritzen: `12_00:00:00`; `dcmip2012_1_2`: `1_00:00:00` |
| `-i` | history output interval | Lauritzen: `3_00:00:00`; `dcmip2012_1_2`: `0_01:00:00` (DCMIP Table VII) |
| `-C` | tracers per chunk in the FCT stage, `config_fct_tracer_chunk`; 0 uses the original routine | not set (original routine) |
| `-o` | run directory | `./runs/<case>_<mesh>_N<ntracers>_rk<order>_lim<on/off>[_C<n>]` |

The script writes the namelists and streams from the case templates, runs
`init_atmosphere_model` (init case 14) and then `atmosphere_model`, and runs the
checker. The exit code is the checker's: 0 if all checks pass.

## Suites

```
./run_suite.sh quick        -d /path/to/meshes      # all three cases, 240 km
./run_suite.sh convergence  -d /path/to/meshes      # all three cases, 240, 120, 60 km
./run_suite.sh scaling      -d /path/to/meshes      # lauritzen_nondiv, 120 km, N = 5, 10, 20, 50, 100, 200, 1 day
./run_suite.sh ranks        -d /path/to/meshes      # all three cases, 240 km, 1 and 4 ranks: bitwise identical
```

`-d` here is the parent directory holding one subdirectory per mesh; `-b`, `-n`, `-k`
and `-L` are passed on to every run. The scaling suite uses whole blocks of the five
shapes and covers the range that matters operationally: about 11 tracers in a default
MPAS run (N = 10 here), and up to 200 for chemistry. The `dcmip2012_1_2` suite runs write
12-hourly output (`-i 0_12:00:00`): hourly files of u, ρ_zz and 3 tracers would take about
0.45 GB per run at 240 km L30, 3.5 GB at 120 km L60 and 27 GB at 60 km L120.

## Time integration

The transport uses the three-stage Runge–Kutta scheme of Wicker and Skamarock (2002),
with the stages set by `config_time_integration_order`:

| Order | Stages | Accuracy |
|---|---|---|
| 2 (MPAS default) | Δt/2, Δt/2, Δt | second order; larger stable step |
| 3 | Δt/3, Δt/2, Δt | third order for linear problems |

See Section 3.3 of the MPAS-A v8 description (Skamarock et al.). The tests use 2 unless
`-k 3` is given, so that they measure the transport as MPAS runs it. Skamarock and
Gassmann (2011) used order 3, so comparisons with their published errors need `-k 3`.
The template always writes the order into the namelist, so a run never depends on the
default.

## What the checker tests

`check_kinematic_test.py <run_dir>` reads the history files and the namelist.

| Check | Pass condition |
|---|---|
| Uniform tracer | max \|q − 1\| < 1e-12 at every output time |
| Mass conservation | relative change in Σ A ρ_zz q Δζ < 1e-12 for every tracer (pairwise summation; checker noise about 1e-16) |
| Shape preservation (if `config_monotonic`) | no value outside the initial global range by more than 1e-12 of the range |
| Identical copies | tracer n bitwise equal to tracer n − (number of shapes) at every output time |
| Identical levels (2D cases) | every level bitwise equal to level 1 |
| Density constant (`lauritzen_nondiv`, `dcmip2012_1_2`) | ρ_zz bitwise equal to its t = 0 value (1 for `lauritzen_nondiv`, where the t = 0 file holds the init density) |
| Exact-solution levels (`dcmip2012_1_2`) | the run's layer midpoints match a level set (L30, L60, L120) of the check file |
| Regression (if a reference exists, and the run reaches T) | l2 errors at T (and at T/2 against the exact solution), and the total mixing at T/2, within 1% of `references.json` |

It also prints, without pass/fail:

- the normalised l1, l2 and l∞ errors at T, volume-weighted over all levels (Lauritzen
  et al. 2012, eqs. 26–28; DCMIP eq. 14). For the 2D cases these equal the level-1 norms;
- for `dcmip2012_1_2`, the same norms of q1 and ξ at T/2 against the exact solution;
- for `lauritzen_div`, the maximum and rms error of the density at T (exact value 1);
- the mixing diagnostics ℓr ("real" mixing), ℓu (range-preserving unmixing) and ℓo
  (overshooting) for the correlated pair (Lauritzen et al. 2012, Sect. 3.5 and Appendix C):
  for `lauritzen_nondiv` on level 1 at T/2; for `dcmip2012_1_2` at T/2 and T, for the
  levels with midpoints in 4400–5400 m (headline) and by level, with χ_min = 0. The band
  data at T/2 are saved to `mixing_T_half.npz` for correlation plots;
- the transport cost from the MPAS timers;

and writes everything to `results.json` in the run directory.

`check_kinematic_test.py --compare <run_a> <run_b>` checks that two runs have bitwise
identical history files; the `ranks` suite uses it.

The regression references are keyed by case, mesh, time step, RK order and limiter
setting; runs with no matching reference skip that check. They come from our runs with
gfortran on Apple M3.
Other compilers and machines change the results at roundoff level, well inside the 1%
tolerance; a larger change means the transport or the test setup has changed. They
are not accuracy targets.

## Notes

- Init case 14 writes a flat, isothermal state at rest. It is only for runs with
  `config_kinematic_case` set; with it `'off'` the model would run the full dynamics
  from an unbalanced state. The model cannot detect this mismatch.
- `config_kinematic_case` must have the same value in both namelists.
- The 2D cases use 5 levels with zero vertical flux, so each level is a copy of the 2D
  problem.
- `dcmip2012_1_2` extensions to DCMIP 2012: the q1 norms at T/2 against the exact
  solution and the mixing diagnostics. At T/2 the flow has no net translation in
  latitude–height, so the errors at T partly cancel; the T/2 norms show the error before
  the reversal. DCMIP defines both diagnostics only for test 1-1.
- The exact solution at T/2 is `dcmip2012_1_2_exact_T_half.npz`: q1 by backward RK4
  trajectories (1000 steps) on 9001 latitudes (poles included) for the layer midpoints of
  L30, L60 and L120 with uniform levels and ztop = 12 km, stored as float32. Its
  interpolation error is about 2.4e-5. Regenerate it with
  `python3 dcmip2012_1_2_exact.py --make-table` only if the script changes.
