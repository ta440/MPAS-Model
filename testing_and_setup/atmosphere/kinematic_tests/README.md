# Kinematic transport tests for MPAS-Atmosphere

These tests run the production MPAS-A scalar transport with the dynamics switched off
and the mass fluxes prescribed from an analytic flow, so that results can be compared
with known exact solutions. Everything else (mesh, decomposition, halos, I/O, timers)
is the normal model.

Location in the fork: `testing_and_setup/atmosphere/kinematic_tests/`.

## Cases

| `config_kinematic_case` | Flow | Exact solution | Status |
|---|---|---|---|
| `lauritzen_nondiv` | Lauritzen et al. (2012, GMD) Case 4: non-divergent deformation with a background solid-body rotation, T = 12 days | the initial state at t = T | available |
| `lauritzen_div` | Lauritzen et al. (2012) Case 3: divergent deformation | initial state at T | planned |
| `dcmip_1_2` | DCMIP-2012 1-2: Hadley-like overturning | initial state at T | planned |
| `dcmip_1_1` | DCMIP-2012 1-1: 3D deformation | initial state at T | planned |
| `terminator` | Lauritzen et al. (2015) terminator toy chemistry | conservation of Cl + 2 Cl2 | planned |

Each case has a directory with namelist and streams templates. The tracers are the
passive-tracer array (`config_num_passive_tracers`); tracer `n` has shape `mod(n-1,4)`:
0 uniform (q = 1), 1 Gaussian hills, 2 cosine bells, 3 slotted cylinders. Tracers `n`
and `n+4` are identical copies, which must stay bitwise identical.

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

  | Mesh | Spacing | Default `dt` (s) |
  |---|---|---|
  | `x1.10242` | 240 km | 900 |
  | `x1.40962` | 120 km | 450 |
  | `x1.163842` | 60 km | 225 |

  Unpack each into its own directory, for example `meshes/x1.10242/`. If a partition
  file for your rank count is missing, the run script makes one with `gpmetis` when it
  is on the path.
- Python 3 with `numpy` and `netCDF4` for the checker.

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
| `-t` | time step (s) | by mesh, table above |
| `-k` | Runge–Kutta order, `config_time_integration_order` (2 or 3) | 2 |
| `-L` | FCT limiter, `config_monotonic` (`on` or `off`) | `on` |
| `-l` | vertical levels | 5 |
| `-N` | passive tracers | 4 |
| `-r` | run duration | `12_00:00:00` |
| `-i` | history output interval | `3_00:00:00` |
| `-o` | run directory | `./runs/<case>_<mesh>_N<ntracers>_rk<order>_lim<on/off>` |

The script writes the namelists and streams from the case templates, runs
`init_atmosphere_model` (init case 14) and then `atmosphere_model`, and runs the
checker. The exit code is the checker's: 0 if all checks pass.

## Suites

```
./run_suite.sh quick        -d /path/to/meshes      # 240 km, 4 tracers, 12 days; about 30 s on 4 cores
./run_suite.sh convergence  -d /path/to/meshes      # 240, 120, 60 km
./run_suite.sh scaling      -d /path/to/meshes      # 120 km, N = 4, 11, 50, 100, 200, 1 day
```

`-d` here is the parent directory holding one subdirectory per mesh; `-b`, `-n`, `-k`
and `-L` are passed on to every run. The scaling suite
covers the range that matters operationally: about 11 tracers in a default MPAS run, and
up to 200 for chemistry.

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
| Uniform tracer | max \|q − 1\| < 1e-13 at every output time |
| Mass conservation | relative change in Σ A ρ_zz q Δζ < 1e-12 for every tracer |
| Shape preservation (if `config_monotonic`) | no value outside the initial global range by more than 1e-13 |
| Identical copies | tracer n bitwise equal to tracer n−4 at every output time |
| Identical levels (2D cases) | every level bitwise equal to level 1 |
| Density (cases with ρ_zz = 1) | ρ_zz exactly 1 after t = 0 (the t = 0 file holds the init density) |
| Regression (if a reference exists, and the run reaches T) | l2 error of each shape within 1% of `references.json` |

It also prints the normalised l1, l2 and l∞ errors at T (Lauritzen et al. 2012,
eqs. 26–28) and the transport cost from the MPAS timers, and writes everything to
`results.json` in the run directory.

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
