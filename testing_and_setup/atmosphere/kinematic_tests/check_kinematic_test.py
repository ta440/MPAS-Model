#!/usr/bin/env python3
"""Check an MPAS-A kinematic transport test run.

Usage: check_kinematic_test.py <run_dir>

Reads namelist.atmosphere, the history files and log.atmosphere.0000.out in <run_dir>.
Prints a table of error norms and checks, writes <run_dir>/results.json, and exits with
status 0 if every check passes and 1 otherwise. Needs numpy and netCDF4 only.

Tracer n (1-based) has shape mod(n-1, 4): 0 uniform, 1 Gaussian hills, 2 cosine bells,
3 slotted cylinders (Lauritzen et al. 2012). Error norms are eqs. (26)-(28) of that paper,
area-weighted, on level 1.
"""
import glob
import json
import os
import re
import sys

import numpy as np
from netCDF4 import Dataset

HERE = os.path.dirname(os.path.abspath(__file__))

SHAPES = ["uniform", "gaussian_hills", "cosine_bells", "slotted_cylinders"]

# Per case: period T in days (the exact solution at T is the initial state); whether all
# levels must stay bitwise identical (2D cases run as copies on every level); and whether
# the model holds rho_zz = 1. In that case the t = 0 history file still has the density of
# the init file, because kinematic mode sets rho_zz = 1 just before the first transport
# step, so the mass check uses rho_zz = 1 and separately checks that it is 1 at t > 0.
CASES = {
    "lauritzen_nondiv": {"period_days": 12.0, "levels_identical": True, "rho_zz_one": True},
}

TOL_UNIFORM = 1.0e-13      # max |q - 1| for the uniform tracer
TOL_MASS = 1.0e-12         # relative change in tracer mass (roundoff grows with the step count: 1.7e-13 at 60 km, 12 days)
TOL_BOUNDS = 1.0e-13       # overshoot allowed beyond the initial range, relative to the range
TOL_REGRESSION = 0.01      # relative change in l2 allowed against references.json


def read_namelist(path):
    """Return the options of a Fortran namelist file as a dict of strings (lower-case keys)."""
    opts = {}
    for line in open(path):
        m = re.match(r"\s*(config_\w+)\s*=\s*(.+?)\s*$", line)
        if m:
            opts[m.group(1).lower()] = m.group(2).strip("'\"")
    return opts


def duration_days(s):
    """'DD_hh:mm:ss' or 'hh:mm:ss' -> days."""
    days, _, hms = s.rpartition("_")
    h, m, sec = (float(x) for x in hms.split(":"))
    return (float(days) if days else 0.0) + (h + m / 60.0 + sec / 3600.0) / 24.0


def elapsed_days(file_a, file_b):
    """Days between the xtime of two history files (both within the same year)."""
    def day_of_year(path):
        with Dataset(path) as d:
            s = b"".join(d["xtime"][0]).decode().strip()     # 'YYYY-MM-DD_hh:mm:ss'
        y, mo, rest = int(s[0:4]), int(s[5:7]), s[8:]
        dd, hms = rest.split("_")
        h, m, sec = (float(x) for x in hms.split(":"))
        leap = (y % 4 == 0 and y % 100 != 0) or y % 400 == 0
        months = [31, 29 if leap else 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        return sum(months[:mo - 1]) + int(dd) - 1 + (h + m / 60.0 + sec / 3600.0) / 24.0
    return day_of_year(file_b) - day_of_year(file_a)


def norms(q, qe, w):
    d = q - qe
    l1 = np.sum(w * np.abs(d)) / np.sum(w * np.abs(qe))
    l2 = np.sqrt(np.sum(w * d**2) / np.sum(w * qe**2))
    linf = np.max(np.abs(d)) / np.max(np.abs(qe))
    return float(l1), float(l2), float(linf)


def read_timers(log_path):
    """Total times (s) of the timers of interest from the MPAS timer table of task 0."""
    want = {"total time", "time integration", "passive_tracer_transport"}
    out = {}
    if not os.path.exists(log_path):
        return out
    pat = re.compile(r"^\s*\d+\s+(.+?)\s{2,}([\d.]+)\s+(\d+)")
    for line in open(log_path):
        m = pat.match(line)
        if m and m.group(1).strip() in want and m.group(1).strip() not in out:
            out[m.group(1).strip()] = float(m.group(2))
    return out


def main(run_dir):
    nml = read_namelist(os.path.join(run_dir, "namelist.atmosphere"))
    case = nml.get("config_kinematic_case", "off")
    if case not in CASES:
        sys.exit(f"Unknown or missing config_kinematic_case '{case}' in namelist.atmosphere")
    spec = CASES[case]
    monotonic = nml.get("config_monotonic", "true").lower() in ("true", ".true.", "t")
    dt = float(nml["config_dt"])
    rk_order = int(nml.get("config_time_integration_order", "2"))   # MPAS default is 2
    ntracers = int(nml["config_num_passive_tracers"])
    run_days = duration_days(nml["config_run_duration"])

    files = sorted(glob.glob(os.path.join(run_dir, "history.*.nc")))
    if len(files) < 2:
        sys.exit(f"Need at least two history files in {run_dir}, found {len(files)}")

    with Dataset(files[0]) as g:
        area = g["areaCell"][:]
        ncells = area.size
        nlev = g.dimensions["nVertLevels"].size
        # The model conserves sum(area * rho_zz * q * dzeta), with dzeta = 1/rdzw the
        # thickness in the vertical coordinate (the same in every column). Without rdzw in
        # the output, equal level weights are used, which is exact when no vertical flux.
        dz_level = 1.0 / g["rdzw"][:] if "rdzw" in g.variables else np.ones(nlev)
    dzeta = np.broadcast_to(dz_level, (ncells, nlev))
    mesh = os.path.basename(glob.glob(os.path.join(run_dir, "*.grid.nc"))[0]).replace(".grid.nc", "")

    def load(fn):
        with Dataset(fn) as d:
            return d["passive_tracers"][0], d["rho_zz"][0]     # (nCells, nLev, N), (nCells, nLev)

    q0, r0 = load(files[0])
    if spec["rho_zz_one"]:
        r0 = np.ones_like(r0)
    mass0 = np.einsum("c,ck,ckn->n", area, r0 * dzeta, q0)
    lo0, hi0 = q0.min(axis=(0, 1)), q0.max(axis=(0, 1))

    checks = {}

    def record(name, ok, value):
        prev = checks.get(name)
        if prev is None or (prev["pass"] and not ok) or (ok == prev["pass"] and value > prev["value"]):
            checks[name] = {"pass": bool(ok), "value": float(value)}

    for i, fn in enumerate(files):
        q, r = load(fn)
        if spec["rho_zz_one"]:
            if i > 0:
                e = float(np.max(np.abs(r - 1.0)))
                record("rho_zz_equals_1", e == 0.0, e)
            r = np.ones_like(r)
        uni = [n for n in range(ntracers) if n % 4 == 0]
        record("uniform_tracer", True, 0.0)
        if uni:
            e = float(np.max(np.abs(q[:, :, uni] - 1.0)))
            record("uniform_tracer", e < TOL_UNIFORM, e)
        mass = np.einsum("c,ck,ckn->n", area, r * dzeta, q)
        rel = float(np.max(np.abs(mass - mass0) / np.abs(mass0)))
        record("mass_conservation", rel < TOL_MASS, rel)
        if monotonic:
            rng = np.maximum(hi0 - lo0, 1.0)            # the uniform tracer has zero range
            under = np.maximum(lo0 - q.min(axis=(0, 1)), 0.0) / rng
            above = np.maximum(q.max(axis=(0, 1)) - hi0, 0.0) / rng
            over = float(np.max(np.maximum(under, above)))
            record("shape_preservation", over < TOL_BOUNDS, over)
        if ntracers > 4:
            same = all(np.array_equal(q[:, :, n], q[:, :, n - 4]) for n in range(4, ntracers))
            record("identical_copies", same, 0.0 if same else 1.0)
        if spec["levels_identical"] and nlev > 1:
            same = all(np.array_equal(q[:, k, :], q[:, 0, :]) for k in range(1, nlev))
            record("identical_levels", same, 0.0 if same else 1.0)

    # Error norms at T, if the run reaches T and the last history file is at T
    errors = {}
    at_T = abs(elapsed_days(files[0], files[-1]) - spec["period_days"]) < 1e-9
    if at_T:
        qT, _ = load(files[-1])
        for s in range(1, min(4, ntracers)):
            errors[SHAPES[s]] = dict(zip(("l1", "l2", "linf"), norms(qT[:, 0, s], q0[:, 0, s], area)))

    # Regression against references.json
    regression = None
    ref_path = os.path.join(HERE, "references.json")
    if at_T and os.path.exists(ref_path):
        key = f"{case}/{mesh}/dt{dt:g}/rk{rk_order}/{'mono' if monotonic else 'nomono'}"
        ref = json.load(open(ref_path)).get("runs", {}).get(key)
        if ref:
            worst = max(abs(errors[s]["l2"] - ref[s]["l2"]) / ref[s]["l2"] for s in errors)
            record("regression_l2", worst < TOL_REGRESSION, worst)
            regression = key

    # Cost from the MPAS timers (task 0)
    timers = read_timers(os.path.join(run_dir, "log.atmosphere.0000.out"))
    nranks = 1
    log0 = os.path.join(run_dir, "log.atmosphere.0000.out")
    if os.path.exists(log0):
        m = re.search(r"MPI task count:\s*(\d+)", open(log0).read())
        nranks = int(m.group(1)) if m else 1
    nsteps = int(round(run_days * 86400.0 / dt))
    cost = None
    if "passive_tracer_transport" in timers and nsteps > 0:
        cost = timers["passive_tracer_transport"] * nranks / (ncells * nlev * ntracers * nsteps) * 1e6

    # Report
    print(f"\n{case} on {mesh}: dt = {dt:g} s, RK order {rk_order}, {nlev} levels, {ntracers} tracers, "
          f"{run_days:g} days, limiter {'on' if monotonic else 'off'}")
    if errors:
        print(f"\n  {'shape':<20}{'l1':>11}{'l2':>11}{'linf':>11}")
        for s, e in errors.items():
            print(f"  {s:<20}{e['l1']:11.3e}{e['l2']:11.3e}{e['linf']:11.3e}")
    else:
        print("\n  (run does not end at T: no error norms or regression check)")
    print(f"\n  {'check':<22}{'value':>11}  result")
    for name, c in checks.items():
        print(f"  {name:<22}{c['value']:11.2e}  {'pass' if c['pass'] else 'FAIL'}")
    if regression:
        print(f"  (regression reference: {regression})")
    if cost is not None:
        print(f"\n  passive-tracer transport: {timers['passive_tracer_transport']:.2f} s on task 0, "
              f"{cost:.4f} core-us per cell x level x tracer x step ({nranks} ranks)")

    ok = all(c["pass"] for c in checks.values())
    print(f"\n  {'ALL CHECKS PASSED' if ok else 'SOME CHECKS FAILED'}\n")

    results = {"case": case, "mesh": mesh, "dt": dt, "rk_order": rk_order, "nlev": nlev, "ntracers": ntracers,
               "run_days": run_days, "monotonic": monotonic, "nranks": nranks,
               "errors_at_T": errors, "checks": checks, "regression_key": regression,
               "timers_s": timers, "cost_core_us_per_cell_level_tracer_step": cost, "pass": ok}
    with open(os.path.join(run_dir, "results.json"), "w") as f:
        json.dump(results, f, indent=2)
    return 0 if ok else 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1]))
