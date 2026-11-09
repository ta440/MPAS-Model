#!/usr/bin/env python3
"""Check an MPAS-A kinematic transport test run.

Usage: check_kinematic_test.py <run_dir>
       check_kinematic_test.py --compare <run_dir_a> <run_dir_b>

The first form reads namelist.atmosphere, the history files and log.atmosphere.0000.out
in <run_dir>, prints error norms, mixing diagnostics and checks, writes
<run_dir>/results.json, and exits with status 0 if every check passes and 1 otherwise.
The second form checks that two runs (for example on 1 and 4 MPI ranks) have bitwise
identical history files in every field. Needs numpy and netCDF4 only.

Tracer n (1-based) has shape mod(n-1, 5): 0 uniform, 1 Gaussian hills, 2 cosine bells,
3 slotted cylinders, 4 correlated cosine bells (Lauritzen et al. 2012, GMD 5, 887-901).
Error norms are their eqs. (26)-(28), area-weighted, on level 1. The mixing diagnostics
(their Sect. 3.5 and Appendix C) use the pair (cosine bells, correlated cosine bells).
"""
import glob
import json
import os
import re
import sys

import numpy as np
from netCDF4 import Dataset

HERE = os.path.dirname(os.path.abspath(__file__))

SHAPES = ["uniform", "gaussian_hills", "cosine_bells", "slotted_cylinders", "correlated_cosine_bells"]
NSHAPES = len(SHAPES)
CHI, XI = 2, 4                 # shapes of the correlated pair

# Per case:
#   period_days       T; the exact solution at T is the initial state
#   levels_identical  2D cases run as copies on every level, which must stay bitwise identical
#   rho0_one          the model sets rho_zz = 1 before the first transport step; the t = 0
#                     history file still holds the density of the init file, so the mass
#                     baseline uses rho_zz = 1
#   rho_const         rho_zz must stay exactly 1 (non-divergent flow)
#   rho_exact_at_T    exact rho_zz at T, for the density error (divergent flow)
#   mixing            compute the mixing diagnostics at T/2
CASES = {
    "lauritzen_nondiv": {"period_days": 12.0, "levels_identical": True, "rho0_one": True,
                         "rho_const": True, "rho_exact_at_T": None, "mixing": True},
    "lauritzen_div":    {"period_days": 12.0, "levels_identical": True, "rho0_one": True,
                         "rho_const": False, "rho_exact_at_T": 1.0, "mixing": False},
}

TOL_UNIFORM = 1.0e-12      # max |q - 1| for the uniform tracer (roundoff grows with the step count: 1.3e-13 for lauritzen_div at 60 km)
TOL_MASS = 1.0e-12         # relative change in tracer mass (roundoff grows with the step count: 1.7e-13 at 60 km, 12 days)
TOL_BOUNDS = 1.0e-12       # overshoot allowed beyond the initial range, relative to the range
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


def file_day(path):
    """Day of the year (fractional, from 0) of the xtime of a history file."""
    with Dataset(path) as d:
        s = b"".join(d["xtime"][0]).decode().strip()     # 'YYYY-MM-DD_hh:mm:ss'
    y, mo, rest = int(s[0:4]), int(s[5:7]), s[8:]
    dd, hms = rest.split("_")
    h, m, sec = (float(x) for x in hms.split(":"))
    leap = (y % 4 == 0 and y % 100 != 0) or y % 400 == 0
    months = [31, 29 if leap else 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    return sum(months[:mo - 1]) + int(dd) - 1 + (h + m / 60.0 + sec / 3600.0) / 24.0


def norms(q, qe, w):
    d = q - qe
    l1 = np.sum(w * np.abs(d)) / np.sum(w * np.abs(qe))
    l2 = np.sqrt(np.sum(w * d**2) / np.sum(w * qe**2))
    linf = np.max(np.abs(d)) / np.max(np.abs(qe))
    return float(l1), float(l2), float(linf)


def mixing_diagnostics(chi, xi, w, a=-0.8, b=0.9, chi_min=0.1, chi_max=1.0):
    """Mixing diagnostics l_r, l_u, l_o of Lauritzen et al. (2012), Appendix C.

    chi, xi: cosine bells and correlated cosine bells per cell; w: cell areas.
    The pre-existing relation is xi = psi(chi) = a chi^2 + b on [chi_min, chi_max].
    d_k is the shortest normalised distance from (chi_k, xi_k) to that curve (their eqs.
    B1-B5), found from the real roots of d/dchi of the squared distance (a cubic) and the
    two end points.
    """
    psi = lambda c: a * c**2 + b
    xi_min, xi_max = psi(chi_max), psi(chi_min)
    r_chi, r_xi = chi_max - chi_min, xi_max - xi_min

    # d/dc [((c - chi)/r_chi)^2 + ((a c^2 + b - xi)/r_xi)^2] = 0
    #   <=>  c^3 + p c + q = 0,  p = (2/r_chi^2 + 4a(b - xi)/r_xi^2) / (4a^2/r_xi^2),  q = (-2 chi/r_chi^2) / (4a^2/r_xi^2)
    c3 = 4.0 * a**2 / r_xi**2
    p = (2.0 / r_chi**2 + 4.0 * a * (b - xi) / r_xi**2) / c3
    q = (-2.0 * chi / r_chi**2) / c3
    comp = np.zeros((chi.size, 3, 3))
    comp[:, 1, 0] = 1.0
    comp[:, 2, 1] = 1.0
    comp[:, 0, 2] = -q
    comp[:, 1, 2] = -p
    roots = np.linalg.eigvals(comp)
    real = np.where(np.abs(roots.imag) < 1e-9, roots.real, np.nan)
    real = np.where((real >= chi_min) & (real <= chi_max), real, np.nan)
    cand = np.concatenate([real, np.full((chi.size, 1), chi_min), np.full((chi.size, 1), chi_max)], axis=1)
    dist2 = ((cand - chi[:, None]) / r_chi)**2 + ((psi(cand) - xi[:, None]) / r_xi)**2
    d = np.sqrt(np.nanmin(dist2, axis=1))

    in_range = (chi >= chi_min) & (chi <= chi_max)
    chord = xi_max + (xi_min - xi_max) * (chi - chi_min) / (chi_max - chi_min)
    in_A = in_range & (xi >= chord) & (xi <= psi(chi))                       # "real" mixing
    in_box = in_range & (xi >= xi_min) & (xi <= xi_max)
    in_B = in_box & ~in_A                                                    # range-preserving unmixing
    in_O = ~in_box                                                           # overshooting
    total = np.sum(w)
    return {"l_r": float(np.sum(w * d * in_A) / total),
            "l_u": float(np.sum(w * d * in_B) / total),
            "l_o": float(np.sum(w * d * in_O) / total)}


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


def compare(run_a, run_b):
    """Bitwise comparison of all floating-point fields of the history files of two runs."""
    files_a = sorted(glob.glob(os.path.join(run_a, "history.*.nc")))
    if not files_a:
        sys.exit(f"No history files in {run_a}")
    differing = []
    for fa in files_a:
        fb = os.path.join(run_b, os.path.basename(fa))
        if not os.path.exists(fb):
            differing.append(f"{os.path.basename(fa)}: missing in {run_b}")
            continue
        with Dataset(fa) as a, Dataset(fb) as b:
            a.set_auto_mask(False)
            b.set_auto_mask(False)
            for v in a.variables:
                if v in b.variables and a[v].dtype.kind == "f" and not np.array_equal(a[v][:], b[v][:]):
                    differing.append(f"{os.path.basename(fa)}: {v}")
    print(f"\n  compared {len(files_a)} history files of\n    {run_a}\n    {run_b}")
    if differing:
        print("  NOT IDENTICAL:\n    " + "\n    ".join(differing) + "\n")
        return 1
    print("  all fields bitwise identical\n")
    return 0


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
    days = [file_day(f) - file_day(files[0]) for f in files]

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
    if spec["rho0_one"]:
        r0 = np.ones_like(r0)
    mass0 = np.einsum("c,ck,ckn->n", area, r0 * dzeta, q0)
    lo0, hi0 = q0.min(axis=(0, 1)), q0.max(axis=(0, 1))

    checks = {}

    def record(name, ok, value):
        prev = checks.get(name)
        if prev is None or (prev["pass"] and not ok) or (ok == prev["pass"] and value > prev["value"]):
            checks[name] = {"pass": bool(ok), "value": float(value)}

    uniform = [n for n in range(ntracers) if n % NSHAPES == 0]
    for i, fn in enumerate(files):
        q, r = load(fn)
        if i == 0 and spec["rho0_one"]:
            r = r0
        if spec["rho_const"] and i > 0:
            e = float(np.max(np.abs(r - 1.0)))
            record("rho_zz_equals_1", e == 0.0, e)
        if uniform:
            e = float(np.max(np.abs(q[:, :, uniform] - 1.0)))
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
        if ntracers > NSHAPES:
            same = all(np.array_equal(q[:, :, n], q[:, :, n - NSHAPES]) for n in range(NSHAPES, ntracers))
            record("identical_copies", same, 0.0 if same else 1.0)
        if spec["levels_identical"] and nlev > 1:
            same = all(np.array_equal(q[:, k, :], q[:, 0, :]) for k in range(1, nlev))
            record("identical_levels", same, 0.0 if same else 1.0)

    period = spec["period_days"]

    # Error norms at T (and the density error, where the exact density at T is known)
    errors, rho_error = {}, None
    at_T = abs(days[-1] - period) < 1e-9
    if at_T:
        qT, rT = load(files[-1])
        for s in range(1, min(NSHAPES, ntracers)):
            errors[SHAPES[s]] = dict(zip(("l1", "l2", "linf"), norms(qT[:, 0, s], q0[:, 0, s], area)))
        if spec["rho_exact_at_T"] is not None:
            d = rT[:, 0] - spec["rho_exact_at_T"]
            rho_error = {"max_abs": float(np.max(np.abs(d))),
                         "rms": float(np.sqrt(np.sum(area * d**2) / np.sum(area)))}

    # Mixing diagnostics at T/2 for the correlated pair (first copy)
    mixing = None
    half = [f for f, t in zip(files, days) if abs(t - 0.5 * period) < 1e-9]
    if spec["mixing"] and half and ntracers > XI:
        qh, _ = load(half[0])
        mixing = mixing_diagnostics(qh[:, 0, CHI], qh[:, 0, XI], area)

    # Regression against references.json
    regression = None
    ref_path = os.path.join(HERE, "references.json")
    if at_T and os.path.exists(ref_path):
        key = f"{case}/{mesh}/dt{dt:g}/rk{rk_order}/{'mono' if monotonic else 'nomono'}"
        ref = json.load(open(ref_path)).get("runs", {}).get(key)
        if ref:
            common = [s for s in errors if s in ref]
            worst = max(abs(errors[s]["l2"] - ref[s]["l2"]) / ref[s]["l2"] for s in common)
            record("regression_l2", worst < TOL_REGRESSION, worst)
            regression = key

    # Cost from the MPAS timers (task 0)
    log0 = os.path.join(run_dir, "log.atmosphere.0000.out")
    timers = read_timers(log0)
    nranks = 1
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
        print(f"\n  {'shape at T':<26}{'l1':>11}{'l2':>11}{'linf':>11}")
        for s, e in errors.items():
            print(f"  {s:<26}{e['l1']:11.3e}{e['l2']:11.3e}{e['linf']:11.3e}")
    else:
        print("\n  (run does not end at T: no error norms or regression check)")
    if rho_error:
        print(f"\n  density at T (exact {spec['rho_exact_at_T']:g}): max error {rho_error['max_abs']:.3e}, "
              f"rms error {rho_error['rms']:.3e}")
    if mixing:
        print(f"\n  mixing at T/2, pair (cosine bells, correlated cosine bells): "
              f"l_r {mixing['l_r']:.3e}  l_u {mixing['l_u']:.3e}  l_o {mixing['l_o']:.3e}")
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
               "errors_at_T": errors, "density_error_at_T": rho_error, "mixing_at_T_half": mixing,
               "checks": checks, "regression_key": regression,
               "timers_s": timers, "cost_core_us_per_cell_level_tracer_step": cost, "pass": ok}
    with open(os.path.join(run_dir, "results.json"), "w") as f:
        json.dump(results, f, indent=2)
    return 0 if ok else 1


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "--compare":
        sys.exit(compare(sys.argv[2], sys.argv[3]))
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1]))
