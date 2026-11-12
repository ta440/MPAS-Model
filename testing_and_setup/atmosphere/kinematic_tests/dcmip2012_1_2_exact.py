"""Exact q1 for DCMIP 2012 test 1-2 (Hadley-like meridional circulation) by backward trajectories.

Usage: dcmip2012_1_2_exact.py [--levels 60] [--hours 12] [--nlat 360] [--steps 2000] [--save out.npz]
       dcmip2012_1_2_exact.py --make-table [dcmip2012_1_2_exact_T_half.npz]

--make-table writes the check file read by check_kinematic_test.py: q1 at T/2 = 12 h on
TABLE_NLAT latitudes for each level set in TABLE_LEVELS, with TABLE_STEPS RK4 steps, as float32
(about 3.5 min). Regenerate it only if this script changes.

The flow (DCMIP 2012 test case document v1.7, eqs. 36-38) has no longitude dependence apart from
the solid-body zonal wind, which only shifts q1 in longitude, and q1 (eq. 43) is zonally uniform.
So q1(lat, z, t) = q1_initial(z at the departure point) of the backward trajectory in (lat, z)
from t to 0, integrated with RK4. Points are the layer midpoints of an equidistant grid with
ztop = 12 km (L30/L60/L120 as in DCMIP Table VI) and nlat equally spaced latitudes.

Prints, per level: the fraction of latitudes with q1 > 1e-3, max q1, and max |dq1/dz|.
--save writes lat, z and q1 for plotting exact latitude-height sections.
"""
import argparse
import os
import sys

import numpy as np

A = 6371229.0          # MPAS sphere radius (m); q1(lat, z) does not depend on it
G = 9.80616            # MPAS gravity (m s-2)
RD = 287.0             # MPAS dry-air gas constant (J kg-1 K-1)
T0 = 300.0
H = RD * T0 / G        # scale height (m)
ZTOP = 12000.0
W0 = 0.15
K = 5
TAU = 86400.0
Z1, Z2 = 2000.0, 5000.0
Z0 = 0.5 * (Z1 + Z2)


def velocity(phi, z, t):
    """(dlat/dt, dz/dt) from eqs. (37)-(38); rho0/rho = exp(z/H)."""
    r = np.exp(z / H)
    c = np.cos(np.pi * t / TAU)
    v = -A * W0 * np.pi / (K * ZTOP) * r * np.cos(phi) * np.sin(K * phi) * np.cos(np.pi * z / ZTOP) * c
    w = W0 / K * r * (-2.0 * np.sin(K * phi) * np.sin(phi) + K * np.cos(phi) * np.cos(K * phi)) \
        * np.sin(np.pi * z / ZTOP) * c
    return v / A, w


def q1_initial(z):
    """Eq. (43)."""
    return np.where((z > Z1) & (z < Z2), 0.5 * (1.0 + np.cos(2.0 * np.pi * (z - Z0) / (Z2 - Z1))), 0.0)


def exact_q1(lat, z, t, steps):
    """q1 at (lat, z) and time t, from RK4 trajectories integrated back to t = 0."""
    p, zz = lat.copy(), z.copy()
    h = -t / steps
    for _ in range(steps):
        k1 = velocity(p, zz, t)
        k2 = velocity(p + 0.5 * h * k1[0], zz + 0.5 * h * k1[1], t + 0.5 * h)
        k3 = velocity(p + 0.5 * h * k2[0], zz + 0.5 * h * k2[1], t + 0.5 * h)
        k4 = velocity(p + h * k3[0], zz + h * k3[1], t + h)
        p += h / 6.0 * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0])
        zz += h / 6.0 * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1])
        t += h
    return q1_initial(zz), zz


TABLE_LEVELS = (30, 60, 120)     # DCMIP Table VI, with the 240 / 120 / 60 km meshes
TABLE_NLAT = 9001                # 0.02 deg; max interpolation error 2.4e-5 (step-3 draft, 2026-11-10)
TABLE_STEPS = 1000               # agrees with 2000 steps to 4e-12
TABLE_HOURS = 12.0


def make_table(path):
    """q1 at T/2 on a fine latitude grid (poles included) for each level set in TABLE_LEVELS."""
    lat = np.radians(np.linspace(-90.0, 90.0, TABLE_NLAT))
    out = {"lat_deg": np.degrees(lat), "hours": TABLE_HOURS, "steps": TABLE_STEPS, "ztop": ZTOP,
           "levels": np.array(TABLE_LEVELS)}
    for nlev in TABLE_LEVELS:
        zc = (np.arange(nlev) + 0.5) * ZTOP / nlev
        P, Z = np.meshgrid(lat, zc, indexing="ij")
        q, _ = exact_q1(P, Z, TABLE_HOURS * 3600.0, TABLE_STEPS)
        out[f"z_L{nlev}"] = zc
        out[f"q1_L{nlev}"] = q.astype(np.float32)
        print(f"L{nlev}: done")
    np.savez_compressed(path, **out)
    print(f"saved {path}")


def main():
    if len(sys.argv) >= 2 and sys.argv[1] == "--make-table":
        default = os.path.join(os.path.dirname(os.path.abspath(__file__)), "dcmip2012_1_2_exact_T_half.npz")
        make_table(sys.argv[2] if len(sys.argv) > 2 else default)
        return
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--levels", type=int, default=60)
    ap.add_argument("--hours", type=float, default=12.0)
    ap.add_argument("--nlat", type=int, default=360)
    ap.add_argument("--steps", type=int, default=2000)
    ap.add_argument("--save", default=None)
    args = ap.parse_args()

    dz = ZTOP / args.levels
    zc = (np.arange(args.levels) + 0.5) * dz
    dlat = 180.0 / args.nlat
    lat = np.radians(np.linspace(-90.0 + 0.5 * dlat, 90.0 - 0.5 * dlat, args.nlat))
    P, Z = np.meshgrid(lat, zc, indexing="ij")

    q, zdep = exact_q1(P, Z, args.hours * 3600.0, args.steps)

    frac = (q > 1e-3).sum(axis=0) / q.shape[0]
    grad = np.abs(np.gradient(q, dz, axis=1)).max(axis=0)
    print(f"Exact q1 at t = {args.hours:g} h, L{args.levels}, {args.nlat} latitudes, {args.steps} RK4 steps")
    print("level  z (km)  frac(q>0)   max q   max|dq/dz| (1/km)")
    for k in range(args.levels):
        print(f"{k+1:5d} {zc[k]/1e3:7.2f} {frac[k]:10.3f} {q[:, k].max():7.3f} {grad[k]*1e3:12.2f}")
    occ = np.where(frac > 0)[0]
    if occ.size:
        print(f"tracer on levels {occ[0]+1}-{occ[-1]+1} ({zc[occ[0]]/1e3:.2f}-{zc[occ[-1]]/1e3:.2f} km)")
    init = np.where(q1_initial(zc) > 0)[0]
    print(f"initial layer on levels {init[0]+1}-{init[-1]+1}")
    print(f"departure heights span {zdep.min():.1f}-{zdep.max():.1f} m")

    if args.save:
        np.savez(args.save, lat=np.degrees(lat), z=zc, q1=q, hours=args.hours)
        print(f"saved {args.save}")


if __name__ == "__main__":
    main()
