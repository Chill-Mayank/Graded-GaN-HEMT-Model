"""lg_snap_lateral.py - Step 13d: where the a200 second turn-on starts, along the device [TCAD].
Run (needs h5py; Python <= 3.14 or the sandbox):  python3 lg_snap_lateral.py [folder with a200s_vg_NNNN_des.tdr]
Reads the a200 snapshots of lg_snap_des.cmd (V_GS = 1.70, 1.75, 1.775, 1.80 V: files 0002, 0003, 0004, 0006) and prints
  1. 2DEG sheet density n_s(y) (eDensity integrated over GaN + AlGaN depth) along the device;
  2. change of trapped charge per layer (NetTrappedCharge integrated over depth) between 1.75 and 1.775 V;
  3. surface-donor occupation (TrapOccupation_0(Do,Le), AlGaN/SiN interface) in both access regions.
Geometry [TCAD, sde_lg.cmd, L_G = 0.2 um]: source n+ y 0-0.5, source access 0.5-1.1, gate 1.1-1.3,
  drain access 1.3-1.9, drain n+ 1.9-2.4 um. Surface donors: 3e13 cm^-2, one level 0.36 eV from midgap [ASSUME, Synopsys example].
"""
import os
import sys
import numpy as np
import h5py

NV = {1: 2, 2: 3, 3: 4}
SNAP = {2: 1.70, 3: 1.75, 4: 1.775, 6: 1.80}


def load(path, want):
    f = h5py.File(path, "r"); g = f["collection/geometry_0"]; V = g["vertex"][:]; X, Y = V["x"], V["y"]
    names = {int(k.split("_")[1]): g[k].attrs["name"].decode() for k in g if k.startswith("region_")}

    def rverts(r):
        e = r["elements_0"][:]; idx, i = [], 0
        while i < len(e):
            n = NV.get(int(e[i]), 1); idx.extend(e[i + 1:i + 1 + n]); i += 1 + n
        return np.unique(np.array(idx))
    RV, out, st = {}, {}, g["state_0"]
    for k in st:
        ds = st[k]
        if not isinstance(ds, h5py.Group) or "quantity" not in ds.attrs:
            continue
        q = ds.attrs["quantity"].decode(); r = int(ds.attrs["region"])
        if q not in want or ds.attrs["location type"] != 0:
            continue
        if r not in RV:
            RV[r] = rverts(g[f"region_{r}"])
        out.setdefault(q, {})[names[r]] = (X[RV[r]], Y[RV[r]], ds["values"][:])
    return out


def column(dat, regions, ys):
    pts = {}
    for rn in regions:
        for x, y, v in zip(*dat.get(rn, ((), (), ()))):
            pts[(round(float(x), 7), round(float(y), 7))] = float(v)
    P = np.array([(x, y, v) for (x, y), v in pts.items()])
    res = []
    for y0 in ys:
        c = P[np.abs(P[:, 1] - y0) < 1e-6]
        c = c[np.argsort(c[:, 0])]
        res.append(np.trapezoid(c[:, 2], c[:, 0] * 1e-4) if len(c) > 1 else np.nan)
    return np.array(res)


if __name__ == "__main__":
    folder = sys.argv[1] if len(sys.argv) > 1 else "."
    want = {"eDensity", "NetTrappedCharge", "TrapOccupation_0(Do,Le)"}
    D = {k: load(os.path.join(folder, f"a200s_vg_{k:04d}_des.tdr"), want) for k in SNAP}
    ys_all = np.unique(np.round(D[2]["eDensity"]["channel"][1], 7))
    near = lambda ps: np.array([ys_all[np.argmin(abs(ys_all - p))] for p in ps])
    zone = lambda y: ("source access" if 0.5 < y < 1.1 else "gate" if 1.1 <= y <= 1.3 else "drain access" if 1.3 < y < 1.9 else "n+")
    ys = near([0.65, 0.8, 0.95, 1.05, 1.15, 1.2, 1.25, 1.35, 1.5, 1.7, 1.85])
    print("1. 2DEG n_s (1e12 cm^-2) along the device, V_GS = " + " / ".join(f"{v}" for v in SNAP.values()) + " V")
    NS = {k: column(D[k]["eDensity"], ["channel", "graded", "barrier"], ys) for k in SNAP}
    for j, y in enumerate(ys):
        print(f"   y {y:5.2f} um  {zone(y):13s} " + " ".join(f"{NS[k][j] / 1e12:6.3f}" for k in SNAP)
              + f"   step 1.75->1.775: {(NS[4][j] - NS[3][j]) / 1e10:+6.1f}e10")
    print("2. trapped charge change 1.75 -> 1.775 V, integrated over depth (q cm^-2), at y = 0.65 / 0.80 / 1.20 / 1.35 um")
    yt = near([0.65, 0.8, 1.2, 1.35])
    for r in ("barrier", "graded", "channel", "buffer"):
        a = column(D[3]["NetTrappedCharge"], [r], yt); b = column(D[4]["NetTrappedCharge"], [r], yt)
        print(f"   {r:8s} " + " ".join(f"{v:+.2e}" for v in b - a))
    print("   ('barrier' holds the AlGaN/SiN surface donors at its top vertices)")
    print("3. surface-donor occupation TrapOccupation_0(Do,Le) along the access regions")
    for side, ps in (("passivation_lump_1+barrier", (0.6, 0.7, 0.8, 0.9, 1.0)), ("barrier+passivation_lump_2", (1.4, 1.5, 1.6, 1.7, 1.8))):
        print(f"   {side} at y = " + " ".join(f"{p:.1f}" for p in ps) + " um")
        for k, v in SNAP.items():
            xx, yy, vv = D[k]["TrapOccupation_0(Do,Le)"][side]; o = np.argsort(yy)
            print(f"     V_GS {v:5.3f}: " + " ".join(f"{np.interp(p, yy[o], vv[o]):.4f}" for p in ps))
