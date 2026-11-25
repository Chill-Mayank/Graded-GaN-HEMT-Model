"""lg_snap.py - Step 13d: where the extra electrons go at the second turn-on of the abrupt control [TCAD].
Run where h5py is installed, in the folder holding the snapshots (sh run_lg.sh snap a200 on the server):
    python3 lg_snap.py a200
  h5py: not on the server's Python, and no build for Python 3.15 on the Mac (25 Nov); Python <= 3.14 or the sandbox.
Reads <name>s_vg_NNNN_des.tdr (HDF5, h5py; written by lg_snap_des.cmd at the V_GS values in VGS below, in order);
prints one line per snapshot and writes <name>s_snap.csv (per-layer sheet densities) and <name>s_cuts.csv
(vertical profiles at the gate centre: depth, electron density, E_C, E_C - E_Fn), small files to send back.
Per layer (GaN "channel", AlGaN "graded" and "barrier"; in the abrupt control both AlGaN layers are x = 0.27):
  n_s = integral of eDensity over the layer depth, at the gate centre and averaged over the gate footprint (cm^-2);
  min(E_C - E_Fn) in the layer at the gate centre (eV), E_Fn = -q eQuasiFermiPotential: below ~0 the layer is
  degenerate, i.e. holds a channel of its own.
Geometry [TCAD, sde_lg.cmd]: gate from y = 1.10 um to 1.10 + L_G; layers from the GaN top (x = -0.94 um) up:
  graded 6 nm, barrier 9 nm, gate on the barrier at x = -0.955 um.
Hypothesis tested [ASSUME]: the second turn-on of a200-a700 is electrons filling the x = 0.27 AlGaN under the gate.
  Result (25 Nov): rejected; AlGaN n_s < 2e8 cm^-2 at all 8 snapshots. See lg_snap_lateral.py for the real cause.
"""
import glob
import sys
import numpy as np
import h5py

# must match lg_snap_des.cmd: V_GS = 0.8 + 1.2 t at Time t
TIMES = [0.16667, 0.66667, 0.75, 0.79167, 0.8125, 0.82083, 0.83333, 0.91667]
VGS = [round(0.8 + 1.2 * t, 3) for t in TIMES]          # 1.00 1.60 1.70 1.75 1.775 1.785 1.80 1.90 V
NV = {1: 2, 2: 3, 3: 4}
LAYERS = ("channel", "graded", "barrier")
YG0 = 1.10                                              # um, gate left edge [TCAD, sde_lg.cmd]


def load(path):
    f = h5py.File(path, "r")
    g = f["collection/geometry_0"]
    V = g["vertex"][:]
    X, Y = V["x"], V["y"]

    def rverts(r):
        e = r["elements_0"][:]
        idx, i = [], 0
        while i < len(e):
            n = NV[int(e[i])]
            idx.extend(e[i + 1:i + 1 + n])
            i += 1 + n
        return np.unique(np.array(idx))

    names = {int(k.split("_")[1]): g[k].attrs["name"].decode()
             for k in g if k.startswith("region_") and "material" in g[k].attrs}
    RV = {r: rverts(g[f"region_{r}"]) for r in names}
    st = g["state_0"]
    fields = {}
    for k in st:
        ds = st[k]
        if not isinstance(ds, h5py.Group) or "quantity" not in ds.attrs:
            continue
        r = int(ds.attrs["region"])
        if r in names and ds.attrs["location type"] == 0:
            fields.setdefault(ds.attrs["quantity"].decode(), {})[r] = ds["values"][:]
    return X, Y, names, RV, fields


def cut(D, q, yc, layer):
    X, Y, names, RV, fields = D
    rows = {}
    for r, vals in fields[q].items():
        if names[r] != layer:
            continue
        s = abs(Y[RV[r]] - yc) < 1e-6
        for x, v in zip(X[RV[r]][s], vals[s]):
            rows.setdefault(round(float(x), 7), float(v))
    xs = np.array(sorted(rows))
    return xs, np.array([rows[x] for x in xs])


def sheet(D, yc, layer):
    xs, n = cut(D, "eDensity", yc, layer)
    return float(np.trapezoid(n, xs * 1e-4)) if len(xs) > 1 else 0.0


if __name__ == "__main__":
    name = sys.argv[1] if len(sys.argv) > 1 else "a200"
    prefix = sys.argv[2] if len(sys.argv) > 2 else f"{name}s_vg"
    lg = int(name[1:]) / 1000 if name[1:].isdigit() else 0.05
    files = sorted(glob.glob(f"{prefix}_[0-9]*_des.tdr"))
    if not files:
        sys.exit(f"no {prefix}_NNNN_des.tdr files in this folder")
    if len(files) != len(VGS):
        print(f"note: {len(files)} snapshots, {len(VGS)} V_GS values in the deck; V_GS labels may not match")
    snap = open(f"{name}s_snap.csv", "w")
    snap.write("file,VGS_V," + ",".join(f"ns_{l}_centre_cm-2,ns_{l}_footprint_cm-2,min_EcEf_{l}_eV" for l in LAYERS) + "\n")
    cuts = open(f"{name}s_cuts.csv", "w")
    cuts.write("file,VGS_V,layer,x_um,eDensity_cm-3,Ec_eV,Ec_minus_Efn_eV\n")
    print(f"{name}: L_G {lg} um, gate centre y = {YG0 + lg / 2:.3f} um")
    print(f"{'V_GS':>6s}  " + "  ".join(f"{l + ' n_s':>16s} {'E_C-E_F':>7s}" for l in LAYERS) + "   (n_s cm^-2 centre, E_C-E_Fn eV)")
    for k, fn in enumerate(files):
        D = load(fn)
        X, Y = D[0], D[1]
        Yall = np.unique(np.round(Y, 7))
        yc = Yall[np.argmin(abs(Yall - (YG0 + lg / 2)))]
        foot = [y for y in Yall if YG0 - 1e-9 <= y <= YG0 + lg + 1e-9]
        vg = VGS[k] if k < len(VGS) else float("nan")
        row, txt = [fn, f"{vg}"], []
        for l in LAYERS:
            nc = sheet(D, yc, l)
            nf = float(np.mean([sheet(D, y, l) for y in foot]))
            xe, ec = cut(D, "ConductionBandEnergy", yc, l)
            xq, phi = cut(D, "eQuasiFermiPotential", yc, l)
            d = ec + np.interp(xe, xq, phi)                   # E_C - E_Fn, E_Fn = -q phi_n
            xs, n = cut(D, "eDensity", yc, l)
            for x, v in zip(xs, n):
                cuts.write(f"{fn},{vg},{l},{x:.7f},{v:.6e},{np.interp(x, xe, ec):.6f},{np.interp(x, xe, d):.6f}\n")
            row += [f"{nc:.6e}", f"{nf:.6e}", f"{d.min():.4f}"]
            txt.append(f"{nc:9.3e} ({nf:9.3e}) {d.min():+7.3f}")
        snap.write(",".join(row) + "\n")
        print(f"{vg:+6.3f}  " + "  ".join(txt))
    snap.close(); cuts.close()
    print(f"wrote {name}s_snap.csv and {name}s_cuts.csv")
    print("columns per layer: n_s at the gate centre (footprint mean), min E_C - E_Fn at the centre")
