"""n8_ns.py - Step 6b-5: electron sheet density under the gate of Device 2 from the n8 TDR snapshots.
Run on the server, in GaN_Dev2:   python3 n8_ns.py
Reads n8_vg*.tdr (HDF5) with h5py; writes two small text files to send back:
  n8_ns.csv    one row per snapshot: file, n_s at the gate centre (cm^-2), n_s averaged over the gate footprint (cm^-2)
  n8_cuts.csv  vertical profiles at the gate centre: file, depth x (um), electron density (cm^-3), E_C (eV)
Geometry [TCAD, n5 deck]: gate contact y = 1.100-1.150 um, on the barrier top at x = -0.955 um; the semiconductor
runs to x = 0 (buffer from x = -0.9 um). n_s integrates the electron density from the barrier top to x = -0.85 um,
as for the n3-n5 equilibrium numbers (tdr_extract.py)."""
import glob, sys
import numpy as np
import h5py

NV = {1: 2, 2: 3, 3: 4}
SEMI = {"channel", "graded", "barrier", "buffer"}
Y_GATE = (1.100, 1.150)     # um, gate footprint [TCAD]
X_STOP = -0.85              # um, integrate electrons above this depth


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
        q = ds.attrs["quantity"].decode()
        r = int(ds.attrs["region"])
        if r in names and ds.attrs["location type"] == 0:
            fields.setdefault(q, {})[r] = ds["values"][:]
    return X, Y, names, RV, fields


def cut(X, Y, names, RV, fields, q, yc):
    rows = {}
    for r, vals in fields[q].items():
        if names[r] not in SEMI:
            continue
        s = abs(Y[RV[r]] - yc) < 1e-6
        for x, v in zip(X[RV[r]][s], vals[s]):
            rows.setdefault(round(float(x), 7), float(v))
    xs = np.array(sorted(rows))
    return xs, np.array([rows[x] for x in xs])


def ns_of(xs, n):
    m = xs < X_STOP
    return np.trapezoid(n[m], xs[m] * 1e-4) if hasattr(np, "trapezoid") else np.trapz(n[m], xs[m] * 1e-4)


files = sorted(glob.glob("n8_vg*.tdr"))
if not files:
    sys.exit("no n8_vg*.tdr files in this folder")
out = open("n8_ns.csv", "w")
out.write("file,ns_centre_cm-2,ns_footprint_mean_cm-2,n_mesh_columns\n")
cuts = open("n8_cuts.csv", "w")
cuts.write("file,x_um,eDensity_cm-3,Ec_eV\n")
for fn in files:
    X, Y, names, RV, fields = load(fn)
    Yall = np.unique(np.round(Y, 7))
    yc = Yall[np.argmin(abs(Yall - 0.5 * sum(Y_GATE)))]
    xs, n = cut(X, Y, names, RV, fields, "eDensity", yc)
    xe, ec = cut(X, Y, names, RV, fields, "ConductionBandEnergy", yc)
    cols = [y for y in Yall if Y_GATE[0] - 1e-9 <= y <= Y_GATE[1] + 1e-9]
    ns_cols = [ns_of(*cut(X, Y, names, RV, fields, "eDensity", y)) for y in cols]
    out.write(f"{fn},{ns_of(xs, n):.6e},{np.mean(ns_cols):.6e},{len(cols)}\n")
    for x, v in zip(xs, n):
        if x < X_STOP:
            cuts.write(f"{fn},{x:.6f},{v:.6e},{np.interp(x, xe, ec):.6f}\n")
    print(f"{fn}: n_s centre {ns_of(xs, n):.3e} cm^-2, footprint mean {np.mean(ns_cols):.3e} ({len(cols)} columns)")
out.close()
cuts.close()
print("wrote n8_ns.csv and n8_cuts.csv")
