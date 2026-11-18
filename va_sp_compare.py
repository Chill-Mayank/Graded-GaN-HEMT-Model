"""
va_sp_compare.py - Step 10c: compare Spectre's S-parameters of graded_hemt.va (va/sp_device1.s2p) with rf_twin.py.
Run:  python3 va_sp_compare.py
This checks the Verilog-A charges and extrinsic elements in AC (Step 10a checked DC only).
Pass criterion [ASSUME, set before the run]: H21 and U within 0.05 dB, |Y| within 0.5 % at every frequency.
Touchstone 2-port data order is S11 S21 S12 S22 [Touchstone 1.1 specification, EIA/IBIS].
"""
import numpy as np
from rf_twin import RF
from va_netlists import SP1


def read_s2p(fn):
    unit, fmt, z0 = "GHZ", "MA", 50.0
    mult = {"HZ": 1, "KHZ": 1e3, "MHZ": 1e6, "GHZ": 1e9}
    rows = []
    for ln in open(fn):
        ln = ln.split("!")[0].strip()
        if not ln:
            continue
        if ln.startswith("#"):
            t = ln[1:].upper().split()
            for i, w in enumerate(t):
                if w in mult:
                    unit = w
                if w in ("MA", "DB", "RI"):
                    fmt = w
                if w == "R":
                    z0 = float(t[i + 1])
            continue
        rows += [float(x) for x in ln.split()]
    a = np.array(rows).reshape(-1, 9)
    f = a[:, 0] * mult[unit]
    p = a[:, 1:].reshape(-1, 4, 2)
    if fmt == "RI":
        v = p[..., 0] + 1j * p[..., 1]
    elif fmt == "MA":
        v = p[..., 0] * np.exp(1j * np.deg2rad(p[..., 1]))
    else:
        v = 10 ** (p[..., 0] / 20) * np.exp(1j * np.deg2rad(p[..., 1]))
    S = np.zeros((len(f), 2, 2), complex)
    S[:, 0, 0], S[:, 1, 0], S[:, 0, 1], S[:, 1, 1] = v[:, 0], v[:, 1], v[:, 2], v[:, 3]
    return f, S, z0


def s_to_y(S, z0):
    I = np.eye(2)
    return np.array([np.linalg.solve(I + s, I - s) / z0 for s in S])


if __name__ == "__main__":
    f, S, z0 = read_s2p("va/sp_device1.s2p")
    Ysp = s_to_y(S, z0)
    r = RF(1, rg=SP1["rg"], cgs=SP1["cgsfr"], cgd=SP1["cgdfr"], cds=SP1["cdsfr"])
    Ytw, I = r.y2(f, SP1["vgs"], SP1["vds"])
    hs, us = RF.gains(Ysp)
    ht, ut = RF.gains(Ytw)
    dh = np.max(np.abs(20 * np.log10(hs / ht)))
    du = np.max(np.abs(10 * np.log10(us / ut)))
    dy = np.max(np.abs(Ysp - Ytw) / np.abs(Ytw).max(axis=(1, 2))[:, None, None])
    print(f"Device 1, V_GS = {SP1['vgs']} V, V_DS = {SP1['vds']} V, {len(f)} frequencies {f[0] / 1e9:.1f}-{f[-1] / 1e9:.1f} GHz")
    print(f"  H21: worst difference {dh:.4f} dB   (criterion 0.05 dB)")
    print(f"  U  : worst difference {du:.4f} dB   (criterion 0.05 dB)")
    print(f"  Y  : worst difference {100 * dy:.3f} % of max |Y|   (criterion 0.5 %)")
    for k in (0, len(f) // 2, len(f) - 1):
        print(f"    {f[k] / 1e9:6.2f} GHz: H21 Spectre {20 * np.log10(hs[k]):6.2f} twin {20 * np.log10(ht[k]):6.2f} dB | "
              f"U Spectre {10 * np.log10(us[k]):6.2f} twin {10 * np.log10(ut[k]):6.2f} dB")
