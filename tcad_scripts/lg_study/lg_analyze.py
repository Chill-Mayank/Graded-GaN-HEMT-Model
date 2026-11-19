"""
lg_analyze.py - Step 13b: gate-length study on Device 2, TCAD transfer curves at V_DS = 5 V [TCAD].
Run from the repo root:   python3 tcad_scripts/lg_study/lg_analyze.py [folder with IdVg_*_des.plt]
  default folder: tcad_scripts/lg_study/plt
Reads every IdVg_<name>_des.plt (name g050..g700 graded [P], a050..a700 abrupt control [ASSUME]);
the n7 run (IdVg_n7_des.plt, same folder, optional) is the reference for the g050 check.
Outputs (next to this script): lg_study.png, lg_study.json.

Definitions (as in fom.py, Step 13a):
  I_D in A/mm (AreaFactor 1e3 um in the decks); g_m = dI_D/dV_GS on a uniform 5 mV grid (cubic interpolation
  of the solver points); g_m3 = (1/6) d3I/dV3, from g_m smoothed with a Savitzky-Golay filter (0.1 V window,
  cubic) [NUM: TCAD points are up to 21 mV apart, so the third derivative needs smoothing].
  GVS80 / GVS90: V_GS range with g_m >= 80 % / 90 % of the peak; V_IP3 = sqrt(4 g_m / 3|g_m3|).
  V_th: linear extrapolation of I_D at the g_m peak; V_OFF: V_GS at I_D = 1 mA/mm (same criterion as n7);
  SS: minimum dV_GS/dlog10(I_D) for 1e-7 < I_D < 1e-3 A/mm.
"""
import glob
import json
import os
import re
import sys
import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))


def read_plt(path):
    """DF-ISE xyplot -> dict of dataset name -> array."""
    t = open(path).read()
    names = re.findall(r'"([^"]+)"', t[t.index("datasets"):t.index("functions")])
    vals = np.array(t[t.index("Data {") + 6:t.rindex("}")].split(), float)
    vals = vals[: len(vals) // len(names) * len(names)].reshape(-1, len(names))
    return {n: vals[:, i] for i, n in enumerate(names)}


def transfer(path):
    d = read_plt(path)
    v, i = d["gate OuterVoltage"], d["drain TotalCurrent"]
    k = np.argsort(v); v, i = v[k], i[k]
    keep = np.r_[True, np.diff(v) > 1e-6]
    return v[keep], i[keep]


def gvs(v, g, frac):
    k = int(np.argmax(g)); thr = frac * g[k]; i = j = k
    while i > 0 and g[i - 1] >= thr:
        i -= 1
    while j < len(g) - 1 and g[j + 1] >= thr:
        j += 1
    lo = v[i] if i == 0 else v[i - 1] + (thr - g[i - 1]) * (v[i] - v[i - 1]) / (g[i] - g[i - 1])
    hi = v[j] if j == len(g) - 1 else v[j] + (thr - g[j]) * (v[j + 1] - v[j]) / (g[j + 1] - g[j])
    return hi - lo, lo, hi, (j == len(g) - 1)


def metrics(v0, i0):
    h = 0.005
    v = np.arange(v0[0], v0[-1] + 1e-9, h)
    i = interp1d(v0, i0, kind="cubic")(v)
    gm = np.gradient(i, h)
    gms = savgol_filter(gm, 21, 3)                         # 0.1 V window
    gm3 = savgol_filter(gm, 21, 3, deriv=2, delta=h) / 6.0
    k = int(np.argmax(gms))
    vth = v[k] - i[k] / gms[k]
    on = i > 1e-3
    voff = float(np.interp(1e-3, i[~on][-50:] if (~on).sum() > 50 else i[~on], v[~on][-50:] if (~on).sum() > 50 else v[~on])) if on.any() and (~on).any() else float("nan")
    sub = (i > 1e-7) & (i < 1e-3)
    ss = float(np.min(np.diff(v[sub]) / np.diff(np.log10(i[sub]))) * 1e3) if sub.sum() > 5 else float("nan")
    w8, lo8, hi8, e8 = gvs(v, gms, 0.8); w9, lo9, hi9, e9 = gvs(v, gms, 0.9)
    vip3 = np.sqrt(4 * np.abs(gms) / (3 * np.maximum(np.abs(gm3), 1e-12)))
    return dict(gm_peak_mS_mm=float(gms[k] * 1e3), V_gm_peak_V=float(v[k]), Vth_V=float(vth), VOFF_1mA_V=voff,
                SS_mV_dec=ss, ID_at_VGS0_A_mm=float(np.interp(0.0, v, i)),
                GVS80_V=w8, range80_V=[lo8, hi8], GVS80_hits_sweep_end=bool(e8), GVS90_V=w9,
                VIP3_at_peak_V=float(vip3[k])), (v, i, gms, vip3)


if __name__ == "__main__":
    folder = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "plt")
    files = sorted(glob.glob(os.path.join(folder, "IdVg_*_des.plt")))
    runs = {re.search(r"IdVg_(.+)_des\.plt", f).group(1): f for f in files}
    print(f"{len(runs)} runs in {folder}: {', '.join(runs)}")
    out, curves = {}, {}
    for name, f in runs.items():
        v0, i0 = transfer(f)
        out[name], curves[name] = metrics(v0, i0)
        out[name]["sweep_V"] = [float(v0[0]), float(v0[-1])]
    if "n7" in runs and "g050" in runs:                        # the check: g050 must reproduce n7
        va, ia = transfer(runs["n7"]); vb, ib = transfer(runs["g050"])
        vv = np.arange(max(va[0], vb[0]) + 0.05, min(va[-1], vb[-1]) - 0.05, 0.01)
        A, B = np.interp(vv, va, ia), np.interp(vv, vb, ib)
        on = A > 1e-3
        rel = float(np.max(np.abs(B[on] / A[on] - 1)) * 100) if on.any() else float("nan")
        print(f"CHECK g050 vs n7 over {vv[0]:+.2f} .. {vv[-1]:+.2f} V: max |dI|/I = {rel:.3f} % where I_D > 1 mA/mm"
              f" (expect < 0.5 %)")
        out["check_g050_vs_n7_maxrel_pct"] = rel
    hdr = f"{'run':5s} {'L_G um':>6s} {'g_m,pk':>7s} {'@V':>6s} {'V_th':>6s} {'V_OFF':>6s} {'SS':>5s} {'GVS80':>6s} {'GVS90':>6s} {'V_IP3':>6s}"
    print(hdr + "\n" + "-" * len(hdr))
    for name in sorted(runs, key=lambda n: (n[0] != "g", n)):
        m = out[name]; lg = int(name[1:]) / 1000 if name[0] in "ga" and name[1:].isdigit() else float("nan")
        print(f"{name:5s} {lg:6.3f} {m['gm_peak_mS_mm']:7.1f} {m['V_gm_peak_V']:+6.2f} {m['Vth_V']:+6.2f} {m['VOFF_1mA_V']:+6.2f}"
              f" {m['SS_mV_dec']:5.0f} {m['GVS80_V']:6.2f}{'*' if m['GVS80_hits_sweep_end'] else ' '}{m['GVS90_V']:5.2f} {m['VIP3_at_peak_V']:6.2f}")
    print("units: g_m mS/mm, voltages V, SS mV/dec; * = GVS80 range reaches the end of the sweep (+1 V), a lower bound")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    cmap = plt.get_cmap("viridis")
    study = [n for n in runs if n[0] in "ga" and n[1:].isdigit()]
    lgs = sorted({int(n[1:]) for n in study})
    for n in study:
        v, i, g, _ = curves[n]
        c = cmap(lgs.index(int(n[1:])) / max(1, len(lgs) - 1))
        ax[0].plot(v - out[n]["Vth_V"], g * 1e3, "-" if n[0] == "g" else "--", c=c, lw=1.6,
                   label=f"{int(n[1:])} nm" if n[0] == "g" else None)
    ax[0].set(xlabel="V_GS - V_th (V)", ylabel="g_m (mS/mm)", title="V_DS = 5 V; solid graded [P], dashed abrupt [ASSUME]")
    ax[0].legend(frameon=False, fontsize=8, title="L_G")
    for kind, mk, lab in (("g", "o-", "graded"), ("a", "s--", "abrupt")):
        L = [x for x in lgs if f"{kind}{x:03d}" in out]
        if not L:
            continue
        ax[1].plot(L, [out[f"{kind}{x:03d}"]["GVS80_V"] for x in L], mk, c="#2a78d6", label=f"GVS80, {lab}")
        ax[1].plot(L, [out[f"{kind}{x:03d}"]["GVS90_V"] for x in L], mk, c="#eb6834", label=f"GVS90, {lab}")
        ax[2].plot(L, [out[f"{kind}{x:03d}"]["gm_peak_mS_mm"] for x in L], mk, c="#1baf7a", label=f"peak g_m, {lab}")
    ax[1].set(xscale="log", xlabel="L_G (nm)", ylabel="GVS (V)", title="Gate voltage swing")
    ax[2].set(xscale="log", xlabel="L_G (nm)", ylabel="peak g_m (mS/mm)", title="Peak transconductance")
    for a in ax:
        a.grid(c="#e1e0d9"); a.spines[["top", "right"]].set_visible(False)
    for a in ax[1:]:
        a.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(HERE, "lg_study.png"), dpi=140)
    json.dump(out, open(os.path.join(HERE, "lg_study.json"), "w"), indent=2)
    print("saved lg_study.png, lg_study.json")
