"""
lg_ext_analyze.py - Step 13c: gate-length study with the sweep extended to +2 V [TCAD].
Run from the repo root:   python3 tcad_scripts/lg_study/lg_ext_analyze.py [folder]
  default folder: tcad_scripts/lg_study/plt   (holds IdVg_<name>_des.plt from Step 13b and IdVg_<name>x_des.plt
  from Step 13c, for the six extended runs g500 g700 a200 a350 a500 a700)
Outputs (next to this script): lg_study_ext.png, lg_study_ext.json.

Merge: Step 13b curve for V_GS < +0.9 V, Step 13c curve from +0.9 V (both runs have points at 0.8-1.0 V).
Check 1: over +0.80 .. +1.00 V, max |dI_D|/I_D between the two runs (log interpolation, as Step 13b); expect < 0.5 %.
Check 2: forward gate current. V_IG = first V_GS where |I_G| > 1 % of I_D [ASSUME criterion]; metrics are given
  on the full merged curve and, if V_IG lies inside the sweep, on the curve cut at V_IG ("valid").
Metrics: lg_analyze.metrics (Step 13b definitions, unchanged). Runs without an extension keep their Step 13b curve.
Main peak (24 Nov, user choice A): the abrupt runs a200-a700 turn on a second time at +1.3..+1.8 V (Step 13c), so the
  global g_m maximum is no longer the channel's peak. Main peak = first local maximum of the smoothed g_m that exceeds
  25 % of the global maximum; the curve is cut at the first local minimum after it that is followed by a rise of more
  than 2 % of that peak (the Step 13b GVS80b rule, applied from the first peak instead of the global one); metrics are
  then recomputed on the cut curve ("main"). Check: for g050 the main GVS80 must equal Step 13b GVS80b (1.44 V).
"""
import glob
import json
import os
import re
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lg_analyze import read_plt, metrics, gvs

HERE = os.path.dirname(os.path.abspath(__file__))
VMERGE, LIMIT = 0.9, 0.01


def curves(path):
    d = read_plt(path)
    v, i, g = d["gate OuterVoltage"], d["drain TotalCurrent"], d["gate TotalCurrent"]
    k = np.argsort(v); v, i, g = v[k], i[k], g[k]
    keep = np.r_[True, np.diff(v) > 1e-6]
    return v[keep], i[keep], g[keep]


def main_peak(cv, base):
    """Main-peak metrics on the smoothed g_m of the whole curve (as Step 13b GVS80b): first local maximum above
    25 % of the global maximum, cut at the first later local minimum followed by a rise > 2 % of that peak."""
    v, i, gms, vip3 = cv
    loc = [j for j in range(1, len(gms) - 1) if gms[j] >= gms[j - 1] and gms[j] >= gms[j + 1] and gms[j] > 0.25 * gms.max()]
    k = loc[0] if loc else int(np.argmax(gms))
    cut = len(gms) - 1
    for j in range(k + 1, len(gms) - 1):
        if gms[j] <= gms[j - 1] and gms[j] <= gms[j + 1] and gms[j:].max() - gms[j] > 0.02 * gms[k]:
            cut = j
            break
    w8, lo8, hi8, e8 = gvs(v[:cut + 1], gms[:cut + 1], 0.8)
    w9 = gvs(v[:cut + 1], gms[:cut + 1], 0.9)[0]
    return dict(base, gm_peak_mS_mm=float(gms[k] * 1e3), V_gm_peak_V=float(v[k]), Vth_V=float(v[k] - i[k] / gms[k]),
                GVS80_V=w8, range80_V=[lo8, hi8], GVS80_hits_sweep_end=bool(e8), GVS90_V=w9,
                VIP3_at_peak_V=float(vip3[k])), (float(v[cut]) if cut < len(gms) - 1 else None)


def overlap(va, ia, vb, ib, lo=0.80, hi=1.00):
    vv = np.arange(max(lo, va[0], vb[0]), min(hi, va[-1], vb[-1]) + 1e-9, 0.005)
    A = np.exp(np.interp(vv, va, np.log(np.abs(ia)))); B = np.exp(np.interp(vv, vb, np.log(np.abs(ib))))
    return float(np.max(np.abs(B / A - 1)) * 100), float(vv[0]), float(vv[-1])


if __name__ == "__main__":
    folder = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "plt")
    base = {re.search(r"IdVg_(.+)_des\.plt", f).group(1): f for f in glob.glob(os.path.join(folder, "IdVg_*_des.plt"))}
    ext = {n[:-1]: f for n, f in base.items() if n.endswith("x")}
    study = sorted([n for n in base if re.fullmatch(r"[ga]\d{3}", n)], key=lambda n: (n[0] != "g", n))
    print(f"{len(study)} Step 13b runs, {len(ext)} extensions ({', '.join(sorted(ext))}) in {folder}")
    out, plot = {}, {}
    for n in study:
        v, i, g = curves(base[n]); rec = {"extended": n in ext}
        if n in ext:
            vx, ix, gx = curves(ext[n])
            rec["check_overlap_maxrel_pct"], a, b = overlap(v, i, vx, ix)
            print(f"CHECK {n}: Step 13b vs 13c over {a:+.2f} .. {b:+.2f} V: max |dI|/I = "
                  f"{rec['check_overlap_maxrel_pct']:.3f} % (expect < 0.5 %)")
            m = v < VMERGE; mx = vx >= VMERGE
            v, i, g = np.r_[v[m], vx[mx]], np.r_[i[m], ix[mx]], np.r_[g[m], gx[mx]]
        on = np.abs(i) > 1e-3
        bad = on & (np.abs(g) > LIMIT * np.abs(i))
        vig = float(v[np.argmax(bad)]) if bad.any() else None
        rec["V_IG_1pct_V"] = vig
        rec["IG_over_ID_at_end"] = float(abs(g[-1] / i[-1]))
        rec["sweep_V"] = [float(v[0]), float(v[-1])]
        rec["all"], cv = metrics(v, i)
        if vig is not None and vig < v[-1] - 0.1:
            k = v < vig
            rec["valid"], cv = metrics(v[k], i[k])
            v, i = v[k], i[k]
        rec["main"], vc = main_peak(cv, rec.get("valid", rec["all"]))
        rec["V_main_cut_V"] = vc
        out[n], plot[n] = rec, (cv, vig, vc)

    hdr = (f"{'run':5s} {'L_G um':>6s} {'end V':>5s} {'V_IG':>6s} {'V_cut':>6s} {'g_m,pk':>7s} {'@V':>6s} {'V_th':>6s} "
           f"{'GVS80':>7s} {'GVS90':>6s} {'V_IP3':>6s} | {'global pk':>9s} {'@V':>6s}")
    print(hdr + "\n" + "-" * len(hdr))
    for n in study:
        r = out[n]; m = r["main"]; a = r.get("valid", r["all"]); f = lambda x: "-" if x is None else format(x, "+6.2f")
        print(f"{n:5s} {int(n[1:]) / 1000:6.3f} {r['sweep_V'][1]:5.1f} {f(r['V_IG_1pct_V']):>6s} {f(r['V_main_cut_V']):>6s} "
              f"{m['gm_peak_mS_mm']:7.1f} {m['V_gm_peak_V']:+6.2f} {m['Vth_V']:+6.2f} "
              f"{m['GVS80_V']:6.2f}{'*' if m['GVS80_hits_sweep_end'] else ' '}{m['GVS90_V']:6.2f} {m['VIP3_at_peak_V']:6.2f} | "
              f"{a['gm_peak_mS_mm']:9.1f} {a['V_gm_peak_V']:+6.2f}")
    print("main-peak metrics (g_m mS/mm, V); V_cut = dip after the main peak where the curve is cut ('-' = no second rise);\n"
          "V_IG = first V_GS with |I_G| > 1 % of I_D ('-' = never); * = GVS80 range reaches the end of the (cut) curve;\n"
          "global pk = maximum g_m of the whole curve, for reference")
    print("graded / abrupt at the same L_G, main peak (GVS80, GVS90 in V):")
    for L in sorted({int(n[1:]) for n in study}):
        g_, a_ = out.get(f"g{L:03d}"), out.get(f"a{L:03d}")
        if g_ and a_:
            G = g_["main"]; A = a_["main"]
            print(f"  {L:3d} nm: GVS80 {G['GVS80_V']:.2f} / {A['GVS80_V']:.2f}   GVS90 {G['GVS90_V']:.2f} / {A['GVS90_V']:.2f}"
                  f"   peak g_m {G['gm_peak_mS_mm']:.0f} / {A['gm_peak_mS_mm']:.0f} mS/mm")

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    cmap = plt.get_cmap("viridis"); lgs = sorted({int(n[1:]) for n in study})
    for n in study:
        (v, i, gms, _), vig, vc = plot[n]
        c = cmap(lgs.index(int(n[1:])) / max(1, len(lgs) - 1)); ls = "-" if n[0] == "g" else "--"
        ax[0].plot(v, gms * 1e3, ls, c=c, lw=1.6, label=f"{int(n[1:])} nm" if n[0] == "g" else None)
        if vig is not None:
            ax[0].axvline(vig, c=c, ls=":", lw=0.8)
        if vc is not None:
            ax[0].plot(vc, np.interp(vc, v, gms) * 1e3, "x", c=c, ms=7)
    ax[0].set(xlabel="V_GS (V)", ylabel="g_m (mS/mm)", xlim=(-3.5, 2.05),
              title="V_DS = 5 V; solid graded [P], dashed abrupt [ASSUME]; x: main-peak cut")
    ax[0].legend(frameon=False, fontsize=8, title="L_G")
    for kind, mk, lab in (("g", "o-", "graded"), ("a", "s--", "abrupt")):
        L = [x for x in lgs if f"{kind}{x:03d}" in out]
        f = lambda x, q: out[f"{kind}{x:03d}"]["main"][q]
        ax[1].plot(L, [f(x, "GVS80_V") for x in L], mk, c="#2a78d6", label=f"GVS80, {lab}")
        ax[1].plot(L, [f(x, "GVS90_V") for x in L], mk, c="#eb6834", label=f"GVS90, {lab}")
    ax[1].set(xscale="log", xlabel="L_G (nm)", ylabel="GVS (V)", title="Gate voltage swing (main peak), sweep to +2 V")
    for a in ax:
        a.grid(c="#e1e0d9"); a.spines[["top", "right"]].set_visible(False)
    ax[1].legend(frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(HERE, "lg_study_ext.png"), dpi=140)
    json.dump(out, open(os.path.join(HERE, "lg_study_ext.json"), "w"), indent=2)
    print("saved lg_study_ext.png, lg_study_ext.json")
