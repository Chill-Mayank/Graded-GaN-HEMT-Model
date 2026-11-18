"""
fom.py - Step 13a: linearity figures of merit for both devices (beyond the paper).
Run:  python3 fom.py      -> fom_linearity.png, fom_results.json
Definitions (the literature uses both GVS thresholds, so both are given):
  GVS80 = V_GS range where g_m >= 0.8 g_m,max  [review, J. Semicond. 44(12), 121801, 2023]
  GVS90 = V_GS range where g_m >= 0.9 g_m,max  ["10 % drop": Hsu et al., Electronics 12(20), 4336, 2023]
  V_IP3 = sqrt(4 g_m / (3 |g_m3|)), g_m3 = (1/6) d^3 I_D / dV_GS^3: amplitude where the extrapolated IM3 equals the
          fundamental for a weakly nonlinear transfer y = a1 x + a3 x^3 (IM3/fund = (3/4)(a3/a1) A^2) [standard,
          e.g. S. Arar, All About Circuits, "A guide to calculating IM3 and IP3", 2025]
  f_T(V_GS) from the small-signal two-port (rf_twin.py), read at |H21| = 1; GVS_fT80 = range with f_T >= 0.8 f_T,max.
Sources compared: this work (Verilog-A twin), the paper's model lines and the measured points [P, Figs. 6(a), 6(b), 6(c)].
Biases: Device 1 V_DS = 10 V (Fig. 6(a)); Device 2 V_DS = 5 V (Fig. 6(b)); f_T at the Figs. 9 / 11 bias (3 V / 10 V).
"""
import json
import os
os.environ.setdefault("GAN_DEVICE", "2")            # fit2_dibl.py (Device 2 RF card) requires this
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from va_twin import Twin, card
from rf_twin import RF

L = lambda f: np.loadtxt(f, delimiter=",", skiprows=3, ndmin=2)
DEV = {1: dict(vds=10.0, v=np.arange(-2.0, 3.0 + 1e-9, 0.01), gm="paper_data/fig6a_{}_gm.csv"),
       2: dict(vds=5.0, v=np.arange(-2.0, 1.0 + 1e-9, 0.005), gm="paper_data/fig6b_{}_gm.csv")}


def gvs(v, g, frac):
    """Width (V) of the contiguous range around the peak where g >= frac * max, by linear interpolation."""
    k = int(np.argmax(g)); thr = frac * g[k]
    i = k
    while i > 0 and g[i - 1] >= thr:
        i -= 1
    j = k
    while j < len(g) - 1 and g[j + 1] >= thr:
        j += 1
    lo = v[i] if i == 0 else v[i - 1] + (thr - g[i - 1]) * (v[i] - v[i - 1]) / (g[i] - g[i - 1])
    hi = v[j] if j == len(g) - 1 else v[j] + (thr - g[j]) * (v[j + 1] - v[j]) / (g[j + 1] - g[j])
    return hi - lo, lo, hi


def model_curves(dev):
    d = DEV[dev]
    tw = Twin({**card(dev), "SIGMA": 0.0})          # the paper's model: no DIBL
    I = np.array([tw.id_ext(x, d["vds"]) for x in d["v"]])
    if dev == 2:
        I = I / (card(2)["W"] * 1e3)                   # A -> A/mm, as Fig. 6(b) is per mm [P]
    g1 = np.gradient(I, d["v"])
    g3 = np.gradient(np.gradient(g1, d["v"]), d["v"]) / 6.0
    return I, g1, g3


if __name__ == "__main__":
    out, fig = {}, plt.figure(figsize=(13, 8))
    for col, dev in enumerate((1, 2)):
        d = DEV[dev]
        I, g1, g3 = model_curves(dev)
        vip3 = np.sqrt(4 * np.abs(g1) / (3 * np.maximum(np.abs(g3), 1e-12)))
        res = {}
        print(f"===== Device {dev}, V_DS = {d['vds']} V")
        for name, (v, g) in {"this work": (d["v"], g1 * 1e3),
                             "paper's model": tuple(L(d["gm"].format("paper_model")).T),
                             "measured": tuple(L(d["gm"].format("measured")).T)}.items():
            w8, lo8, hi8 = gvs(v, g, 0.8); w9, lo9, hi9 = gvs(v, g, 0.9)
            res[name] = {"GVS80_V": w8, "GVS90_V": w9, "range80_V": [lo8, hi8], "gm_peak_V": float(v[np.argmax(g)])}
            print(f"  {name:14s} GVS80 {w8:5.2f} V ({lo8:+.2f} .. {hi8:+.2f})   GVS90 {w9:5.2f} V   peak g_m at {v[np.argmax(g)]:+.2f} V")
        # V_IP3: sweet spot (g_m3 = 0 above threshold) and values at the g_m peak
        on = g1 > 0.2 * g1.max()
        on[:4] = on[-4:] = False                         # drop finite-difference end points
        zc = [d["v"][i] for i in range(len(g3) - 1) if on[i] and np.sign(g3[i]) != np.sign(g3[i + 1])]
        kp = int(np.argmax(g1))
        res["this work"].update({"VIP3_at_gm_peak_V": float(vip3[kp]), "gm3_zero_crossings_V": [float(z) for z in zc]})
        print(f"  V_IP3 at the g_m peak ({d['v'][kp]:+.2f} V): {vip3[kp]:.2f} V; g_m3 = 0 at {', '.join(f'{z:+.2f}' for z in zc)} V")
        if dev == 1:
            vq = -0.33532
            kq = int(np.argmin(abs(d["v"] - vq)))
            res["this work"]["VIP3_at_fig10_bias_V"] = float(vip3[kq])
            print(f"  V_IP3 at the Fig. 10 bias (V_GS = {vq} V, here at V_DS = 10 V): {vip3[kq]:.2f} V")
        if dev == 2:
            # measured and paper V_IP3 from the printed Fig. 6(c) G_M and G_M3 (G_M3 ~ d2G_M/dV2 smoothed over +-0.23 V,
            # printed item 14, so these V_IP3 are high by roughly sqrt(2) near the G_M3 peaks)
            for name, tag in (("measured", "measured"), ("paper's model", "paper_model")):
                a, b = L(f"paper_data/fig6c_{tag}_GM.csv"), L(f"paper_data/fig6c_{tag}_GM3.csv")
                vv = np.arange(-1.3, 0.8, 0.01)
                G = np.interp(vv, *a.T); G3 = np.interp(vv, *b.T) / 6.0
                vi = np.sqrt(4 * np.abs(G) / (3 * np.maximum(np.abs(G3), 1e-9)))
                kk = int(np.argmax(G)); res[name]["VIP3_at_gm_peak_V_from_fig6c"] = float(vi[kk])
                print(f"  {name}: V_IP3 at its g_m peak from printed Fig. 6(c): {vi[kk]:.2f} V (smoothed G_M3)")
        # f_T against V_GS at the RF bias of each device
        if dev == 1:
            r = RF(1, rg=17.5, cgd=16.7e-15); vds_rf, vgrid = 10.0, np.arange(-1.0, 3.01, 0.1)
        else:
            from fit2_dibl import BASE, KEY
            P = json.load(open("params_device2_dibl.json"))["fitted [CALIBRATE]"]; R = json.load(open("rf_device2.json"))
            ov = {**{k: BASE[k] for k in ("VOFF", "SIGMA", "VDSREF")}, **{KEY[n]: x for n, x in P.items()},
                  "DVOFF": R["DVOFF_V"], "U0": R["U0_m2Vs"]}
            r = RF(2, rg=R["RG_ohm"], cgs=R["CGSFR_F"], cgd=R["CGDFR_F"], cds=R["CDSFR_F"], over=ov)
            vds_rf, vgrid = 3.0, np.arange(-2.4, 1.01, 0.05)
        ft = np.array([r.ft_fmax(x, vds_rf)[0] / 1e9 for x in vgrid])
        ok = np.isfinite(ft)
        w, lo, hi = gvs(vgrid[ok], ft[ok], 0.8)
        res["this work"].update({"fT_max_GHz": float(np.nanmax(ft)), "GVS_fT80_V": w})
        print(f"  f_T (|H21| = 1) at V_DS = {vds_rf} V: max {np.nanmax(ft):.1f} GHz; f_T >= 80 % over {w:.2f} V ({lo:+.2f} .. {hi:+.2f})")
        out[f"device{dev}"] = res

        a1 = fig.add_subplot(2, 2, 1 + col); a2 = fig.add_subplot(2, 2, 3 + col)
        m, p = L(d["gm"].format("measured")), L(d["gm"].format("paper_model"))
        sc = 1.0
        a1.plot(m[:, 0], m[:, 1], "o", mfc="none", c="#52514e", ms=5, label="measured [P]")
        a1.plot(*p.T, "-", c="#2a78d6", lw=1.5, label="paper's model [P]")
        a1.plot(d["v"], g1 * 1e3 * sc, "--", c="#eb6834", lw=2, label="this work")
        for name, c in (("measured", "#52514e"), ("this work", "#eb6834")):
            lo, hi = res[name]["range80_V"]
            a1.axvspan(lo, hi, color=c, alpha=0.08)
        a1.set(xlabel="V_GS (V)", ylabel=f"g_m ({'mS' if dev == 1 else 'mS/mm'})",
               title=f"Device {dev}, V_DS = {d['vds']:g} V: g_m and GVS80 (shaded)")
        a1.legend(frameon=False, fontsize=8)
        a2.semilogy(d["v"][on], vip3[on], "-", c="#eb6834", lw=2, label="this work")
        a2.axvline(d["v"][kp], c="#52514e", lw=1, ls=":", label="g_m peak")
        if dev == 1:
            a2.axvline(-0.33532, c="#2a78d6", lw=1, ls="--", label="Fig. 10 bias")
        a2.legend(frameon=False, fontsize=8)
        a2.set(xlabel="V_GS (V)", ylabel="V_IP3 (V)", ylim=(0.3, 100),
               title=f"Device {dev}: V_IP3 = sqrt(4 g_m / 3|g_m3|), this work")
        for a in (a1, a2):
            a.grid(c="#e1e0d9"); a.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    plt.savefig("fom_linearity.png", dpi=140)
    json.dump(out, open("fom_results.json", "w"), indent=2)
    print("saved fom_linearity.png, fom_results.json")
