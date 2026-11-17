"""
fig6c.py - Step 9: Device 2 G_M, G_M2 and G_M3 at V_DS = 5 V against the paper's Fig. 6(c).
Run:  GAN_DEVICE=2 python3 fig6c.py   (needs params_device2.json from fit2.py, Step 8c)
      -> fig6c_device2.png, fig6c_gm3_definition.png

Nothing is fitted here: the model is the Step 8c Device 2 model, fitted to Fig. 6(b) only.

Definitions. The paper plots G_M, G_M2 and G_M3 on one axis and does not define G_M2 or G_M3 [P, Fig. 6(c)].
Its own curves fix them:
  G_M2 = dG_M/dV_GS           the paper's G_M2 line is the derivative of its G_M line (checked below)
  G_M3 = central difference of G_M2 over +-H, with H inferred from the paper's own G_M2 and G_M3 lines [CALC]
         (the printed G_M3 is not the exact d^2 G_M/dV_GS^2: its peaks are about 0.6x the exact ones; printed item 14)
This work is shown both ways: processed like the paper (compared with the data) and exact (for reference).
Units: per mm of gate width (printed axis "S"; printed item 11), G_M2 in S/(mm V), G_M3 in S/(mm V^2).
"""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

assert os.environ.get("GAN_DEVICE") == "2", "run as: GAN_DEVICE=2 python3 fig6c.py"
import fit2 as f2

VDS = 5.0                     # V [P, Fig. 6(c)]
W_MM = f2.W_MM                # mm, 0.075 [P]
STEP = 0.005                  # V, V_GS grid for the model derivatives [ASSUME: numerical step]
V_LO = -1.8                   # V, comparison window starts where the paper's G_M and G_M2 lines start [ASSUME]


def load(name):
    return f2.load(name)


def line(name):
    """A paper model line, sorted and without repeated V_GS, as (V, y)."""
    a = load(name)
    o = np.argsort(a[:, 0])
    v, i = np.unique(a[o, 0], return_index=True)
    return v, a[o, 1][i]


# ---- 1. Fix the paper's definitions from its own model lines (no model of ours involved)
pg1, pg2, pg3 = line("fig6c_paper_model_GM"), line("fig6c_paper_model_GM2"), line("fig6c_paper_model_GM3")
vchk = np.arange(-1.75, 0.70 + 1e-9, 0.005)
P1 = lambda v: np.interp(v, *pg1, left=0.0)
P2 = lambda v: np.interp(v, *pg2, left=0.0)
vi = np.arange(-1.835, 0.70 + 1e-9, 0.005)    # integrate G_M2 from where the G_M line starts (extraction noise
cum = P1(vi[0]) + np.r_[0.0, np.cumsum(0.5 * (P2(vi[1:]) + P2(vi[:-1])) * np.diff(vi))]   # is smoothed, not amplified)
rms_gm2_def = np.sqrt(np.mean((cum - P1(vi)) ** 2))
scan = np.round(np.arange(0.10, 0.35 + 1e-9, 0.005), 3)
err = [np.sqrt(np.mean(((P2(vchk + h) - P2(vchk - h)) / (2 * h) - np.interp(vchk, *pg3)) ** 2)) for h in scan]
H = float(scan[int(np.argmin(err))])          # V [CALC from the paper's own lines]
exact3 = (P2(vchk + 0.02) - P2(vchk - 0.02)) / 0.04   # near-exact derivative of the paper's G_M2 line (+-20 mV)

# ---- 2. This work: Step 8c model, I_D per mm on a fine grid, then derivatives
p = list(json.load(open("params_device2.json"))["fitted [CALIBRATE]"].values())
V = np.arange(-3.0 - H - 0.02, 1.0 + H + 0.02 + 1e-9, STEP)
I = np.array([f2.model_id(p, v, VDS) for v in V]) * 1e-3 / W_MM     # A/mm
GM = np.gradient(I, V)                          # S/mm
GM2 = np.gradient(GM, V)                        # S/(mm V)
GM3x = np.gradient(GM2, V)                      # S/(mm V^2), exact
k = int(round(H / STEP))
GM3p = np.full_like(GM2, np.nan)
GM3p[k:-k] = (GM2[2 * k:] - GM2[:-2 * k]) / (V[2 * k:] - V[:-2 * k])   # processed like the paper
keep = (V >= -3.0 - 1e-9) & (V <= 1.0 + 1e-9)
V, GM, GM2, GM3x, GM3p = V[keep], GM[keep], GM2[keep], GM3x[keep], GM3p[keep]

# ---- 3. Compare at the measured points
meas = {q: load(f"fig6c_measured_{q}") for q in ("GM", "GM2", "GM3")}
ours = {"GM": GM, "GM2": GM2, "GM3": GM3p}
paper = {"GM": pg1, "GM2": pg2, "GM3": pg3}


def rms_at(q, f):
    m = meas[q][(meas[q][:, 0] >= V_LO) & (meas[q][:, 0] <= 1.0)]
    return np.sqrt(np.mean((f(m[:, 0]) - m[:, 1]) ** 2)), len(m)


def features(v, g2, g3):
    """G_M2 peak, G_M3 peaks, and the G_M3 zero crossing between its two peaks."""
    i2, ia, ib = np.argmax(g2), np.nanargmax(g3), np.nanargmin(g3)
    s = slice(ia, ib + 1)
    vs, gs = v[s], g3[s]
    j = np.where(np.sign(gs[:-1]) != np.sign(gs[1:]))[0][0]
    vz = vs[j] - gs[j] * (vs[j + 1] - vs[j]) / (gs[j + 1] - gs[j])
    return g2[i2], v[i2], g3[ia], v[ia], g3[ib], v[ib], vz


print(f"Fig. 6(c), Device 2, V_DS = {VDS:.0f} V (per mm, W = {W_MM} mm)")
print(f"  definitions from the paper's own lines: G_M2 = dG_M/dV_GS (integral of G_M2 gives G_M to rms"
      f" {rms_gm2_def:.4f} S/mm); G_M3 = central difference of G_M2 over +-H, H = {H:.3f} V"
      f" (rms {min(err):.3f} S/(mm V^2))")
print(f"  derivative of the paper's G_M2 line (+-20 mV) peaks at {exact3.max():.2f} / {exact3.min():.2f};"
      f" printed G_M3 at {np.interp(vchk, *pg3).max():.2f} / {np.interp(vchk, *pg3).min():.2f}")
print()
print(f"rms against the measured points, {V_LO} V <= V_GS <= 1 V")
print("  quantity   unit         this work   paper's model   points")
for q, u in (("GM", "S/mm"), ("GM2", "S/(mm V)"), ("GM3", "S/(mm V^2)")):
    r1, n = rms_at(q, lambda x, q=q: np.interp(x, V, ours[q]))
    r2, _ = rms_at(q, lambda x, q=q: np.interp(x, *paper[q], left=0.0))
    print(f"  {q:<9}  {u:<11}  {r1:9.3f}   {r2:13.3f}   {n:6d}")
print()
fm = features(meas["GM2"][:, 0], meas["GM2"][:, 1], meas["GM3"][:, 1])
vp = np.arange(-1.8, 1.0 + 1e-9, 0.005)
fp = features(vp, np.interp(vp, *pg2, left=0.0), np.interp(vp, *pg3))
fo = features(V, GM2, GM3p)
fx = features(V, GM2, GM3x)
print("features                     measured   paper's model   this work   this work (exact G_M3)")
names = ["G_M2 peak (S/(mm V))", "  at V_GS (V)", "G_M3 max (S/(mm V^2))", "  at V_GS (V)",
         "G_M3 min (S/(mm V^2))", "  at V_GS (V)", "G_M3 zero crossing (V)"]
for i, nm in enumerate(names):
    print(f"  {nm:<26} {fm[i]:8.2f}   {fp[i]:13.2f}   {fo[i]:9.2f}   {fx[i]:9.2f}")
print(f"G_M peak: measured {meas['GM'][:, 1].max():.3f} S/mm at {meas['GM'][np.argmax(meas['GM'][:, 1]), 0]:.2f} V;"
      f" this work {GM.max():.3f} S/mm at {V[np.argmax(GM)]:.2f} V")

# ---- 4. Figure 1: one panel per quantity (measured, paper's model, this work), axes as the paper (-3 to 1 V)
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.3))
units = {"GM": "G_M (S/mm)", "GM2": "G_M2 (S/(mm V))", "GM3": "G_M3 (S/(mm V$^2$))"}
for a, q in zip(ax, ("GM", "GM2", "GM3")):
    a.plot(meas[q][:, 0], meas[q][:, 1], "o", mfc="none", c="#52514e", ms=6, label="measured [P]")
    a.plot(*paper[q], "-", c="#2a78d6", lw=2, label="paper's model [P]")
    a.plot(V, ours[q], "--", c="#eb6834", lw=2, label="this work (Step 8c, not refitted)")
    a.axhline(0, c="#b5b3aa", lw=0.8)
    a.set(xlabel="V_GS (V)", ylabel=units[q], xlim=(-3.0, 1.0))
    a.grid(c="#e1e0d9")
    a.spines[["top", "right"]].set_visible(False)
ax[0].set_title("Device 2, V_DS = 5 V", loc="left", fontsize=10)
ax[1].set_title("G_M2 = dG_M/dV_GS", loc="left", fontsize=10)
ax[2].set_title(f"G_M3: central difference of G_M2, +-{H:.2f} V (as the paper)", loc="left", fontsize=10)
ax[0].legend(frameon=False, fontsize=8, loc="lower right")
plt.tight_layout()
plt.savefig("fig6c_device2.png", dpi=150)

# ---- 5. Figure 2: evidence for the G_M3 definition (paper's own lines only; printed item 14)
fig, a = plt.subplots(figsize=(6.5, 4.3))
a.plot(vchk, exact3, ":", c="#2a78d6", lw=1.5, label="derivative of the paper's G_M2 line (+-20 mV)")
a.plot(*pg3, "-", c="#2a78d6", lw=2.5, label="paper's G_M3 line, as printed [P]")
a.plot(vchk, (P2(vchk + H) - P2(vchk - H)) / (2 * H), "--", c="#14213d", lw=1.2,
       label=f"central difference of the paper's G_M2, +-{H:.2f} V [CALC]")
a.axhline(0, c="#b5b3aa", lw=0.8)
a.set(xlabel="V_GS (V)", ylabel="G_M3 (S/(mm V$^2$))", xlim=(-2.0, 0.5))
a.set_title("Fig. 6(c): the printed G_M3 is a smoothed derivative", loc="left", fontsize=10)
a.grid(c="#e1e0d9"); a.spines[["top", "right"]].set_visible(False)
a.legend(frameon=False, fontsize=8, loc="upper right")
plt.tight_layout()
plt.savefig("fig6c_gm3_definition.png", dpi=150)
print("\nsaved fig6c_device2.png, fig6c_gm3_definition.png")
