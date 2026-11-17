"""
docs/make_figures.py - every figure shown in README.md, regenerated from this repository.
Run from anywhere:   python3 docs/make_figures.py        (about 2-3 min; writes docs/figures/*.png, prints the quoted numbers)

Sources of each curve
  - the compact model in this repo (Device 1: params_device1.json; Device 2: params_device2.json, run in a
    subprocess with GAN_DEVICE=2 because device.py picks the device at import),
  - the paper's plotted data, read from the PDF's vector drawing (paper_data/*.csv),
  - our Sentaurus TCAD runs (docs/data/dev2_*: cuts and sweeps of runs n4-n8),
  - the Step 5d closed-form comparison (docs/data/closed_vs_numeric.npz: Device 1, closed form Eqs. (9)-(15) and the
    numerical solution, each refitted to the measured data).
Style: measured = open circles, the paper's model = solid blue, this work = dashed orange, our TCAD = aqua.
"""
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT, DATA = ROOT / "docs" / "figures", ROOT / "docs" / "data"

# ---------------------------------------------------------------- Device 2 curves (child process, GAN_DEVICE=2)
if len(sys.argv) == 3 and sys.argv[1] == "--device2":
    sys.path.insert(0, str(ROOT)); os.chdir(ROOT)
    import numpy as np
    import fit2 as f2
    p = list(json.load(open("params_device2.json"))["fitted [CALIBRATE]"].values())
    V = np.linspace(-1.8, 1.0, 141)
    I = np.array([f2.model_id(p, v, 5.0) for v in V])
    G = np.array([f2.model_gm(p, v) for v in V])
    np.savez(sys.argv[2], V=V, I=I, G=G)
    sys.exit(0)

os.environ["GAN_DEVICE"] = "1"                     # everything below is Device 1 (device.py default)
sys.path.insert(0, str(ROOT)); os.chdir(ROOT)
OUT.mkdir(parents=True, exist_ok=True)
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from charge import sheet_exact, n_sub, n_mid, n_on, n_on_printed
import gate as g
from access import id_ext
import terminal_charges as tc
import capacitances as cp

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRID, SURF, BAND = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#ffffff", "#f1efe8"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 10.5, "axes.labelsize": 10,
                     "legend.fontsize": 8.8, "xtick.labelsize": 9, "ytick.labelsize": 9, "axes.edgecolor": "#8f8e88",
                     "axes.labelcolor": INK, "xtick.color": "#8f8e88", "ytick.color": "#8f8e88", "xtick.labelcolor": INK2,
                     "ytick.labelcolor": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False, "figure.facecolor": SURF, "axes.facecolor": SURF, "lines.linewidth": 1.7,
                     "axes.linewidth": 0.8, "savefig.dpi": 200, "mathtext.default": "regular"})
MEAS = dict(ls="none", marker="o", mfc="none", mec=INK2, ms=4.6, mew=1.0)


def load(name):
    return np.loadtxt(f"paper_data/{name}.csv", delimiter=",", skiprows=3)


def loadd(name):
    return np.loadtxt(DATA / name, delimiter=",", skiprows=3)


def grid(a):
    a.grid(True, color=GRID, lw=0.6)
    a.set_axisbelow(True)


def ttl(a, t):
    a.set_title(t, loc="left", color=INK, pad=6)


def save(fig, name):
    fig.savefig(OUT / f"{name}.png")
    plt.close(fig)
    print("saved docs/figures/" + name + ".png")


rms = lambda x: float(np.sqrt(np.mean(np.asarray(x) ** 2)))

# ================================================================ Device 1 model (shared by several figures)
params = json.load(open("params_device1.json"))
P1 = {**params["fixed"], **params["fitted [CALIBRATE]"]}
for k in ("n_acc_cm2", "l_sg", "l_gd"):
    P1.pop(k, None)
model1 = lambda vg, vd: id_ext(vg, vd, **P1) * 1e3          # mA, Device 1 (W = 150 um)
idm, idp, om, op = (load(n) for n in ["fig6a_measured_ID", "fig6a_paper_model_ID", "fig7a_measured", "fig7a_paper_model"])
gmm, gmp = load("fig6a_measured_gm"), load("fig6a_paper_model_gm")
Vg1 = np.linspace(-3, 3, 601)
I1 = np.array([model1(v, 10.0) for v in Vg1])
gm1 = np.gradient(I1, Vg1)

# ================================================================ hero: both devices, per mm of gate width
with tempfile.TemporaryDirectory() as tmp:
    f = Path(tmp) / "dev2.npz"
    subprocess.run([sys.executable, str(Path(__file__).resolve()), "--device2", str(f)], check=True,
                   env={**os.environ, "GAN_DEVICE": "2"})
    d2 = dict(np.load(f))
W1, W2 = 0.150, 0.075                                       # mm [P]
i2m, i2p = load("fig6b_measured_ID"), load("fig6b_paper_model_ID")
g2m, g2p = load("fig6b_measured_gm"), load("fig6b_paper_model_gm")   # already per mm (printed item 11)
fig, ax = plt.subplots(2, 2, figsize=(9.2, 5.6), constrained_layout=True, sharex="col")
for col, (title, vlim) in enumerate([(r"Device 1: long-channel PolFET, $L_G$ = 0.7 µm, $V_{DS}$ = 10 V", (-3, 3)),
                                     (r"Device 2: short-channel PolFET, $L_G$ = 50 nm, $V_{DS}$ = 5 V", (-2, 1))]):
    a, b = ax[0, col], ax[1, col]
    if col == 0:
        a.plot(idm[:, 0], idm[:, 1] / W1, **MEAS, label="measured [paper]")
        a.plot(idp[:, 0], idp[:, 1] / W1, color=BLUE, label="paper's model")
        a.plot(Vg1, I1 / W1, color=ORANGE, ls="--", label="this work")
        b.plot(gmm[:, 0], gmm[:, 1] / W1, **MEAS); b.plot(gmp[:, 0], gmp[:, 1] / W1, color=BLUE)
        b.plot(Vg1, gm1 / W1, color=ORANGE, ls="--")
    else:
        a.plot(i2m[:, 0], i2m[:, 1] / W2, **MEAS); a.plot(i2p[:, 0], i2p[:, 1] / W2, color=BLUE)
        a.plot(d2["V"], d2["I"] / W2, color=ORANGE, ls="--")
        b.plot(g2m[:, 0], g2m[:, 1], **MEAS); b.plot(g2p[:, 0], g2p[:, 1], color=BLUE)
        b.plot(d2["V"], d2["G"] / W2, color=ORANGE, ls="--")
    a.set(xlim=vlim, ylabel=r"$I_D$ (mA/mm)"); b.set(xlim=vlim, ylim=(0, None), xlabel=r"$V_{GS}$ (V)",
                                                    ylabel=r"$g_m$ (mS/mm)")
    a.set_ylim(0, None); ttl(a, title)
    grid(a); grid(b)
ax[0, 0].legend(loc="upper left")
save(fig, "hero_two_devices")

# ================================================================ charge model against the paper's Fig. 4
p4 = {k: load(f"fig4a_{k}") for k in ["numerical", "off", "middle", "deep_on"]}
pb = load("fig4b_markers")
eta = np.linspace(-4, 10, 701); ep = np.linspace(-4, 10, 22); nx = np.array([sheet_exact(e) for e in ep])
fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.9), constrained_layout=True)
a = ax[0]
a.semilogy(eta, n_sub(eta), color=ORANGE, lw=1.3, label="OFF, Eq. (9)")
a.semilogy(eta, n_mid(eta), color=AQUA, ls="-.", lw=1.3, label="Middle, Eq. (13)")
a.semilogy(eta, n_on(eta), color=BLUE, label="Deep ON, Eq. (10) corrected")
a.semilogy(eta, n_on_printed(eta), color=BLUE, ls="--", lw=1.3, label="Deep ON, Eq. (10) as printed")
a.semilogy(ep, nx, **{**MEAS, "ms": 4.2}, label=r"exact $(N_C F_{1/2})^{2/3}$")
for k in ["off", "middle", "deep_on"]:
    a.semilogy(p4[k][::8, 0], p4[k][::8, 1], ".", color=INK, ms=2.2)
a.semilogy(p4["numerical"][:, 0], p4["numerical"][:, 1], ".", color=INK, ms=2.2, label="paper, Fig. 4(a)")
a.set(xlim=(-4, 10), ylim=(1e11, 3e15), xlabel=r"$\eta_F$", ylabel=r"n (cm$^{-2}$)"); grid(a)
a.legend(loc="upper left", handlelength=2.2, fontsize=8.3, labelspacing=0.35)
ttl(a, "(a) Regional solutions, paper Fig. 4(a)")
V = np.linspace(-5, 0, 501); nex = np.array([g.n_exact(v)[0] for v in V]); nt = g.n_total(V); ok = nt > 0
a = ax[1]
a.semilogy(V, nex, color=MUTED, lw=1.4, label="numerical solution")
a.semilogy(V[ok], nt[ok], color=BLUE, label="corrected closed form")
a.semilogy(pb[:, 0], pb[:, 1], **MEAS, label="paper, Fig. 4(b)")
a.axvline(g.V_OFF, color=GRID, lw=1.0)
a.text(g.V_OFF - 0.08, 1.2e13, r"$V_{OFF}$ = −3 V", color=INK2, fontsize=8.8, ha="right", va="center")
a.set(xlim=(-5, 0), ylim=(1e-1, 1e14), xlabel=r"$V_{GS}$ (V)", ylabel=r"$n_{total}$ (cm$^{-2}$)"); grid(a)
a.legend(loc="lower right", bbox_to_anchor=(1.0, 0.08))
ttl(a, "(b) Charge against gate voltage, paper Fig. 4(b)")
save(fig, "charge_model")
sel = (pb[:, 0] >= -2.5) & (pb[:, 0] <= 0)
rr = np.array([g.n_exact(v)[0] for v in pb[sel, 0]]) / pb[sel, 1] - 1
print(f"  Fig. 4(b): numerical solution vs paper points, -2.5..0 V: rms {100 * rms(rr):.2f} % ({sel.sum()} points)")

# ================================================================ closed form against the numerical solution (Device 1)
d = np.load(DATA / "closed_vs_numeric.npz")
fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.8), constrained_layout=True)
a = ax[0]; vnr = float(d["v_noroot"])
a.axvspan(-1.8, vnr, color=ORANGE, alpha=0.10, lw=0)
a.semilogy(d["V"], d["nE"], color=INK2, lw=2.0, label="numerical solution")
a.semilogy(d["V"], d["nA"], color=ORANGE, label="closed form, Eqs. (9)–(15)")
a.text(-1.74, 1.35e13, "Eq. (14): no real root\nfor $V_{GS}$ < " + f"{vnr:.2f} V".replace("-", "−"), color=INK2,
       fontsize=8.8, va="top")
a.set(xlim=(-1.8, 3), ylim=(1e9, 2e13), xlabel=r"$V_{GS}$ (V), intrinsic", ylabel=r"n (cm$^{-2}$)"); grid(a)
a.legend(loc="lower right")
ttl(a, "(a) Device 1 channel charge")
a = ax[1]
a.plot(gmm[:, 0], gmm[:, 1], **MEAS, label="measured")
a.plot(d["Vg"], d["gmE"], color=INK2, lw=2.0, label="numerical solution, refit")
a.plot(d["Vg"], d["gmA"], color=ORANGE, label="closed form, refit")
nv, ng = float(d["notch_v"]), float(d["notch_gm"])
a.annotate(f"notch at {nv:.2f} V".replace("-", "−"), (nv, ng), (0.25, 21.5), color=INK2, fontsize=8.8,
           arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
a.set(xlim=(-1.7, 3), ylim=(0, 50), xlabel=r"$V_{GS}$ (V)", ylabel=r"$g_m$ (mS)"); grid(a); a.legend(loc="lower right")
ttl(a, r"(b) $g_m$ at $V_{DS}$ = 10 V")
save(fig, "closed_vs_numerical")

# ================================================================ Device 1 against the measurement and the paper's model
fam = np.arange(-1.0, 3.01, 0.5); VD = np.linspace(0, 10, 101)
Iout = {vg: np.array([model1(vg, vd) for vd in VD]) for vg in fam}
fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.4), constrained_layout=True)
a = ax[0]
a.plot(idm[:, 0], idm[:, 1], **MEAS, label="measured"); a.plot(idp[:, 0], idp[:, 1], color=BLUE, label="paper's model")
a.plot(Vg1, I1, color=ORANGE, ls="--", label="this work")
a.set(xlim=(-3, 3), ylim=(0, 160), xlabel=r"$V_{GS}$ (V)", ylabel=r"$I_D$ (mA)"); a.legend(loc="upper left")
ttl(a, r"(a) Transfer, $V_{DS}$ = 10 V")
a = ax[1]
a.plot(gmm[:, 0], gmm[:, 1], **MEAS); a.plot(gmp[:, 0], gmp[:, 1], color=BLUE); a.plot(Vg1, gm1, color=ORANGE, ls="--")
a.set(xlim=(-3, 3), ylim=(0, 50), xlabel=r"$V_{GS}$ (V)", ylabel=r"$g_m$ (mS)"); ttl(a, r"(b) $g_m$, $V_{DS}$ = 10 V")
a = ax[2]
for vg in fam:
    s = om[om[:, 0] == vg]; a.plot(s[:, 1], s[:, 2], **{**MEAS, "ms": 3.4, "mew": 0.8})
    s = op[op[:, 0] == vg]; a.plot(s[:, 1], s[:, 2], color=BLUE, lw=1.2)
    a.plot(VD, Iout[vg], color=ORANGE, ls="--", lw=1.2)
a.text(10.25, Iout[3.0][-1], "3 V", va="center", color=INK2, fontsize=8.6)
a.text(10.25, Iout[-1.0][-1], "−1 V", va="center", color=INK2, fontsize=8.6)
a.set(xlim=(0, 11.9), ylim=(0, 160), xlabel=r"$V_{DS}$ (V)", ylabel=r"$I_D$ (mA)"); ttl(a, "(c) Output")
for a in ax:
    grid(a)
save(fig, "device1_dc")
id6 = idm[idm[:, 0] > -1.4]; gm6 = gmm[(gmm[:, 0] > -1.4) & (gmm[:, 0] < 2.95)]; out7 = om[om[:, 1] > 0]
h = 0.01
print(f"  Device 1 vs measured: I_D {rms([model1(v, 10.0) - i for v, i in id6]):.3f} mA, "
      f"g_m {rms([(model1(v + h, 10.0) - model1(v - h, 10.0)) / (2 * h) - y for v, y in gm6]):.3f} mS, "
      f"output {rms([model1(v, vd) - i for v, vd, i in out7]):.3f} mA")

# ================================================================ gate capacitances against drain bias (Device 1)
kw = tc.model_params(); fF = 1e15
VDs = np.logspace(np.log10(4e-4), 1, 70)
res = {}
for vgs in (0.0, -1.0):
    C = np.array([cp.cap_matrix(vgs, vd, h=min(1e-3, vd / 4), **kw) for vd in VDs]) * fF
    a_ = cp.cgg_at_vds0(vgs, **kw) / tc.QG_SCALE
    res[vgs] = dict(cgs=C[:, 0, 2], cgd=C[:, 0, 1], a=a_, lim=1 / (2 * a_ - 1), cgg0=cp.cgg_at_vds0(vgs, **kw) * fF)
    print(f"  V_GS {vgs:+.1f}: a = {a_:.4f}, limit 1/(2a-1) = {1 / (2 * a_ - 1):.4f}")
fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.7), constrained_layout=True)
a = ax[0]; r = res[0.0]
a.semilogx(VDs, r["cgs"], color=BLUE, label=r"$C_{gs}$"); a.semilogx(VDs, r["cgd"], color=ORANGE, ls="--", label=r"$C_{gd}$")
a.axhline(r["cgg0"] / 2, color=MUTED, lw=0.9, ls=":")
a.annotate("symmetric model:\nboth → " + r"$C_{gg}$/2 = " + f"{r['cgg0'] / 2:.0f} fF", (3e-2, r["cgg0"] / 2),
           (1.2e-1, 302), color=INK2, fontsize=8.6, va="center", ha="right",
           arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
a.text(4.5e-4, r["cgd"][0] + 9, r"$C_{gd}$ → " + f"{r['cgd'][0]:.0f} fF", color=INK2, fontsize=8.6)
a.text(4.5e-4, r["cgs"][0] - 25, r"$C_{gs}$ → " + f"{r['cgs'][0]:.0f} fF", color=INK2, fontsize=8.6)
a.set(xlim=(4e-4, 10), ylim=(0, 430), xlabel=r"$V_{DS}$ (V)", ylabel="capacitance (fF)"); grid(a)
a.legend(loc="center left", bbox_to_anchor=(0.0, 0.22))
ttl(a, r"(a) $C_{gs}$, $C_{gd}$ at $V_{GS}$ = 0 V (Device 1, intrinsic)")
a = ax[1]
for vgs, col, ls in ((0.0, BLUE, "-"), (-1.0, ORANGE, "--")):
    r = res[vgs]; m = VDs <= 0.1
    a.semilogx(VDs[m], (r["cgd"] / r["cgs"])[m], color=col, ls=ls, label=r"$V_{GS}$ = " + f"{vgs:.0f} V".replace("-", "−"))
    a.axhline(r["lim"], color=col, lw=0.8, ls=":")
    a.text(4.6e-4, r["lim"] + 0.011, f"limit 1/(2a − 1) = {r['lim']:.3f}", color=INK2, fontsize=8.6, ha="left")
a.axhline(1.0, color=MUTED, lw=0.9); a.text(4.6e-4, 1.011, "symmetric model: 1", color=INK2, fontsize=8.6, ha="left")
a.set(xlim=(4e-4, 0.1), ylim=(0.95, 1.40), xlabel=r"$V_{DS}$ (V)", ylabel=r"$C_{gd}$ / $C_{gs}$"); grid(a)
a.legend(loc="lower left", bbox_to_anchor=(0.0, 0.14))
ttl(a, r"(b) Ratio as $V_{DS}$ → 0")
save(fig, "capacitances")

# ================================================================ Device 2 TCAD: strain reference (runs n4, n5)
n4, n5 = np.load(DATA / "dev2_n4.npz"), np.load(DATA / "dev2_n5.npz")
dep = lambda x: (x + 0.955) * 1e3                           # nm below the barrier top (gate at x = -0.955 um)
fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.8), constrained_layout=True)
layers = [(0, 9, "barrier"), (9, 15, "grade"), (15, 55, "GaN"), (55, 110, r"Al$_{0.04}$GaN buffer")]
for a in ax:
    for k, (x0, x1, lab) in enumerate(layers):
        if k % 2 == 0:
            a.axvspan(x0, x1, color=BAND, lw=0)
a = ax[0]
a.plot(dep(n4["gate_x_ec"]), n4["gate_ec"], color=BLUE, label="GaN reference (default), run n4")
a.plot(dep(n5["gate_x_ec"]), n5["gate_ec"], color=ORANGE, ls="--", label="buffer reference, run n5")
a.axhline(0, color=MUTED, lw=0.8); a.text(108, 0.05, r"$E_F$", color=INK2, fontsize=8.6, ha="right")
for x0, x1, lab in layers:
    a.text((x0 + x1) / 2 if x1 < 100 else 82, 2.03, lab, color=INK2, fontsize=8.2, ha="center",
           rotation=90 if x1 <= 15 else 0, va="top")
a.set(xlim=(0, 110), ylim=(-0.2, 2.05), xlabel="depth below the gate (nm)", ylabel=r"$E_C − E_F$ (eV)"); grid(a)
a.legend(loc="upper right", bbox_to_anchor=(1.0, 0.9))
ttl(a, "(a) Conduction band under the gate")
a = ax[1]
a.plot(dep(n4["gate_x_n"]), n4["gate_n"] / 1e18, color=BLUE); a.plot(dep(n5["gate_x_n"]), n5["gate_n"] / 1e18, color=ORANGE, ls="--")
a.text(24, 6.3, r"run n4: $n_s$ = " + f"{n4['gate_ns'] / 1e12:.2f}" + r"×10$^{12}$ cm$^{-2}$", color=INK2, fontsize=8.8)
a.text(24, 4.9, r"run n5: $n_s$ = " + f"{n5['gate_ns'] / 1e12:.2f}" + r"×10$^{12}$ cm$^{-2}$", color=INK2, fontsize=8.8)
for x0, x1, lab in layers[:3]:
    a.text((x0 + x1) / 2, 8.12, lab, color=INK2, fontsize=8.2, ha="center", rotation=90 if x1 <= 15 else 0, va="top")
a.set(xlim=(0, 60), ylim=(0, 8.2), xlabel="depth below the gate (nm)", ylabel=r"n (10$^{18}$ cm$^{-3}$)"); grid(a)
ttl(a, r"(b) Electron density, $V_{GS}$ = $V_{DS}$ = 0")
save(fig, "device2_strain_reference")

# ================================================================ Device 2 TCAD against the measured Fig. 6(b) (runs n6, n7)
t6, t7 = loadd("dev2_n6_IdVg.csv"), loadd("dev2_n7_IdVg.csv")
fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.7), constrained_layout=True)
a = ax[0]
a.plot(i2m[:, 0], i2m[:, 1] / W2, **MEAS, label=r"measured, $V_{DS}$ = 5 V [paper]")
a.plot(t7[:, 0], t7[:, 2] * 1e3, color=AQUA, label=r"our TCAD, $V_{DS}$ = 5 V (run n7)")
a.set(xlim=(-2, 1), ylim=(0, 900), xlabel=r"$V_{GS}$ (V)", ylabel=r"$I_D$ (mA/mm)"); grid(a); a.legend(loc="upper left")
ttl(a, "(a) Gate work function 5.45 eV, nothing else tuned")
a = ax[1]
pos = i2m[:, 1] / W2 > 0.1                                  # above the plotted zero line (about 0.01 mA/mm)
a.semilogy(i2m[pos, 0], i2m[pos, 1] / W2 * 1e-3, **MEAS, label="measured (5 V)")
a.semilogy(t7[:, 0], np.abs(t7[:, 2]), color=AQUA, label="TCAD n7 (5 V)")
a.semilogy(t6[:, 0], np.abs(t6[:, 2]), color=AQUA, ls=":", label="TCAD n6 (0.05 V)")
a.axhline(1e-3, color=MUTED, lw=0.8, ls="--"); a.text(0.95, 1.4e-3, "1 mA/mm", color=INK2, fontsize=8.6, ha="right")
a.set(xlim=(-3, 1), ylim=(1e-12, 3), xlabel=r"$V_{GS}$ (V)", ylabel=r"$I_D$ (A/mm)"); grid(a)
a.legend(loc="lower right")
ttl(a, "(b) Log scale; measured shown above 0.1 mA/mm")
save(fig, "device2_tcad_vs_measured")
v1 = lambda t: np.interp(-3.0, np.log10(np.abs(t[:, 2]) + 1e-30), t[:, 0])
print(f"  TCAD V_GS at 1 mA/mm: n7 (5 V) {v1(t7):.3f} V, n6 (0.05 V) {v1(t6):.3f} V; "
      f"I_D(0 V, 5 V) TCAD {np.interp(0, t7[:, 0], t7[:, 2]) * 1e3:.0f} vs measured {np.interp(0, i2m[:, 0], i2m[:, 1]) / W2:.0f} mA/mm")

# ================================================================ Device 2 TCAD: C_g'' from n_s(V_GS) (run n8)
ns = loadd("dev2_n8_ns.csv")
EPS = 9.5 * 8.8541878128e-14                                 # F/cm, eps_r = 9.5 [ASSUME, as in the model]
Q = 1.602176634e-19
u_print, u_dt = EPS / 24e-7 / Q, EPS / 15e-7 / Q             # eps/(2d+t), eps/(d+t): d = 9 nm, t = 6 nm [P]
res_fit = least_squares(lambda p: np.log([g.n_exact(v, p[0], p[1])[0] for v in ns[:, 0]]) - np.log(ns[:, 1]),
                        [-0.8, 3e12], x_scale=[0.3, 1e12])
voff, ufit = res_fit.x
vplot = np.linspace(-1.1, 0.5, 161)
fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.7), constrained_layout=True)
a = ax[0]
a.plot(ns[:, 0], ns[:, 1] / 1e12, **{**MEAS, "mec": AQUA, "ms": 6, "mew": 1.4}, label="our TCAD (run n8)")
for u, col, ls, lab in ((ufit, ORANGE, "--", f"charge equation, fitted: $C_g''/q$ = {ufit / 1e12:.2f}"),
                        (u_dt, INK2, "-", f"ε/(d + t) = 15 nm: {u_dt / 1e12:.2f}"),
                        (u_print, BLUE, "-", f"ε/(2d + t) = 24 nm (printed): {u_print / 1e12:.2f}")):
    vo = least_squares(lambda p: np.log([g.n_exact(v, p[0], u)[0] for v in ns[:, 0]]) - np.log(ns[:, 1]), [-0.7]).x[0] \
        if u is not ufit else voff
    a.plot(vplot, [g.n_exact(v, vo, u)[0] / 1e12 for v in vplot], color=col, ls=ls, label=lab)
a.set(xlim=(-1.1, 0.5), ylim=(0, 4.6), xlabel=r"$V_{GS}$ (V), $V_{DS}$ = 0.05 V", ylabel=r"$n_s$ (10$^{12}$ cm$^{-2}$)")
grid(a); a.legend(loc="upper left", fontsize=8.0, title=r"$C_g''/q$ in 10$^{12}$ cm$^{-2}$/V", title_fontsize=8.0)
ttl(a, "(a) Sheet density under the gate (each line: its best V$_{OFF}$)")
a = ax[1]
a.plot(ns[:, 0], ns[:, 2], **{**MEAS, "mec": AQUA, "ms": 6, "mew": 1.4}, label="charge centroid, TCAD")
a.axhline(15, color=INK2, lw=1.2); a.text(-0.6, 15.35, "d + t = 15 nm (GaN / grade interface)", color=INK2, fontsize=8.6, va="top")
a.axhline(24, color=BLUE, lw=1.2); a.text(-1.08, 23.65, r"2d + t = 24 nm, the printed $C_g''$", color=BLUE, fontsize=8.6, va="bottom")
a.axhspan(9, 15, color=BAND, lw=0); a.text(0.48, 10.2, "graded layer (9-15 nm)", color=INK2, fontsize=8.6, ha="right", va="center")
a.set(xlim=(-1.1, 0.5), ylim=(8, 27), xlabel=r"$V_{GS}$ (V)", ylabel="depth below the gate (nm)"); grid(a)
a.invert_yaxis(); a.legend(loc="lower right")
ttl(a, "(b) Where the channel charge sits")
save(fig, "device2_gate_capacitance")
s = (ns[5, 1] - ns[4, 1]) / (ns[5, 0] - ns[4, 0])
print(f"  n8: fitted V_OFF {voff:+.3f} V, C_g''/q {ufit:.3e} (eps/(d+t) {u_dt:.3e}, printed {u_print:.3e}); "
      f"slope -0.4..0 V gives d_eff = {EPS / (Q * s) * 1e7:.1f} nm")
print("done")
