"""
fig11.py - Step 10c: Device 1 small-signal RF against Fig. 11 of Nazir et al. 2025.
Run:  python3 fig11.py      -> fig11_device1.png
Model: graded_hemt.va through its Python twin (rf_twin.py), common source, V_DS = 10 V.
Extrinsic values fitted to Fig. 11(a) only (V_GS = 0 V) [CALIBRATE, user choice A, 18 Nov]:
  R_G = 17.5 ohm, C_GD,Fr = 16.7 fF, C_GS,Fr = 0 (fit floor), C_DS,Fr = 0 (no effect on H21 or U) [ASSUME].
Fig. 11(b) is then predicted with nothing refitted.
Fig. 11(b) extraction: -20 dB/decade from f0 = 2 GHz, as the paper's own lines are (printed item 15) [CALC].
Our U, like the paper's U line, flattens above ~30 GHz (Ward-Dutton drain transcapacitance), so the U = 1 crossing
lies far above the extrapolated f_MAX; both are printed.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rf_twin import RF

RG, CGS, CGD = 17.5, 0.0, 16.7e-15
VDS = 10.0
F0 = 2e9          # Hz, extraction frequency for Fig. 11(b) [ASSUME; the paper's lines give the same result for 1-5 GHz]
L = lambda f: np.loadtxt(f, delimiter=",", skiprows=3)
P = "paper_data/fig11"

r = RF(1, rg=RG, cgs=CGS, cgd=CGD)
f = np.logspace(9, np.log10(60e9), 200)
Y, I0 = r.y2(f, 0.0, VDS)
H, U = r.gains(Y)
rows = np.array([[*r.ft_fmax(v, VDS, mode="extrap", f0=F0)] for v in np.arange(-1.1, 3.01, 0.05)])
rows = np.column_stack([rows[:, 2] * 1e3 / 0.15, rows[:, 0] / 1e9, rows[:, 1] / 1e9])   # mA/mm, GHz, GHz

m = {q: L(f"{P}a_measured_{q}.csv") for q in ("H21_dB", "U_dB")}
pm = {q: L(f"{P}a_paper_model_{q}.csv") for q in ("H21_dB", "U_dB")}
mb = {q: L(f"{P}b_measured_{q}.csv") for q in ("fT_GHz", "fMAX_GHz")}
pb = {q: L(f"{P}b_paper_model_{q}.csv") for q in ("fT_GHz", "fMAX_GHz")}

ours_a = {"H21_dB": lambda x: np.interp(np.log(x * 1e9), np.log(f), 20 * np.log10(H)),
          "U_dB": lambda x: np.interp(np.log(x * 1e9), np.log(f), 10 * np.log10(U))}
ours_b = {"fT_GHz": lambda x: np.interp(x, rows[:, 0], rows[:, 1]),
          "fMAX_GHz": lambda x: np.interp(x, rows[:, 0], rows[:, 2])}
rms = lambda a, b: np.sqrt(np.mean((a - b) ** 2))
print(f"Device 1, V_DS = {VDS} V, V_GS = 0 V: I_D = {I0 * 1e3:.2f} mA ({I0 * 1e3 / 0.15:.0f} mA/mm)")
print("rms against the measured points          this work   paper's model")
for q in ("H21_dB", "U_dB"):
    print(f"  Fig. 11(a) {q:<10} (dB)            {rms(ours_a[q](m[q][:, 0]), m[q][:, 1]):8.2f}   "
          f"{rms(np.interp(m[q][:, 0], *pm[q].T), m[q][:, 1]):10.2f}")
for q in ("fT_GHz", "fMAX_GHz"):
    print(f"  Fig. 11(b) {q:<10} (GHz)           {rms(ours_b[q](mb[q][:, 0]), mb[q][:, 1]):8.2f}   "
          f"{rms(np.interp(mb[q][:, 0], *pb[q].T), mb[q][:, 1]):10.2f}")
ft0, fm0, _ = r.ft_fmax(0.0, VDS, mode="extrap", f0=F0)
ftc, fmc, _ = r.ft_fmax(0.0, VDS)
print(f"  at V_GS = 0 V: extrapolated f_T {ft0 / 1e9:.1f} GHz, f_MAX {fm0 / 1e9:.1f} GHz (measured 22.9, 42.9 near 260 mA/mm);"
      f" U = 1 crossing {fmc / 1e9:.1f} GHz, |H21| = 1 crossing {ftc / 1e9:.1f} GHz")

fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
for q, c, lab in (("H21_dB", "#a94b17", "H21"), ("U_dB", "#2a78d6", "U")):
    ax[0].semilogx(m[q][:, 0], m[q][:, 1], "o", mfc="none", c=c, ms=6, label=f"{lab} measured [P]")
    ax[0].semilogx(*pm[q].T, "-", c=c, lw=1.2, alpha=0.6, label=f"{lab} paper's model [P]")
    ax[0].semilogx(f / 1e9, ours_a[q](f / 1e9), "--", c=c, lw=2, label=f"{lab} this work")
ax[0].set(xlabel="Frequency (GHz)", ylabel="Gain (dB)", xlim=(1, 60), ylim=(0, 45))
ax[0].set_title("(a) V_GS = 0 V, V_DS = 10 V; R_G, C_GD,Fr fitted here", loc="left", fontsize=10)
for q, c, lab in (("fT_GHz", "#a94b17", "f_T"), ("fMAX_GHz", "#2a78d6", "f_MAX")):
    ax[1].plot(mb[q][:, 0], mb[q][:, 1], "o", mfc="none", c=c, ms=6, label=f"{lab} measured [P]")
    ax[1].plot(*pb[q].T, "-", c=c, lw=1.2, alpha=0.6, label=f"{lab} paper's model [P]")
    ax[1].plot(rows[:, 0], rows[:, 1 if q == "fT_GHz" else 2], "--", c=c, lw=2, label=f"{lab} this work")
ax[1].set(xlabel="I_DS (mA/mm)", ylabel="f_T, f_MAX (GHz)", xlim=(0, 620), ylim=(0, 60))
ax[1].set_title("(b) V_DS = 10 V; predicted, not refitted (-20 dB/dec from 2 GHz)", loc="left", fontsize=10)
for a in ax:
    a.grid(c="#e1e0d9"); a.spines[["top", "right"]].set_visible(False)
ax[0].legend(frameon=False, fontsize=7.5)
ax[1].legend(frameon=False, fontsize=7.5, ncol=2, loc="lower center")
plt.tight_layout()
plt.savefig("fig11_device1.png", dpi=150)
print("saved fig11_device1.png")
