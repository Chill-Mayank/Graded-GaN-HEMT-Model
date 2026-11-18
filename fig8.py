"""
fig8.py - Step 12b: Device 2 S-parameters against Fig. 8 of Nazir et al. 2025 (V_DS = 3 V, 1-110 GHz).
Run:  GAN_DEVICE=2 python3 fig8.py      -> fig8_device2.png
Same model and calibration as fig9.py (Step 11: DIBL D5, RF-device shift D6, R_G and fringing C from Fig. 9);
nothing is refitted here. The extracted points carry no frequency, so the comparison metric is the distance in the
complex plane from each measured point to the nearest point of a model curve (ours, and the paper's line) [CALC].
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rf_twin import RF
from fit2_dibl import BASE, KEY

VG, VDS = (-3, -2, -1, 0), 3.0
COL = {-3: "#a94b17", -2: "#14213d", -1: "#b8860b", 0: "#1d5b3a"}
RIM = {"S11": 1.0, "S12": 0.6, "S21": 3.0, "S22": 1.0}
IDX = {"S11": (0, 0), "S12": (0, 1), "S21": (1, 0), "S22": (1, 1)}
L = lambda f: np.loadtxt(f, delimiter=",", skiprows=3, ndmin=2)
cx = lambda a: a[:, 0] + 1j * a[:, 1]

P = json.load(open("params_device2_dibl.json"))["fitted [CALIBRATE]"]
RFD = json.load(open("rf_device2.json"))
OV = {**{k: BASE[k] for k in ("VOFF", "SIGMA", "VDSREF")}, **{KEY[n]: v for n, v in P.items()},
      "DVOFF": RFD["DVOFF_V"], "U0": RFD["U0_m2Vs"]}
r = RF(2, rg=RFD["RG_ohm"], cgs=RFD["CGSFR_F"], cgd=RFD["CGDFR_F"], cds=RFD["CDSFR_F"], over=OV)
FS = np.logspace(9, np.log10(110e9), 300)
ours = {}
for v in VG:
    Y, _ = r.y2(FS, v, VDS)
    ours[v] = np.array([np.linalg.solve(np.eye(2) + 50 * y, np.eye(2) - 50 * y) for y in Y])

dist = lambda pts, curve: np.array([np.min(np.abs(curve - p)) for p in pts])
print("mean distance of the measured points to the curve (|S| units; rim in brackets)")
print("            this work   paper's model")
for q, (i, j) in IDX.items():
    do, dp = [], []
    for v in VG:
        tag = f"VGS{v}" if v else "VGS0"
        m = cx(L(f"paper_data/fig8_measured_{q}_{tag}.csv"))
        do += list(dist(m, ours[v][:, i, j]))
        dp += list(dist(m, cx(L(f"paper_data/fig8_paper_model_{q}_{tag}.csv"))))
    print(f"  {q} ({RIM[q]:.1f})   {np.mean(do):8.3f}   {np.mean(dp):10.3f}")

fig, ax = plt.subplots(1, 4, figsize=(17, 4.6))
for a, (q, (i, j)) in zip(ax, IDX.items()):
    R = RIM[q]
    t = np.linspace(0, 2 * np.pi, 400)
    a.plot(R * np.cos(t), R * np.sin(t), c="#222", lw=1)
    if q in ("S11", "S22"):
        for rr in (0.2, 0.5, 1.0, 2.0, 5.0):                    # constant-resistance circles of a Smith chart
            c, rad = rr / (1 + rr), 1 / (1 + rr)
            a.plot(c + rad * np.cos(t), rad * np.sin(t), c="#ccc", lw=0.6)
        a.axhline(0, c="#ccc", lw=0.6)
    else:
        for f in (1 / 3, 2 / 3):
            a.plot(f * R * np.cos(t), f * R * np.sin(t), c="#ccc", lw=0.6)
        for ang in range(0, 360, 30):
            a.plot([0, R * np.cos(np.radians(ang))], [0, R * np.sin(np.radians(ang))], c="#ddd", lw=0.5)
    for v in VG:
        tag = f"VGS{v}" if v else "VGS0"
        m = L(f"paper_data/fig8_measured_{q}_{tag}.csv"); pm = L(f"paper_data/fig8_paper_model_{q}_{tag}.csv")
        a.plot(m[:, 0], m[:, 1], "o", mfc="none", c=COL[v], ms=5)
        a.plot(pm[:, 0], pm[:, 1], "-", c=COL[v], lw=1, alpha=0.5)
        s = ours[v][:, i, j]
        a.plot(s.real, s.imag, "--", c=COL[v], lw=2, label=f"V_GS = {v} V")
    a.set_aspect("equal"); a.axis("off")
    a.set_title(f"{q} (rim {R:g})", fontsize=10)
h, l = ax[0].get_legend_handles_labels()
fig.legend(h, l, frameon=False, fontsize=9, loc="lower center", ncol=4)
fig.suptitle("Fig. 8, Device 2, V_DS = 3 V, 1-110 GHz: measured (circles), paper's model (thin), this work (dashed); "
             "not refitted (Fig. 9 calibration)", fontsize=10, y=0.99)
plt.tight_layout(rect=(0, 0.06, 1, 0.93))
plt.savefig("fig8_device2.png", dpi=130)
print("saved fig8_device2.png")
