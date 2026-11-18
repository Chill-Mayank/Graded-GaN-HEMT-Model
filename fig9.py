"""
fig9.py - Step 11b: Device 2 small-signal RF against Fig. 9 of Nazir et al. 2025 (V_DS = 3 V).
Run:  GAN_DEVICE=2 python3 fig9.py      -> fig9_device2.png, rf_device2.json
Model: Verilog-A twin (rf_twin.py) with the Step 11a DC parameters (params_device2_dibl.json: DIBL from TCAD, D5).
RF-device calibration [CALIBRATE, user choice 1B, 18 Nov; DEVIATION D6]:
  1. DVOFF (V_OFF shift) and U0 (access mobility) to the low-frequency g_m = Re(Y21) of Fig. 9(b) at 1-3 GHz:
     the S-parameter device is on at V_GS = -2 V while Fig. 6(b) is off below -1.5 V (candidate printed item 17).
  2. R_G and constant C_GS,Fr, C_GD,Fr, C_DS,Fr [ASSUME constant] to Fig. 9(a) |S21| and Fig. 9(d) H21, all biases.
  g_ds (Fig. 9(c)) and g_m above 3 GHz are predicted, not fitted.
"""
import json
import numpy as np
from scipy.optimize import least_squares
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rf_twin import RF
from fit2_dibl import BASE, KEY

L = lambda f: np.loadtxt(f, delimiter=",", skiprows=3)
VG, VDS = (-3, -2, -1, 0), 3.0
Q = {"a": "S21_dB", "b": "gm_mS", "c": "gds_mS", "d": "H21_dB"}
M = {p: {v: L(f"paper_data/fig9{p}_measured_{q}_VGS{v}.csv") for v in VG} for p, q in Q.items()}
PM = {p: {v: L(f"paper_data/fig9{p}_paper_model_{q}_VGS{v}.csv") for v in VG} for p, q in Q.items()}
P = json.load(open("params_device2_dibl.json"))["fitted [CALIBRATE]"]
OV = {**{k: BASE[k] for k in ("VOFF", "SIGMA", "VDSREF")}, **{KEY[n]: v for n, v in P.items()}}
FS = np.logspace(9, np.log10(110e9), 60)
low = lambda a: np.median(a[a[:, 0] < 3, 1])


def curves(r):
    out = {}
    for v in VG:
        Y, _ = r.y2(FS, v, VDS)
        S = np.array([np.linalg.solve(np.eye(2) + 50 * y, np.eye(2) - 50 * y) for y in Y])
        H, _ = r.gains(Y)
        out[v] = {"a": 20 * np.log10(abs(S[:, 1, 0])), "b": Y[:, 1, 0].real * 1e3, "c": Y[:, 1, 1].real * 1e3,
                  "d": 20 * np.log10(H)}
    return out


def at(c, p, v, f):
    return np.interp(np.log(f), np.log(FS / 1e9), c[v][p])


# 1. RF-device shift and access mobility from low-frequency g_m
gm_t = np.array([low(M["b"][v]) for v in VG])
def gm_low(x):
    r = RF(2, over={**OV, "DVOFF": x[0], "U0": x[1]})
    return np.array([r.y2([2e9], v, VDS)[0][0, 1, 0].real * 1e3 for v in VG])
s1 = least_squares(lambda x: gm_low(x) - gm_t, [-0.97, 0.1], bounds=([-2, 0.005], [0, 1]))
OV.update(DVOFF=s1.x[0], U0=s1.x[1])
r = RF(2, over=OV)

# 2. R_G and fringing capacitors from S21 and H21
def res(x):
    r.rg, r.cgs, r.cgd, r.cds = x[0], x[1] * 1e-15, x[2] * 1e-15, x[3] * 1e-15
    c = curves(r)
    return np.concatenate([at(c, p, v, M[p][v][:, 0]) - M[p][v][:, 1] for p in "ad" for v in VG if not (p == "d" and v == -3)])
s2 = least_squares(res, [5, 5, 5, 5], bounds=([0] * 4, [200] * 4))
res(s2.x)
c = curves(r)
rf = {"DVOFF_V": s1.x[0], "U0_m2Vs": s1.x[1], "RG_ohm": s2.x[0], "CGSFR_F": s2.x[1] * 1e-15, "CGDFR_F": s2.x[2] * 1e-15,
      "CDSFR_F": s2.x[3] * 1e-15}
json.dump(rf, open("rf_device2.json", "w"), indent=2)
print("RF device: DVOFF %.3f V, U0 %.4f m^2/Vs (DC fit %.4f); R_G %.1f ohm, C_GS,Fr %.1f fF, C_GD,Fr %.1f fF, C_DS,Fr %.1f fF"
      % (s1.x[0], s1.x[1], P["u0"], *s2.x))
print("rms against measured, all biases     this work   paper's model")
for p in "abcd":
    vs = [v for v in VG if not (p == "d" and v == -3)]
    o = np.concatenate([at(c, p, v, M[p][v][:, 0]) - M[p][v][:, 1] for v in vs])
    pp = np.concatenate([np.interp(M[p][v][:, 0], *PM[p][v].T) - M[p][v][:, 1] for v in vs])
    print(f"  Fig. 9({p}) {Q[p]:<8}               {np.sqrt(np.mean(o ** 2)):8.2f}   {np.sqrt(np.mean(pp ** 2)):10.2f}")
print("low-frequency g_ds (mS), V_GS = -3, -2, -1, 0: ours", np.round([at(c, 'c', v, 2.0) for v in VG], 2),
      " measured", np.round([low(M['c'][v]) for v in VG], 2))

COL = {-3: "#a94b17", -2: "#14213d", -1: "#b8860b", 0: "#1d5b3a"}
fig, ax = plt.subplots(1, 4, figsize=(17, 4.2))
lab = {"a": "|S21| (dB)", "b": "g_m = Re(Y21) (mS)", "c": "g_ds = Re(Y22) (mS)", "d": "H21 (dB)"}
for k, p in enumerate("abcd"):
    for v in VG:
        if p == "d" and v == -3:
            continue
        m = M[p][v]
        ax[k].semilogx(m[:, 0], m[:, 1], "o", mfc="none", c=COL[v], ms=3, alpha=0.6)
        ax[k].semilogx(*PM[p][v].T, "-", c=COL[v], lw=1, alpha=0.5)
        ax[k].semilogx(FS / 1e9, c[v][p], "--", c=COL[v], lw=2, label=f"V_GS = {v} V")
    ax[k].set(xlabel="Frequency (GHz)", ylabel=lab[p], xlim=(1, 115))
    ax[k].set_title(f"({p})", loc="left", fontsize=10)
    ax[k].grid(c="#e1e0d9")
ax[0].legend(frameon=False, fontsize=8)
fig.suptitle("Device 2, V_DS = 3 V: measured (circles), paper's model (thin), this work (dashed); "
             "RF-device V_OFF shift and U0 from low-f g_m, R_G and fringing C from S21/H21", fontsize=10)
plt.tight_layout()
plt.savefig("fig9_device2.png", dpi=130)
print("saved fig9_device2.png, rf_device2.json")
