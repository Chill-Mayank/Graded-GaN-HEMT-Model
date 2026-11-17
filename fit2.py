"""
fit2.py - Step 8b: fit Device 2 (short-channel PolFET) to the measured points of the paper's Fig. 6(b).
Run:  GAN_DEVICE=2 python3 fit2.py        (device2.py supplies W, L, C_g'', V_OFF, mobility, access geometry)
Same method and same five fitted parameters as Device 1 (fit.py, Step 5b); the paper gives no values for either device.
Fitted data: Fig. 6(b) I_D and G_M at V_DS = 5 V [P]. Fig. 7(b) is used as a check, not fitted, because its V_GS
labels are not printed (printed item 12): it is scored with both label sets of fig7b_measured.csv.
Saves params_device2.json.
"""
import json
import os
import time
import numpy as np
from scipy.optimize import least_squares

assert os.environ.get("GAN_DEVICE") == "2", "run as: GAN_DEVICE=2 python3 fit2.py"
import device2 as d2
import current as cur
from access import id_ext, r_access, U0_ACC, N_ACC, L_SG, L_GD

VDS_FIG6B = 5.0                                   # V [P, Fig. 6(b)]
W_MM = d2.W * 1e3                                 # mm, 0.075


def load(name):
    """Read one paper_data CSV (2 comment lines + 1 header line)."""
    return np.loadtxt(f"paper_data/{name}.csv", delimiter=",", skiprows=3)


id6 = load("fig6b_measured_ID")
id6 = id6[id6[:, 0] > -1.5]                      # points on the zero line carry no information (as Device 1)
gm6 = load("fig6b_measured_gm")
gm6 = gm6[gm6[:, 0] > -1.5]
gm6[:, 1] *= W_MM                                # mS/mm -> mS for the 75 um device (printed item 11)
out7 = load("fig7b_measured")
out7 = out7[out7[:, 2] > 0]                      # V_DS = 0 gives I_D = 0 for any parameters

FREE = [("u0",    U0_ACC,         100e-4, 3000e-4, 0.05, "m^2/Vs", "U0, access mobility [T Eq. 5.2]"),
        ("vsat0", cur.VSAT0,      2e4,    1e6,     1e5,  "m/s",    "V_SAT0 [P]"),
        ("th2",   cur.THETA_SAT2, 0.01,   2.0,     0.3,  "-",      "theta_SAT2 [P]"),
        ("delta", cur.DELTA_DS,   1.0,    10.0,    1.0,  "-",      "delta in V_DS,eff [P, no value]"),
        ("beta",  2.0,            0.5,    20.0,    1.0,  "-",      "beta in access R [T Eq. 5.1, no value]")]
NAMES = [f[0] for f in FREE]
FIXED = {"vsat_acc": 1.8e5,     # m/s [TCAD, GaN.par saturation velocity], as Device 1
         "mu_a": 0.0, "theta_sat": 0.0}          # as Device 1 (Step 5b)
W_ID, W_GM = 1.0, 2.0                            # per mA and per mS, as Device 1


def model_id(p, V_GS, V_DS):
    """Model drain current in mA."""
    return id_ext(V_GS, V_DS, **FIXED, **dict(zip(NAMES, p))) * 1e3


def model_gm(p, V_GS, V_DS=VDS_FIG6B, h=0.01):
    """Model G_M in mS (whole device), central difference of +-10 mV."""
    return (model_id(p, V_GS + h, V_DS) - model_id(p, V_GS - h, V_DS)) / (2 * h)


def residuals(p):
    r6 = [model_id(p, v, VDS_FIG6B) - i for v, i in id6]
    rg = [model_gm(p, v) - g for v, g in gm6]
    return np.r_[W_ID * np.array(r6), W_GM * np.array(rg)]


def rms(p):
    r = residuals(p)
    n6 = len(id6)
    return np.sqrt(np.mean((r[:n6] / W_ID) ** 2)), np.sqrt(np.mean((r[n6:] / W_GM) ** 2))


def rms_out(p, col):
    """rms (mA) against Fig. 7(b), with V_GS labels from column col (0: assumed -1.6:0.2:1.0; 1: from Fig. 6(b))."""
    return np.sqrt(np.mean([(model_id(p, row[col], row[2]) - row[3]) ** 2 for row in out7]))


def paper_rms():
    """The paper's own model line on the same points (interpolated), for comparison."""
    li = load("fig6b_paper_model_ID")
    lg = load("fig6b_paper_model_gm")
    ri = np.interp(id6[:, 0], li[:, 0], li[:, 1]) - id6[:, 1]
    rg = np.interp(gm6[:, 0], lg[:, 0], lg[:, 1] * W_MM) - gm6[:, 1]
    return np.sqrt(np.mean(ri ** 2)), np.sqrt(np.mean(rg ** 2))


def run_fit(p0):
    return least_squares(residuals, p0, bounds=([f[2] for f in FREE], [f[3] for f in FREE]),
                         x_scale=[f[4] for f in FREE], max_nfev=300)


if __name__ == "__main__":
    print(f"{d2.NAME}: W = {d2.W * 1e6:.0f} um, L = {d2.L * 1e9:.0f} nm, C_g''/q = {d2.U:.3e} cm^-2/V (printed form), "
          f"V_OFF = {d2.V_OFF:+.2f} V, mu_LF = {d2.MU_LF * 1e4:.0f} cm^2/Vs, n_acc = {N_ACC:.2e} cm^-2")
    print(f"Fitting to {len(id6)} I_D points and {len(gm6)} G_M points of Fig. 6(b) (V_DS = 5 V)")
    p0 = [f[1] for f in FREE]
    print("start:  rms I_D = {:.2f} mA, G_M = {:.2f} mS".format(*rms(p0)))
    t0 = time.time()
    sol = run_fit(p0)
    print("fitted: rms I_D = {:.2f} mA, G_M = {:.2f} mS".format(*rms(sol.x)), f"  ({time.time() - t0:.0f} s)")
    print("paper's model line on the same points: rms I_D = {:.2f} mA, G_M = {:.2f} mS".format(*paper_rms()))
    sol2 = run_fit([0.05, 3e5, 0.1, 6.0, 4.0])
    print()
    print(f"{'parameter':36s}{'start':>11s}{'fitted':>12s}{'other start':>14s}  unit")
    for (name, start, lo, hi, scale, unit, meaning), v, v2 in zip(FREE, sol.x, sol2.x):
        note = "  <- at its bound" if min(v - lo, hi - v) < 1e-3 * scale else ""
        print(f"{meaning:36s}{start:11.4g}{v:12.4g}{v2:14.4g}  {unit}{note}")
    Rs0, Rd0 = r_access(0.0, L_SG, u0=sol.x[0]), r_access(0.0, L_GD, u0=sol.x[0])
    print(f"R_S0 = {Rs0:.2f} ohm ({Rs0 * W_MM:.2f} ohm.mm), R_D0 = {Rd0:.2f} ohm ({Rd0 * W_MM:.2f} ohm.mm)")
    print()
    print(f"Fig. 7(b) check ({len(out7)} points, not fitted): rms {rms_out(sol.x, 0):.2f} mA with V_GS = -1.6:0.2:1.0 V, "
          f"{rms_out(sol.x, 1):.2f} mA with V_GS from Fig. 6(b)")
    print()
    print(" V_GS   I_D model  measured (mA)    G_M model  measured (mS)")
    for v in [-1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.0]:
        print(f"{v:5.1f} {model_id(sol.x, v, 5.0):9.2f} {np.interp(v, id6[:, 0], id6[:, 1]):9.2f}"
              f"      {model_gm(sol.x, v):9.2f} {np.interp(v, gm6[:, 0], gm6[:, 1]):9.2f}")
    result = {"source": "fit2.py, Step 8b: least squares to the measured points of Nazir et al. 2025 Fig. 6(b)",
              "fitted [CALIBRATE]": dict(zip(NAMES, map(float, sol.x))),
              "fixed": {**FIXED, "v_off": d2.V_OFF, "mu_lf": d2.MU_LF, "th1": cur.THETA_SAT1, "n_max": cur.N_MAX,
                        "mu_b": cur.MU_B, "n_acc_cm2": N_ACC, "l_sg": L_SG, "l_gd": L_GD, "cg_form": "eps/(2d+t) [P]"},
              "rms": {"fig6b_ID_mA": float(rms(sol.x)[0]), "fig6b_gM_mS": float(rms(sol.x)[1]),
                      "fig7b_ID_mA_assumed_labels": float(rms_out(sol.x, 0)),
                      "fig7b_ID_mA_labels_from_fig6b": float(rms_out(sol.x, 1))}}
    with open("params_device2.json", "w") as f:
        json.dump(result, f, indent=2)
    print("\nsaved params_device2.json")
