"""
fit2_dibl.py - Step 11a: Device 2 refit to Fig. 6(b) with DIBL [DEVIATION D5, user choice 2B, 18 Nov].
Run:  GAN_DEVICE=2 python3 fit2_dibl.py      -> params_device2_dibl.json
Model: the Verilog-A twin (va_twin.py: Eq. (2) F_1/2, the model Spectre runs), with
  V_OFF(V_DS,i) = V_OFF,n8 - SIGMA (V_DS,i - 0.05 V),  V_OFF,n8 = -0.63 V [TCAD n8], SIGMA = 0.38/4.95 [TCAD n6 -> n7],
  V_DS,i the intrinsic drain voltage. Nothing about the DIBL is fitted.
Same data, weights, free parameters and bounds as fit2.py (Step 8c); Fig. 7(b) is checked, not fitted.
"""
import json
import os
import time
import numpy as np
from scipy.optimize import least_squares
assert os.environ.get("GAN_DEVICE") == "2", "run as: GAN_DEVICE=2 python3 fit2_dibl.py"
import fit2 as F
from va_twin import Twin, card

SIGMA, VOFF_N8, VREF = 0.38 / 4.95, -0.63, 0.05
BASE = {**card(2), "VOFF": VOFF_N8, "SIGMA": SIGMA, "VDSREF": VREF}
KEY = {"u0": "U0", "vsat0": "VSAT0", "th2": "TH2", "delta": "DELTA", "beta": "BETA", "gamma": "GAMMA"}


def twin(p):
    return Twin({**BASE, **{KEY[n]: v for n, v in zip(F.NAMES, p)}})


def model_id(p, vg, vd, t=None):
    return (t or twin(p)).id_ext(vg, vd) * 1e3


def model_gm(p, vg, t=None, h=0.01):
    t = t or twin(p)
    return (t.id_ext(vg + h, F.VDS_FIG6B) - t.id_ext(vg - h, F.VDS_FIG6B)) / (2 * h) * 1e3


def residuals(p):
    t = twin(p)
    r6 = [model_id(p, v, F.VDS_FIG6B, t) - i for v, i in F.id6]
    rg = [model_gm(p, v, t) - g for v, g in F.gm6]
    return np.r_[F.W_ID * np.array(r6), F.W_GM * np.array(rg)]


def rms(p):
    r = residuals(p); n = len(F.id6)
    return np.sqrt(np.mean((r[:n] / F.W_ID) ** 2)), np.sqrt(np.mean((r[n:] / F.W_GM) ** 2))


def rms_out(p, col):
    t = twin(p)
    return np.sqrt(np.mean([(model_id(p, row[col], row[2], t) - row[3]) ** 2 for row in F.out7]))


if __name__ == "__main__":
    p0 = list(json.load(open("params_device2.json"))["fitted [CALIBRATE]"].values())
    print(f"start (Step 8c values, DIBL on): rms I_D = {rms(p0)[0]:.2f} mA, G_M = {rms(p0)[1]:.2f} mS")
    t0 = time.time()
    sol = least_squares(residuals, p0, bounds=([f[2] for f in F.FREE], [f[3] for f in F.FREE]),
                        x_scale=[f[4] for f in F.FREE], max_nfev=200)
    ri, rg = rms(sol.x)
    print(f"fitted: rms I_D = {ri:.2f} mA, G_M = {rg:.2f} mS  ({time.time() - t0:.0f} s); "
          f"paper's model {F.paper_rms()[0]:.2f} / {F.paper_rms()[1]:.2f}; Step 8c without DIBL 0.18 / 0.86")
    for (name, start, lo, hi, scale, unit, meaning), v, v0 in zip(F.FREE, sol.x, p0):
        print(f"  {meaning:36s} Step 8c {v0:10.4g}   with DIBL {v:10.4g} {unit}")
    o0, o1 = rms_out(sol.x, 0), rms_out(sol.x, 1)
    print(f"Fig. 7(b) check (not fitted): {o0:.2f} mA with printed-step labels, {o1:.2f} mA with labels from Fig. 6(b)"
          f"  (Step 8c without DIBL: 7.68 / 3.78)")
    json.dump({"source": "fit2_dibl.py, Step 11a: Fig. 6(b) refit with DIBL, Verilog-A twin (Eq. (2) F_1/2)",
               "fitted [CALIBRATE]": dict(zip(F.NAMES, map(float, sol.x))),
               "dibl [TCAD]": {"v_off_n8": VOFF_N8, "sigma": SIGMA, "vds_ref": VREF},
               "rms": {"fig6b_ID_mA": float(ri), "fig6b_gM_mS": float(rg), "fig7b_ID_mA_assumed_labels": float(o0),
                       "fig7b_ID_mA_labels_from_fig6b": float(o1)}},
              open("params_device2_dibl.json", "w"), indent=2)
    print("saved params_device2_dibl.json")
