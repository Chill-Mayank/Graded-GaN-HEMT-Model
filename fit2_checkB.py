"""fit2_checkB.py - Step 8b check B: as fit2.py, with V_OFF also free [CALIBRATE], bounds -1.6 .. -0.6 V.
Run: GAN_DEVICE=2 python3 fit2_checkB.py"""
import numpy as np, time, os
from scipy.optimize import least_squares
assert os.environ.get("GAN_DEVICE") == "2"
import fit2 as f2
from access import id_ext

def model_id(p, V_GS, V_DS):
    return id_ext(V_GS, V_DS, v_off=p[5], **f2.FIXED, **dict(zip(f2.NAMES, p[:5]))) * 1e3
def gm(p, v, h=0.01):
    return (model_id(p, v + h, 5.0) - model_id(p, v - h, 5.0)) / (2 * h)
def res(p):
    return np.r_[[model_id(p, v, 5.0) - i for v, i in f2.id6], [2.0 * (gm(p, v) - g) for v, g in f2.gm6]]
import json
pA = list(json.load(open("params_device2.json"))["fitted [CALIBRATE]"].values())
lo = [f[2] for f in f2.FREE] + [-1.6]; hi = [f[3] for f in f2.FREE] + [-0.6]
t0 = time.time()
s = least_squares(res, pA + [-1.01], bounds=(lo, hi), x_scale=[f[4] for f in f2.FREE] + [0.2], max_nfev=200)
r = s.fun; n = len(f2.id6)
print(f"B: V_OFF = {s.x[5]:+.3f} V; rms I_D = {np.sqrt(np.mean(r[:n]**2)):.2f} mA, G_M = {np.sqrt(np.mean((r[n:]/2)**2)):.2f} mS ({time.time()-t0:.0f} s)")
print("   fitted:", dict(zip(f2.NAMES, np.round(s.x[:5], 4))))
print("   Fig. 7(b) rms (labels from Fig. 6(b)):", round(float(np.sqrt(np.mean([(model_id(s.x, row[1], row[2]) - row[3])**2 for row in f2.out7]))), 2), "mA")
for v in [-1.2, -1.0, -0.8]:
    print(f"   {v:+.1f} V: I_D {model_id(s.x, v, 5.0):.2f} (meas {np.interp(v, f2.id6[:,0], f2.id6[:,1]):.2f}) mA, G_M {gm(s.x, v):.2f} (meas {np.interp(v, f2.gm6[:,0], f2.gm6[:,1]):.2f}) mS")
