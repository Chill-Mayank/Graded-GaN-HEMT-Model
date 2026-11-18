"""
va_sp2_compare.py - Step 12c: Spectre S-parameters of graded_hemt.va for Device 2 (va/sp_device2_vg0..3.s2p, V_DS = 3 V,
V_GS = -3, -2, -1, 0 V) against rf_twin.py with the same parameters. Checks the DIBL and V_OFF-shift terms in AC.
Run:  GAN_DEVICE=2 python3 va_sp2_compare.py        Pass: H21, |S21| within 0.05 dB; Y within 0.5 % [ASSUME]
"""
import json
import numpy as np
from rf_twin import RF
from fit2_dibl import BASE, KEY
from va_sp_compare import read_s2p, s_to_y
from va_netlists import SP2_VG

P = json.load(open("params_device2_dibl.json"))["fitted [CALIBRATE]"]
R = json.load(open("rf_device2.json"))
OV = {**{k: BASE[k] for k in ("VOFF", "SIGMA", "VDSREF")}, **{KEY[n]: v for n, v in P.items()},
      "DVOFF": R["DVOFF_V"], "U0": R["U0_m2Vs"]}
r = RF(2, rg=R["RG_ohm"], cgs=R["CGSFR_F"], cgd=R["CGDFR_F"], cds=R["CDSFR_F"], over=OV)
print("Device 2, V_DS = 3 V: Spectre (graded_hemt.va) against rf_twin.py")
for k, vg in enumerate(SP2_VG):
    f, S, z0 = read_s2p(f"va/sp_device2_vg{k}.s2p")
    Ys = s_to_y(S, z0)
    Yt, _ = r.y2(f, vg, 3.0)
    St = np.array([np.linalg.solve(np.eye(2) + z0 * y, np.eye(2) - z0 * y) for y in Yt])
    ds21 = np.max(np.abs(20 * np.log10(np.abs(S[:, 1, 0]) / np.abs(St[:, 1, 0]))))
    hs, ht = RF.gains(Ys)[0], RF.gains(Yt)[0]
    dh = np.max(np.abs(20 * np.log10(hs / ht))) if vg > -3 else float("nan")
    dy = np.max(np.abs(Ys - Yt) / np.abs(Yt).max(axis=(1, 2))[:, None, None])
    print(f"  V_GS {vg:+.0f} V: |S21| {ds21:.4f} dB, H21 {dh:.4f} dB, Y {100 * dy:.3f} %   "
          f"(g_m at 1 GHz: Spectre {Ys[0, 1, 0].real * 1e3:.2f}, twin {Yt[0, 1, 0].real * 1e3:.2f} mS)")
