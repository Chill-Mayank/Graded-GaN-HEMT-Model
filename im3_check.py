"""
im3_check.py - Step 10d check: two-tone IM3 with the 1 H bias feeds (va/hb2_device1.raw, baseband sees ~50 ohm)
against 100 nH feeds (va/hb2b_device1.raw, baseband ~6 ohm). Run: python3 im3_check.py
"""
import numpy as np
from fig10 import two_tone
a, b = two_tone("va/hb2_device1.raw"), two_tone("va/hb2b_device1.raw")
m = np.loadtxt("paper_data/fig10b_measured_IM3_dBm.csv", delimiter=",", skiprows=3)
print(" P_DEL/tone   IM3 1 H feeds   IM3 100 nH feeds   measured   (dBm)")
for v in (-1.4, 0.6, 2.6, 4.6, 6.5):
    print(f"   {v:+5.1f}      {np.interp(v, a[:, 0], a[:, 2]):8.1f}        {np.interp(v, b[:, 0], b[:, 2]):8.1f}       "
          f"{np.interp(v, m[:, 0], m[:, 1]):7.1f}")
print("P_FUND at +0.6 dBm: %.2f / %.2f dBm" % (np.interp(0.6, a[:, 0], a[:, 1]), np.interp(0.6, b[:, 0], b[:, 1])))
