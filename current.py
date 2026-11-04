"""
current.py - compact model step 4: drain current I_DS of Device 1, intrinsic (no access resistance yet).
Source [P]: Nazir et al., IEEE TED 72(9), 2025, Eqs. (16)-(17) and the mobility and saturation equations below them.
"""
from functools import lru_cache
import numpy as np
from scipy.interpolate import CubicSpline
import device1 as d1
from gate import Q, PHIT, n_total, n_exact

# ---- Device 1 geometry [P, Sec. III and Fig. 5(a); the text's "LG = 700 um" is read as 0.7 um, as in Fig. 5(a)]
W = 150e-6                      # m, total gate width
L = 0.7e-6                      # m, gate length
CG = d1.CG * 1e4                # F/m^2, C_g'' (device1.py gives F/cm^2)
T_EFF = d1.D_ALGAN + d1.T_GRADE # m, d_AlGaN + t, the distance in E_y,eff [P, below Eq. (17)]

# ---- Device 1 V_OFF, both candidates from device1.py; Fig. 6(a) supports -1.26 V (Step 4a)
V_OFF_PAPER = d1.v_off()                                 # -0.988 V, paper polarization
V_OFF_TCADP = d1.v_off(vt=d1.v_total() * 3.41 / 3.16)    # -1.260 V, polarization scaled to the TCAD value

# ---- starting parameters (all fitted in Step 5)
MU_LF = 910e-4                  # m^2/Vs, low-field mobility at 300 K [R11, Hall]; no self-heating, as in the paper
MU_A, MU_B = 0.0, 0.0           # m/V, m^2/V^2, vertical-field mobility terms [ASSUME: start at 0]
THETA_SAT = 0.0                 # 1/V, lateral-field mobility term in mu_eff [ASSUME: start at 0]
VSAT0 = 1e5                     # m/s       [P, Fig. 7(c) legend, Device 2 values used as a start]
THETA_SAT1 = 0.3                #           [P, Fig. 7(c) legend]
N_MAX = 5e16                    # m^-2      [P, Fig. 7(c) legend]
THETA_SAT2 = 0.3                #           [P, Fig. 7(c) legend gives 0.1 ... 0.5]
DELTA_DS = 2.0                  # smoothing exponent in V_DS,eff [ASSUME: the paper gives no value]
N_FLOOR = 1.0                   # cm^-2, n_total is clamped here: Eq. (15) turns negative below -2.71 V [DEVIATION]

# ---- which channel charge feeds Eq. (17)
CHARGE = "exact"                # "exact": root-finding n [DEVIATION, decided 2026-11-04: smooth g_M for the DC fit]
                                # "compact": closed-form n_total, Eqs. (9), (11), (14), (15) [P]
V_TABLE = np.arange(-4.5, 4.0 + 1e-9, 0.01)   # V_GS grid for the exact-n table, 10 mV steps


@lru_cache(maxsize=8)
def n_exact_spline(v_off):
    """Solve the exact n once on the V_TABLE grid, then join the points with a smooth cubic spline of log(n).
    lru_cache keeps the result, so the slow root finding runs only once for each V_OFF."""
    logn = np.log([n_exact(v, v_off, d1.U)[0] for v in V_TABLE])
    return CubicSpline(V_TABLE, logn)


def n_channel(V_GS, v_off=V_OFF_TCADP, charge=CHARGE):
    """Channel sheet density in cm^-2 that feeds Eq. (17)."""
    if charge == "exact":
        return np.exp(n_exact_spline(float(v_off))(V_GS))
    return np.maximum(n_total(V_GS, v_off=v_off, u=d1.U), N_FLOOR)


def ids(V_GS, V_DS, v_off=V_OFF_TCADP, printed=False, charge=CHARGE, mu_lf=MU_LF, mu_a=MU_A, mu_b=MU_B,
        theta_sat=THETA_SAT, vsat0=VSAT0, th1=THETA_SAT1, th2=THETA_SAT2, n_max=N_MAX, delta=DELTA_DS):
    """Intrinsic drain current in A [P, Eq. (17)]:
    I_DS = (W/L) mu_eff C_g'' (V_g0 - V_m + phi_t) V_DS,eff
    printed=True  : V_DSAT uses V_g0, exactly as printed (breaks below V_OFF, where V_g0 < 0)
    printed=False : V_DSAT uses q n / C_g'' instead of V_g0 [DEVIATION, decided 2026-11-04]"""
    Vg0 = V_GS - v_off
    n = n_channel(V_GS, v_off, charge)                            # cm^-2
    Vq = n / d1.U                                                 # V, q n / C_g'': channel charge as a voltage
    EF = Vg0 - Vq                                                 # V, E_F = V_g0 - q n / C_g''  [P]

    vsat = vsat0 / (th1 + (n * 1e4 / n_max) ** th2)               # m/s, V_SAT(n)  [P]

    # V_DSAT needs a mobility; to stay closed-form we use mu (no theta_sat factor) with E_y at the source end,
    # where V_m = E_F, instead of mu_eff [DEVIATION, decided 2026-11-04; no effect while mu_A = mu_B = 0]
    Ey_s = (Vg0 - EF) / T_EFF
    mu_s = mu_lf / (1 + mu_a * Ey_s + mu_b * Ey_s**2)
    Vd = Vg0 if printed else Vq
    with np.errstate(all="ignore"):                               # the printed form divides by zero near -2.6 V
        Vdsat = Vd / (1 + Vd * mu_s / (2 * vsat * L))             # V  [P]
        Vdeff = V_DS / (1 + (V_DS / Vdsat) ** delta) ** (1 / delta)   # V, smooth min(V_DS, V_DSAT)  [P]

    Vm = (2 * EF + Vdeff) / 2                                     # V, average channel potential  [P]
    Ey = (Vg0 - Vm) / T_EFF                                       # V/m, vertical field  [P]
    mu = mu_lf / (1 + mu_a * Ey + mu_b * Ey**2)                   # m^2/Vs  [P]
    mu_eff = mu / np.sqrt(1 + (theta_sat * Vdeff) ** 2)           # m^2/Vs  [P]
    return W / L * mu_eff * CG * (Vg0 - Vm + PHIT) * Vdeff        # A  [P, Eq. (17)]


def i_diffusion(V_GS, v_off=V_OFF_TCADP, charge=CHARGE, mu_lf=MU_LF):
    """Textbook subthreshold (diffusion) current at large V_DS, (W/L) mu q n phi_t, in A: a check, not part of the model."""
    return W / L * mu_lf * Q * n_channel(V_GS, v_off, charge) * 1e4 * PHIT


if __name__ == "__main__":
    print(f"W = {W * 1e6:.0f} um, L = {L * 1e6:.1f} um, C_g'' = {CG:.3e} F/m^2, "
          f"V_OFF = {V_OFF_PAPER:+.3f} V (paper pol.) or {V_OFF_TCADP:+.3f} V (TCAD pol.)")

    print()
    print("Spline check: exact n from the table against a direct root-finding solve")
    for V in [-2.995, -1.2345, -0.5055, 0.0005, 2.7777]:
        direct = n_exact(V, V_OFF_TCADP, d1.U)[0]
        print(f"  V_GS = {V:+8.4f} V   n = {direct:.5e} cm^-2   table error = {n_channel(V) / direct - 1:+.1e}")

    paper = np.loadtxt("paper_data/fig6a_measured_ID.csv", delimiter=",", skiprows=3)
    print()
    print("V_DS = 10 V, I_D in mA     ------------- V_OFF = -1.26 V -------------     V_OFF = -0.99 V    Fig. 6(a)")
    print("  V_GS    printed V_DSAT   fixed, exact n  fixed, compact n  exact/diffusion   fixed, exact n   measured")
    for V, p in paper[::2]:
        a = ids(V, 10.0, printed=True) * 1e3
        b = ids(V, 10.0) * 1e3
        b2 = ids(V, 10.0, charge="compact") * 1e3
        r = ids(V, 10.0) / i_diffusion(V)
        c = ids(V, 10.0, v_off=V_OFF_PAPER) * 1e3
        print(f"{V:+6.2f}   {a:13.4g}   {b:14.4g}   {b2:15.4g}   {r:14.3f}   {c:14.4g}   {p:8.2f}")

    print()
    print("Output curves, V_OFF = -1.26 V, fixed, exact n: I_D (mA)")
    print("  V_GS   V_DS = 0.5 V    1 V      2 V      5 V     10 V")
    for V in [-1.0, 0.0, 1.0, 2.0, 3.0]:
        row = "  ".join(f"{ids(V, vd) * 1e3:7.2f}" for vd in [0.5, 1.0, 2.0, 5.0, 10.0])
        print(f"{V:+6.1f}    {row}")