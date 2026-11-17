"""
current.py - compact model step 4: drain current I_DS of Device 1, intrinsic (no access resistance yet).
Step 7a: ids() split into bias_state() + ids(), so the terminal charges reuse E_F and V_DS,eff; I_DS is unchanged.
Source [P]: Nazir et al., IEEE TED 72(9), 2025, Eqs. (16)-(17) and the mobility and saturation equations below them.
Channel charge: the simultaneous solution of the paper's two charge equations [P, Sec. II-B: n = N_C F_1/2(eta_F)
with the n^(2/3) transformation, and the gate line of Eq. (8)], i.e. the paper's own "numerical solution" [P, Fig. 4(a)].
Numerics: Brent's root finder [R. P. Brent, Algorithms for Minimization without Derivatives, 1973] (scipy brentq,
in gate.py), tabulated once and joined by a quintic interpolating spline [C. de Boor, A Practical Guide to Splines, 1978]
(scipy make_interp_spline, k = 5), so n, G_M, G_M2 and G_M3 are continuous.
"""
from functools import lru_cache
import numpy as np
from scipy.interpolate import make_interp_spline
from device import dev as d1           # Step 8a: device1 or device2 (GAN_DEVICE)
from gate import Q, PHIT, n_total, n_exact

# ---- Device 1 geometry [P, Sec. III and Fig. 5(a); the text's "LG = 700 um" is read as 0.7 um, as in Fig. 5(a)]
W = d1.W                        # m, total gate width (Device 1: 150 um)
L = d1.L                        # m, gate length (Device 1: 0.7 um)
CG = d1.CG * 1e4                # F/m^2, C_g'' (device1.py gives F/cm^2)
T_EFF = d1.D_ALGAN + d1.T_GRADE # m, d_AlGaN + t, the distance in E_y,eff [P, below Eq. (17)]

# ---- Device 1 V_OFF, both candidates from device1.py; Fig. 6(a) supports -1.26 V (Step 4a)
V_OFF_PAPER = d1.v_off() if hasattr(d1, "v_off") else d1.V_OFF   # Device 1: -0.988 V, paper polarization
V_OFF_TCADP = d1.V_OFF          # Device 1: -1.260 V (TCAD polarization); Device 2: -1.01 V [TCAD n8 + DIBL]

# ---- starting parameters (all fitted in Step 5)
MU_LF = d1.MU_LF                # m^2/Vs, low-field mobility (Device 1: 910, Hall [R11]); no self-heating, as in the paper
MU_A, MU_B = 0.0, 0.0           # m/V, m^2/V^2, vertical-field mobility terms [ASSUME: start at 0]
THETA_SAT = 0.0                 # 1/V, lateral-field mobility term in mu_eff [ASSUME: start at 0]
VSAT0 = 1e5                     # m/s       [P, Fig. 7(c) legend, Device 2 values used as a start]
THETA_SAT1 = 0.3                #           [P, Fig. 7(c) legend]
N_MAX = 5e16                    # m^-2      [P, Fig. 7(c) legend]
THETA_SAT2 = 0.3                #           [P, Fig. 7(c) legend gives 0.1 ... 0.5]
DELTA_DS = 2.0                  # smoothing exponent in V_DS,eff [ASSUME: the paper gives no value]
N_FLOOR = 1.0                   # cm^-2, n_total is clamped here: Eq. (15) turns negative below -2.71 V [DEVIATION]

# ---- which channel charge feeds Eq. (17)
CHARGE = "exact"                # "exact": numerical solution of the two charge equations [P, Sec. II-B; Fig. 4(a)]
                                #   [DEVIATION from Eqs. (9)-(15), decided 2026-11-09: printed Eq. (14) has no real root
                                #   for V_g0,mid < ~0.8 V (Device 1), and the paper's own Fig. 4(b) and Fig. 6(a) lines
                                #   follow the numerical solution there, not the printed closed form]
                                # "compact": closed-form n_total, Eqs. (9), (11), (14), (15) [P], kept for Fig. 4
V_TABLE = np.arange(-4.5, 4.0 + 1e-9, 0.01)   # V_GS grid for the exact-n table, 10 mV steps


GAMMA_EX = getattr(d1, "GAMMA_EX", None)   # Step 8c: gamma of Eq. (9) in the numerical charge; None = m 1 (Device 1)


@lru_cache(maxsize=8)
def n_exact_spline(v_off, gamma=None):
    """Solve the exact n once on the V_TABLE grid, then join the points with a quintic spline of log(n)
    [de Boor 1978]. Quintic (not cubic) keeps the third derivative continuous, which G_M3 needs.
    lru_cache keeps the result, so the slow root finding runs only once for each V_OFF."""
    logn = np.log([n_exact(v, v_off, d1.U, gamma)[0] for v in V_TABLE])
    return make_interp_spline(V_TABLE, logn, k=5)


def n_channel(V_GS, v_off=V_OFF_TCADP, charge=CHARGE, gamma=GAMMA_EX):
    """Channel sheet density in cm^-2 that feeds Eq. (17)."""
    if charge == "exact":
        return np.exp(n_exact_spline(float(v_off), None if gamma is None else float(gamma))(V_GS))
    return np.maximum(n_total(V_GS, v_off=v_off, u=d1.U), N_FLOOR)


def bias_state(V_GS, V_DS, v_off=V_OFF_TCADP, printed=False, charge=CHARGE, mu_lf=MU_LF, mu_a=MU_A, mu_b=MU_B,
               theta_sat=THETA_SAT, vsat0=VSAT0, th1=THETA_SAT1, th2=THETA_SAT2, n_max=N_MAX, delta=DELTA_DS,
               gamma=GAMMA_EX):
    """Every intermediate quantity of Eq. (17) at one intrinsic bias, as a dict:
    Vg0, Vq (= q n / C_g''), EF, Vdsat, Vdeff, Vm in V; n in cm^-2; mu_eff in m^2/Vs; I in A.
    printed=True  : V_DSAT uses V_g0, exactly as printed (breaks below V_OFF, where V_g0 < 0)
    printed=False : V_DSAT uses q n / C_g'' instead of V_g0 [DEVIATION, decided 2026-11-04]
    The channel potentials of Sec. II-C are V_ch,s = E_F and V_ch,d = E_F + V_DS,eff: then
    V_m = (V_ch,s + V_ch,d)/2 is the same V_m as in Eq. (17) [P, both define V_m] (used by terminal_charges.py)."""
    Vg0 = V_GS - v_off
    n = n_channel(V_GS, v_off, charge, gamma)                     # cm^-2
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
    I = W / L * mu_eff * CG * (Vg0 - Vm + PHIT) * Vdeff           # A  [P, Eq. (17)]
    return {"Vg0": Vg0, "n": n, "Vq": Vq, "EF": EF, "Vdsat": Vdsat, "Vdeff": Vdeff, "Vm": Vm, "mu_eff": mu_eff, "I": I}


def ids(V_GS, V_DS, *args, **kw):
    """Intrinsic drain current in A [P, Eq. (17)]:
    I_DS = (W/L) mu_eff C_g'' (V_g0 - V_m + phi_t) V_DS,eff
    Takes the same arguments, in the same order, as before Step 7a (v_off, printed, charge, mu_lf, ...)."""
    return bias_state(V_GS, V_DS, *args, **kw)["I"]


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

    print()
    print("Smoothness check: jump of d3(log n)/dV3 across table knots (0 for a curve with a continuous 3rd derivative)")
    d3 = n_exact_spline(float(V_OFF_TCADP)).derivative(3)
    for V in [-1.0, -0.5, 0.0, 1.0]:
        print(f"  V_GS = {V:+5.2f} V   relative jump = {abs(d3(V + 1e-7) - d3(V - 1e-7)) / abs(d3(V)):.1e}")

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