"""
gate.py - compact model step 3: sheet density n (cm^-2) against gate voltage V_GS.
Source [P]: Nazir et al., IEEE TED 72(9), 2025, Eqs. (8), (11), (12), (14), (15) and Fig. 4(b).
"""
import numpy as np
from scipy.optimize import brentq
from charge import KB, T, NC2, A_AY, B_AY, THETA1, THETA2, sheet_exact

# ---- constants
Q = 1.602176634e-19         # C
EPS0 = 8.8541878128e-14     # F/cm
PHIT = KB * T / Q           # thermal voltage kT/q, V

# ---- Fig. 4(b) device [P, Fig. 4(b) legend]
V_OFF = -3.0                # V
D_ALGAN = 9e-7              # cm, AlGaN barrier thickness (9 nm)
T_GRADE = 6e-7              # cm, graded layer thickness (6 nm)
EPS_R = 9.5                 # [P, fitted to the ON part of Fig. 4(b); the paper does not state epsilon]

CG = EPS_R * EPS0 / (2 * D_ALGAN + T_GRADE)      # F/cm^2, C_g'' = eps / (2 d_AlGaN + t)  [P, below Eq. (8)]
U = CG / Q                                        # cm^-2 per V: sheet charge per volt of gate drive


def n_exact(V_GS, v_off=V_OFF, u=U):
    """Exact sheet density at gate voltage V_GS: solve
    (N_C F_1/2(eta_F))^(2/3) = (C_g''/q) (V_g0 - phi_t eta_F),  V_g0 = V_GS - V_OFF   [P, Eq. (8) + n^(2/3) transformation]
    for eta_F by root finding. Returns (n, eta_F)."""
    Vg0 = V_GS - v_off
    g = lambda eta: sheet_exact(eta) - u * (Vg0 - PHIT * eta)    # zero at the solution
    eta = brentq(g, -200.0, Vg0 / PHIT + 1e-9)                   # g < 0 at -200; g > 0 just past where the gate line hits zero
    return sheet_exact(eta), eta


# ---- closed-form region solutions (no root finding: this is what makes it a compact model)
GAMMA = 0.55                # [P, fitted to the Fig. 4(b) subthreshold slope; Eq. (9) calls gamma a fitting parameter]
K_ON = (2.0 / (4.24 * np.sqrt(np.pi))) ** (2.0 / 3.0)   # corrected Eq. (10) prefactor [DEVIATION, Step 2b]
DELTA_S, ZETA = 0.1, 0.1    # V, V^2 [P, "Delta_s and zeta = 0.1", read as both 0.1]


def n_sub_V(V_GS, v_off=V_OFF, gamma=GAMMA):
    """OFF region against gate voltage [P, Eq. (9), final form]: n_sub = N_C'' exp(gamma V_g0 / phi_t)."""
    return NC2 * np.exp(gamma * (V_GS - v_off) / PHIT)


def n_on_V(V_GS, corrected=True, v_off=V_OFF, u=U):
    """ON region against gate voltage [P, Eqs. (11)-(12)].
    Eq. (11) solves  (C_g''/q)(V_g0 - phi_t eta) = N [b + eta + sqrt((eta - b)^2 + a)]  for eta.
    As printed, N = N_C'' (Eq. (10)'s prefactor left out); corrected, N = K_ON N_C'' [DEVIATION]."""
    N = (K_ON if corrected else 1.0) * NC2
    Vg0 = V_GS - v_off
    al = u * Vg0 / N - B_AY                       # alpha [P, below Eq. (11)]
    be = u * PHIT / N + 1.0                       # beta
    disc = 4 * (al * be - B_AY) ** 2 - 4 * (be**2 - 1) * (al**2 - A_AY - B_AY**2)
    eta = (2 * (al * be - B_AY) - np.sqrt(disc)) / (2 * (be**2 - 1))    # Eq. (11)
    return u * (Vg0 - PHIT * eta)                 # Eq. (12)


def n_mid_V(V_GS, corrected=True, v_off=V_OFF, u=U):
    """Middle region against gate voltage [P, Eq. (14), rederived].
    Eq. (13) = gate line, with exp(theta2 eta) expanded to second order, gives P1 eta^2 + P2 eta + P3 = 0.
    Printed typos fixed: root uses P1*P3, denominator 2*P1; a, b in P1, P2 mean theta1, theta2.
    As printed the 2/sqrt(pi) of Eq. (13) is left out (N = N_C''); corrected, N = (2/sqrt(pi)) N_C'' [DEVIATION]."""
    N = (2.0 / np.sqrt(np.pi) if corrected else 1.0) * NC2
    Vg0 = V_GS - v_off
    Vm = 0.5 * (Vg0 + DELTA_S + np.sqrt((Vg0 - DELTA_S) ** 2 + ZETA))   # V_g0,mid: smooth max(V_g0, Delta_s)
    w = (u * Vm / N) ** 1.5
    P1 = np.sqrt(np.pi) * THETA2**2 / 200 + THETA1
    P2 = 1.5 * w * PHIT / Vm + THETA2 * np.sqrt(np.pi) / 100
    P3 = -(w - np.sqrt(np.pi) / 100 - 1)
    disc = np.maximum(P2**2 - 4 * P1 * P3, 0.0)  # [DEVIATION] clamp: with Delta_s = 0.1 V it still goes negative below V_GS ~ -2 V
    eta = (-P2 + np.sqrt(disc)) / (2 * P1)
    return u * (Vm - PHIT * eta)


DELTA = 1e6                 # cm^-2, smoothing width of Eq. (15) [ASSUME: not given in the paper]


def n_total(V_GS, corrected=True, v_off=V_OFF, u=U, gamma=GAMMA):
    """All regions joined [P, Eq. (15)]:
    n_sm = n_sub n_mid / (n_sub + n_mid)          (the smaller of OFF and Middle, smoothly)
    n_total = 0.5 [(n_sm + n_on) - sqrt((n_sm - n_on)^2 + Delta^2)]   (the smaller of that and ON, smoothly)
    The printed Eq. (15) has (n_sm - n_on) in the first bracket, which makes n_total negative [DEVIATION: sign fixed]."""
    ns = n_sub_V(V_GS, v_off, gamma)
    nm = n_mid_V(V_GS, corrected, v_off, u)
    no = n_on_V(V_GS, corrected, v_off, u)
    nsm = ns * nm / (ns + nm)
    return 0.5 * ((nsm + no) - np.sqrt((nsm - no) ** 2 + DELTA**2))


def n_total_printed_sign(V_GS):
    """Eq. (15) exactly as printed, for comparison only: 0.5 [(n_sm - n_on) - sqrt((n_sm - n_on)^2 + Delta^2)]."""
    ns, nm, no = n_sub_V(V_GS), n_mid_V(V_GS), n_on_V(V_GS)
    nsm = ns * nm / (ns + nm)
    return 0.5 * ((nsm - no) - np.sqrt((nsm - no) ** 2 + DELTA**2))


if __name__ == "__main__":
    print(f"C_g'' = {CG:.3e} F/cm^2  ->  C_g''/q = {U:.3e} cm^-2 per V,  phi_t = {PHIT * 1000:.2f} mV")
    paper = np.loadtxt("paper_data/fig4b_markers.csv", delimiter=",", skiprows=3)
    print()
    print("  V_GS    eta_F      exact n       paper n      ratio")
    for V in [-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0]:
        n, eta = n_exact(V)
        p = paper[np.argmin(abs(paper[:, 0] - V)), 1]
        print(f"{V:+6.2f}  {eta:+7.2f}   {n:.4e}   {p:.4e}   {n / p:6.3f}")

    print()
    print("  V_GS     exact      OFF(9)    Mid(14)pr  Mid(14)co  ON(11)pr   ON(11)co")
    for V in [-4.0, -3.5, -3.0, -2.5, -2.0, -1.0, 0.0]:
        n, _ = n_exact(V)
        print(f"{V:+6.2f}  {n:.3e}  {n_sub_V(V):.3e}  {n_mid_V(V, False):.3e}  {n_mid_V(V):.3e}  "
              f"{n_on_V(V, False):.3e}  {n_on_V(V):.3e}")

    print()
    print("  V_GS     paper      Eq.(15) fixed  ratio    Eq.(15) as printed")
    for V in [-4.3, -4.0, -3.5, -3.2, -3.0, -2.8, -2.5, -2.0, -1.0, 0.0]:
        p = paper[np.argmin(abs(paper[:, 0] - V)), 1]
        nt = n_total(V)
        print(f"{V:+6.2f}  {p:.3e}   {nt:.3e}   {nt / p:6.3f}    {n_total_printed_sign(V):+.3e}")