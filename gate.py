"""
gate.py - compact model step 3: sheet density n (cm^-2) against gate voltage V_GS.
Source [P]: Nazir et al., IEEE TED 72(9), 2025, Eqs. (8), (11), (12), (14), (15) and Fig. 4(b).
"""
import numpy as np
from scipy.optimize import brentq
from charge import KB, T, sheet_exact

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


def n_exact(V_GS):
    """Exact sheet density at gate voltage V_GS: solve
    (N_C F_1/2(eta_F))^(2/3) = (C_g''/q) (V_g0 - phi_t eta_F),  V_g0 = V_GS - V_OFF   [P, Eq. (8) + n^(2/3) transformation]
    for eta_F by root finding. Returns (n, eta_F)."""
    Vg0 = V_GS - V_OFF
    g = lambda eta: sheet_exact(eta) - U * (Vg0 - PHIT * eta)    # zero at the solution
    eta = brentq(g, -200.0, Vg0 / PHIT)                          # g < 0 at -200, g > 0 where the gate line hits zero
    return sheet_exact(eta), eta


if __name__ == "__main__":
    print(f"C_g'' = {CG:.3e} F/cm^2  ->  C_g''/q = {U:.3e} cm^-2 per V,  phi_t = {PHIT * 1000:.2f} mV")
    paper = np.loadtxt("paper_data/fig4b_markers.csv", delimiter=",", skiprows=3)
    print()
    print("  V_GS    eta_F      exact n       paper n      ratio")
    for V in [-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0]:
        n, eta = n_exact(V)
        p = paper[np.argmin(abs(paper[:, 0] - V)), 1]
        print(f"{V:+6.2f}  {eta:+7.2f}   {n:.4e}   {p:.4e}   {n / p:6.3f}")