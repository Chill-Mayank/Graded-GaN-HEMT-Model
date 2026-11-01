"""
charge.py - compact model step 2: sheet carrier density n (cm^-2) against eta_F.
Source [P]: Nazir et al., IEEE TED 72(9), 2025, Eqs. (9), (10), (13) and Fig. 4(a).
"""
import numpy as np
from fermi import fermi_half_numerical

# ---- physical constants (SI) [CODATA]
HBAR = 1.054571817e-34      # J s
M0 = 9.1093837e-31          # kg
KB = 1.380649e-23           # J/K

# ---- device values
T = 300.0                   # K
M_EFF = 0.217               # m*/m0 [P, fitted to the Fig. 4(a) numerical curve; the paper does not state m*]

# effective densities of states [P, text below Eqs. (1) and (9)]
x = M_EFF * M0 * KB * T / (np.pi * HBAR**2)              # m^-2
NC = 0.25 * (2.0 * x) ** 1.5 * 1e-6                      # cm^-3, N_C = 0.25 (2 m* kT / pi hbar^2)^(3/2)
NC2 = 0.8 * x * 1e-4                                     # cm^-2, N_C'' = 0.8 m* kT / (pi hbar^2)


def sheet_exact(eta_F):
    """Exact sheet density: n = (N_C F_1/2(eta_F))^(2/3)  [P, the n^(2/3) transformation]."""
    return (NC * fermi_half_numerical(eta_F)) ** (2.0 / 3.0)


def n_sub(eta_F, gam=2.0 / 3.0):
    """OFF region, eta_F << 0: n_sub = N_C'' exp(gamma eta_F)  [P, Eq. (9)].
    gamma = 2/3 is what the Fig. 4(a) OFF curve uses [P, fitted to the figure]."""
    return NC2 * np.exp(gam * eta_F)


if __name__ == "__main__":
    print(f"m* = {M_EFF} m0  ->  N_C = {NC:.3e} cm^-3,  N_C'' = {NC2:.3e} cm^-2,  N_C^(2/3) = {NC**(2/3):.3e} cm^-2")
    print()
    print(" eta_F     exact       Eq. (9)     error")
    for eta_F in [-10.0, -6.0, -4.0, -2.0, 0.0, 2.0, 5.0, 10.0]:
        exact = sheet_exact(eta_F)
        off = n_sub(eta_F)
        print(f"{eta_F:+6.1f}   {exact:.4e}  {off:.4e}  {100 * (off / exact - 1):+8.2f} %")