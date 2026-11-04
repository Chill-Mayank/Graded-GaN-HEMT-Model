"""
device1.py - Device 1 (long-channel PolFET, paper Fig. 5(a), ref. [11]): numbers for the compact model,
and its turn-off voltage V_OFF from the potential balance.
Sources: [P] Nazir et al., IEEE TED 72(9), 2025, Eqs. (3)-(8); [T] Nazir, PhD thesis, IIT Kanpur, 2025, Eqs. 3.9-3.15;
[TCAD] our calibrated Device 1 deck (n13-n16).
"""
import numpy as np

# ---- constants (SI)
Q = 1.602176634e-19           # C
EPS0_SI = 8.8541878128e-12    # F/m

# ---- structure [P, Fig. 5(a); R11]
T_GRADE = 20e-9               # m, graded AlGaN thickness (20 nm)
X_MAX = 0.4                   # Al fraction at the top of the grade (0 at the GaN)
D_ALGAN = 0.0                 # m, no AlGaN barrier: the gate sits on the grade
EPS_R = 9.5                   # [ASSUME, same as Fig. 4(b) fit]

# ---- band-edge numbers [TCAD]
WF_GATE = 6.415               # eV, gate work function [CALIBRATE, n13-n15]
CHI_TOP = 3.281               # eV, electron affinity at the grade top (x = 0.4) [TCAD]
CHI_GAN = 3.959               # eV, electron affinity of GaN [TCAD]
PHI_B = WF_GATE - CHI_TOP     # eV, Schottky barrier height on the grade top
DEC = CHI_GAN - CHI_TOP       # eV, conduction-band step from the grade top to the GaN


def v_total(x_max=X_MAX, t=T_GRADE, eps_r=EPS_R):
    """Potential across the graded layer from its polarization charge [P, Eqs. (3)-(6)].
    sigma_p(y) = a*y + b with a = -0.0984 (x_max/t)^2, b = -0.0593 (x_max/t)  (SI: C/m^4, C/m^3)
    V_Total = -(1/eps) (a t^3/6 + b t^2/2)."""
    a = -0.0984 * (x_max / t) ** 2
    b = -0.0593 * (x_max / t)
    return -(a * t**3 / 6 + b * t**2 / 2) / (eps_r * EPS0_SI)


def v_off(phi_b=PHI_B, dec=DEC, vt=None):
    """Turn-off voltage with physically consistent signs [T, Eq. 3.9 convention; DEVIATION from the printed Eq. (8)]:
    V_OFF = phi_B/q - dE_C/q - V_Total   (no barrier, so the sigma_AlGaN * d_AlGaN term is zero)."""
    vt = v_total() if vt is None else vt
    return phi_b - dec - vt


CG = EPS_R * EPS0_SI / (2 * D_ALGAN + T_GRADE) * 1e-4   # F/cm^2, C_g'' = eps / (2 d_AlGaN + t)  [P]
U = CG / Q                                               # cm^-2 per V


if __name__ == "__main__":
    vt = v_total()
    print(f"phi_B   = {PHI_B:.3f} eV      dE_C = {DEC:.3f} eV")
    print(f"V_Total = {vt:.3f} V        (TCAD polarization, 7.8% stronger: {vt * 3.41 / 3.16:.3f} V)")
    print(f"C_g''   = {CG:.3e} F/cm^2 ->  C_g''/q = {U:.3e} cm^-2 per V")
    print()
    print(f"V_OFF (paper Eq. (3) polarization) = {v_off():+.3f} V")
    print(f"V_OFF (TCAD polarization)          = {v_off(vt=vt * 3.41 / 3.16):+.3f} V")
    print(f"V_OFF, printed signs of Eq. (8)    = {PHI_B + DEC + vt:+.3f} V   <- impossible for a normally-on device")