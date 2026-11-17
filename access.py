"""
access.py - compact model step 5a: non-linear access resistance and the external drain current of Device 1.
Source: the paper uses "the nonlinear access resistance model described in [25] and [26]" [P, Sec. III];
the equations are from the first author's thesis [T, Eqs. 5.1-5.3], R = R0 / [1 - (I/I_sat)^beta]^(1/beta).
"""
import numpy as np
from scipy.optimize import brentq
from gate import Q
from current import ids, W, V_TABLE
from device import dev                 # Step 8a

# ---- access regions of Device 1 [P, Fig. 5(a)]
L_SG = dev.L_SG                 # m, source-to-gate access length (Device 1: 0.5 um)
L_GD = dev.L_GD                 # m, gate-to-drain access length (Device 1: 3.5 um)
N_ACC = dev.N_ACC               # cm^-2, access sheet density [TCAD] (Device 1: 1.13e13, run n2)
Q_ACC = Q * N_ACC * 1e4         # C/m^2, Q_acc = q n_acc

# ---- starting values (fitted in Step 5b); the same values for the source and the drain side (decided 2026-11-04)
U0_ACC = dev.U0_ACC             # m^2/Vs, access mobility start (Device 1: 910, R11 Hall)
VSAT_ACC = 1e5                  # m/s, access saturation velocity [ASSUME: start]
BETA_ACC = 2.0                  # exponent [ASSUME: start; not given in [T]]


def r_access(I, l_acc, u0=U0_ACC, vsat_acc=VSAT_ACC, beta=BETA_ACC):
    """Access resistance in ohm at current I in A [T, Eqs. 5.1-5.3]:
    R0 = L / (W Q_acc U0), I_sat = W Q_acc V_SAT,acc, R = R0 / [1 - (I/I_sat)^beta]^(1/beta)."""
    R0 = l_acc / (W * Q_ACC * u0)
    I_sat = W * Q_ACC * vsat_acc
    return R0 / (1 - (I / I_sat) ** beta) ** (1 / beta)


def id_ext(V_GS, V_DS, u0=U0_ACC, vsat_acc=VSAT_ACC, beta=BETA_ACC, **kw):
    """External drain current in A: solve I = ids(V_GS - I R_S, V_DS - I (R_S + R_D)) for I.
    kw goes on to ids() (vsat0, th2, mu_a, theta_sat, ...)."""
    if V_DS <= 0:
        return 0.0
    R0_sum = (L_SG + L_GD) / (W * Q_ACC * u0)
    I_sat = W * Q_ACC * vsat_acc

    def f(I):                                         # zero at the solution; rises with I
        Rs = r_access(I, L_SG, u0, vsat_acc, beta)
        Rd = r_access(I, L_GD, u0, vsat_acc, beta)
        vgs_int = max(V_GS - I * Rs, V_TABLE[0])      # gate-source voltage the channel sees (kept inside the n table)
        vds_int = max(V_DS - I * (Rs + Rd), 0.0)      # drain-source voltage the channel sees (never below 0)
        return I - ids(vgs_int, vds_int, **kw)

    I_hi = min(V_DS / R0_sum, I_sat * (1 - 1e-9))     # here the channel sees V_DS <= 0, so f > 0
    return brentq(f, 0.0, I_hi, xtol=1e-15, rtol=1e-12)


if __name__ == "__main__":
    R0_S = r_access(0.0, L_SG)
    R0_D = r_access(0.0, L_GD)
    I_sat = W * Q_ACC * VSAT_ACC
    print(f"Q_acc = {Q_ACC * 1e3:.2f} mC/m^2, R_S0 = {R0_S:.2f} ohm ({R0_S * W * 1e3:.2f} ohm.mm), "
          f"R_D0 = {R0_D:.2f} ohm ({R0_D * W * 1e3:.2f} ohm.mm), I_sat = {I_sat * 1e3:.0f} mA")

    print()
    print("Access resistance against current (starting values)")
    print("  I (mA)    R_S (ohm)   R_D (ohm)")
    for I in [0, 50, 100, 150, 200, 250]:
        print(f"  {I:5.0f}    {r_access(I * 1e-3, L_SG):8.2f}    {r_access(I * 1e-3, L_GD):8.2f}")

    paper = np.loadtxt("paper_data/fig6a_measured_ID.csv", delimiter=",", skiprows=3)
    print()
    print("V_DS = 10 V, I_D (mA)   no access R   with access R   measured")
    for V, p in paper[6::2]:
        print(f"  V_GS = {V:+5.2f}      {ids(V, 10.0) * 1e3:8.2f}      {id_ext(V, 10.0) * 1e3:8.2f}      {p:8.2f}")

    out = np.loadtxt("paper_data/fig7a_measured.csv", delimiter=",", skiprows=3)
    print()
    print("V_DS = 1 V, I_D (mA)    no access R   with access R   measured")
    for V in [0.0, 1.0, 2.0, 3.0]:
        m = out[(out[:, 0] == V) & (out[:, 1] == 1.0), 2][0]
        print(f"  V_GS = {V:+5.2f}      {ids(V, 1.0) * 1e3:8.2f}      {id_ext(V, 1.0) * 1e3:8.2f}      {m:8.2f}")