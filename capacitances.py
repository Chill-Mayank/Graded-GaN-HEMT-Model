"""
capacitances.py - compact model step 7c: intrinsic capacitances of Device 1 from the terminal charges.

Definitions [D. E. Ward and R. W. Dutton, IEEE JSSC 13(5), 703-708, 1978], terminals i, j in {g, d, s}:
  C_ii = dQ_i/dV_i      C_ij = -dQ_i/dV_j  (i != j)
so the capacitances are positive in normal operation, and the matrix is not symmetric (C_dg != C_gd): a
charge-based model has transcapacitances. A circuit simulator forms exactly these from ddt(Q) in Verilog-A.

Charges: terminal_charges.py (Steps 7a-7b), intrinsic device (no access resistance). They depend on
V_GS = V_G - V_S and V_DS = V_D - V_S only, so a change of V_S is a change of both, with the sign reversed.

Derivatives: central differences with step h, improved by one Richardson extrapolation, (4 D(h/2) - D(h))/3
  [W. H. Press et al., Numerical Recipes, 3rd ed., Cambridge, 2007, Sec. 5.7]. The charges are smooth in V_GS
  (quintic spline, C^4) and in V_DS for V_DS > 0.
The model is written for V_DS >= 0 only: V_DS,eff contains (V_DS/V_DSAT)^delta with non-integer delta, which
has no real value for V_DS < 0. Every difference here keeps V_DS > 0. [open item: the Verilog-A port needs
the standard source-drain interchange for V_DS < 0]

Checks run by __main__:
  1. Charge conservation: each column of dQ_i/dV_j sums to zero, because Q_g + Q_d + Q_s = 0 [P, Sec. II-C].
  2. Reference invariance: each row of dQ_i/dV_j sums to zero, because the charges depend only on voltage
     differences. (With the sign convention above, each diagonal C_ii equals the sum of its row's C_ij.)
  3. V_DS = 0: C_gg = q W L dn/dV_GS, since Eq. (19) gives Q_g = q W L n there (Step 7a, check 4); dn/dV_GS is
     taken exactly from the spline of the numerical charge solution [CALC].
  4. AC symmetry at V_DS -> 0 [C. C. McAndrew, IEEE TED 53(9), 2202-2206, 2006]: a source-drain symmetric model
     gives C_gd = C_gs there. In the paper's formulation u_s = q n(V_GS)/C_g'' but u_d = u_s - V_DS,eff, so the
     drain end follows V_D with the full C_g'' while the source end follows V_S with C_g'' a, where
     a = (q/C_g'') dn/dV_GS <= 1 (degenerate electrons, the n^(2/3) form). Then C_gs -> (2a - 1) C_g''WL/2 and
     C_gd -> C_g''WL/2, so C_gd/C_gs -> 1/(2a - 1) [CALC]: equal to 1 only when a = 1. A property of the
     paper's model, not of this code; recorded as an open item.
  5. Long-channel pinch-off limit of Eqs. (19) and the D.6 Q_d (phi_t = 0, u_d = 0), per unit C_g'' W L:
     dQ_g/du_s = 2/3, dQ_g/du_d = 0, -dQ_d/du_s = 4/15, the textbook saturation values
     [Y. Tsividis and C. McAndrew, Operation and Modeling of the MOS Transistor, 3rd ed., Oxford, 2011].
  6. Step size: h = 1 mV and h = 0.5 mV agree after extrapolation.
"""
import numpy as np
import current as cur
import terminal_charges as tc
from gate import Q

NAMES = ("g", "d", "s")


def charges(V_GS, V_DS, **kw):
    """(Q_g, Q_d, Q_s) in C at one intrinsic bias."""
    c = tc.channel_potentials(V_GS, V_DS, **kw)
    return np.array([tc.qg(c["us"], c["ud"]), tc.qd(c["us"], c["ud"]), tc.qs(c["us"], c["ud"])])


def charge_derivatives(V_GS, V_DS, h, **kw):
    """3x3 matrix [i, j] = dQ_i/dV_j (F) by central differences, i in (g, d, s), j in (G, D, S).
    V_G moves V_GS; V_D moves V_DS; V_S moves both, with the opposite sign."""
    cols = []
    for dgs, dds in [(1.0, 0.0), (0.0, 1.0), (-1.0, -1.0)]:
        q_plus = charges(V_GS + dgs * h, V_DS + dds * h, **kw)
        q_minus = charges(V_GS - dgs * h, V_DS - dds * h, **kw)
        cols.append((q_plus - q_minus) / (2 * h))
    return np.array(cols).T


def derivative_matrix(V_GS, V_DS, h=1e-3, **kw):
    """dQ_i/dV_j in F, Richardson-extrapolated central differences, (4 D(h/2) - D(h))/3."""
    return (4 * charge_derivatives(V_GS, V_DS, h / 2, **kw) - charge_derivatives(V_GS, V_DS, h, **kw)) / 3


def cap_matrix(V_GS, V_DS, h=1e-3, **kw):
    """Capacitance matrix in F, C_ii = dQ_i/dV_i and C_ij = -dQ_i/dV_j, rows and columns in the order (g, d, s)."""
    D = derivative_matrix(V_GS, V_DS, h, **kw)
    C = -D
    np.fill_diagonal(C, np.diag(D))
    return C


def cgg_at_vds0(V_GS, h=1e-3, **kw):
    """C_gg = dQ_g/dV_GS at V_DS = 0 (only V_GS is moved, so V_DS stays at 0)."""
    def qg0(v):
        return charges(v, 0.0, **kw)[0]
    d1 = (qg0(V_GS + h) - qg0(V_GS - h)) / (2 * h)
    d2 = (qg0(V_GS + h / 2) - qg0(V_GS - h / 2)) / h
    return (4 * d2 - d1) / 3


if __name__ == "__main__":
    kw = tc.model_params()
    fF = 1e15
    print(f"Device 1: C_g'' W L = {tc.QG_SCALE * fF:.2f} fF, V_OFF = {kw['v_off']:+.3f} V; Step 5e parameters, intrinsic")

    print()
    print("Check 1-2: sums of dQ_i/dV_j over each column (charge conservation) and each row (reference invariance),")
    print("           relative to the largest entry")
    for V, VD in [(-1.0, 0.5), (0.0, 0.05), (0.0, 10.0), (3.0, 10.0)]:
        D = derivative_matrix(V, VD, **kw)
        scale = np.abs(D).max()
        print(f"  V_GS = {V:+5.2f} V, V_DS = {VD:5.2f} V:   max |column sum| = {np.abs(D.sum(axis=0)).max() / scale:.1e}"
              f"   max |row sum| = {np.abs(D.sum(axis=1)).max() / scale:.1e}")

    print()
    print("Check 3: V_DS = 0, C_gg against q W L dn/dV_GS (dn/dV_GS exact from the spline)")
    spl = cur.n_exact_spline(float(kw["v_off"]))
    dlogn = spl.derivative(1)
    for V in [-1.5, -1.0, 0.0, 3.0]:
        n = np.exp(spl(V))
        ref = Q * cur.W * cur.L * n * dlogn(V) * 1e4
        cgg = cgg_at_vds0(V, **kw)
        print(f"  V_GS = {V:+5.2f} V   C_gg = {cgg * fF:10.4f} fF   q W L dn/dV_GS = {ref * fF:10.4f} fF   "
              f"relative difference = {abs(cgg / ref - 1):.1e}")

    print()
    print("Check 4: AC symmetry, C_gd/C_gs as V_DS -> 0, against 1/(2a - 1) with a = (q/C_g'') dn/dV_GS")
    for V in [0.0, -1.0]:
        a = cgg_at_vds0(V, **kw) / tc.QG_SCALE
        print(f"  V_GS = {V:+5.2f} V:  a = {a:.4f}, predicted limit 1/(2a - 1) = {1 / (2 * a - 1):.4f}")
        for VD in [0.01, 0.002, 0.0004]:
            C = cap_matrix(V, VD, h=VD / 4, **kw)
            print(f"     V_DS = {VD * 1e3:4.1f} mV   C_gs = {C[0, 2] * fF:8.3f} fF   C_gd = {C[0, 1] * fF:8.3f} fF   "
                  f"C_gd/C_gs = {C[0, 1] / C[0, 2]:.4f}")

    print()
    print("Check 5: pinch-off limit of the charge formulas, phi_t = 0, u_s = 1 V, u_d = 0, per unit C_g'' W L")
    hu = 1e-4
    dqg_dus = (tc.qg(1 + hu, 0.0, phit=0.0) - tc.qg(1 - hu, 0.0, phit=0.0)) / (2 * hu) / tc.QG_SCALE
    dqg_dud = (tc.qg(1.0, hu, phit=0.0) - tc.qg(1.0, -hu, phit=0.0)) / (2 * hu) / tc.QG_SCALE
    dqd_dus = (tc.qd(1 + hu, 0.0, phit=0.0) - tc.qd(1 - hu, 0.0, phit=0.0)) / (2 * hu) / tc.QG_SCALE
    print(f"  dQ_g/du_s = {dqg_dus:.6f} (2/3 = 0.666667)   dQ_g/du_d = {dqg_dud:.1e} (0)   "
          f"-dQ_d/du_s = {-dqd_dus:.6f} (4/15 = 0.266667)")

    print()
    print("Check 6: step size, C_gg and C_gd at V_GS = 0 V, V_DS = 10 V")
    for h in [2e-3, 1e-3, 5e-4]:
        C = cap_matrix(0.0, 10.0, h=h, **kw)
        print(f"  h = {h * 1e3:3.1f} mV   C_gg = {C[0, 0] * fF:.6f} fF   C_gd = {C[0, 1] * fF:.6f} fF")

    for VD in [10.0, 0.05]:
        print()
        print(f"Capacitances (fF) at V_DS = {VD:g} V, intrinsic")
        print("   V_GS     C_gg      C_gs      C_gd      C_dg      C_dd     C_dg - C_gd")
        for V in [-1.5, -1.0, -0.5, 0.0, 1.0, 2.0, 3.0]:
            C = cap_matrix(V, VD, **kw) * fF
            print(f"  {V:+5.2f}  {C[0, 0]:8.3f}  {C[0, 2]:8.3f}  {C[0, 1]:8.3f}  {C[1, 0]:8.3f}  {C[1, 1]:8.3f}    {C[1, 0] - C[0, 1]:+8.3f}")