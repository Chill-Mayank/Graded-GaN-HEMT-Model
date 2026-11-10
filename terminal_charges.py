"""
terminal_charges.py - compact model step 7: terminal charges of Device 1, intrinsic (no access resistance).
  Step 7a: gate charge Q_g    [P: Nazir et al., IEEE TED 72(9), 2025, Sec. II-C, Eqs. (18)-(19)]
  Step 7b: drain charge Q_d by Ward-Dutton partition [D. E. Ward and R. W. Dutton, IEEE JSSC 13(5), 703-708, 1978],
           "akin to [24]" [P]; [24] = Khandelwal, Chauhan, Fjeldly, IEEE TED 59(10), 2012, Eq. (13) as corrected
           in spec D.6 (printed [24] Eqs. (12)-(13) contain errors: spec P15). Source charge Q_s = -(Q_g + Q_d) [P].

Channel potentials: V_ch,s = E_F and V_ch,d = E_F + V_DS,eff, from current.bias_state().
  Eq. (17) and Sec. II-C both define V_m: (2 E_F + V_DS,eff)/2 there, (V_ch,s + V_ch,d)/2 here, so the two
  readings coincide [P]. Written from the gate, u = V_g0 - V_ch:
  u_s = V_g0 - V_ch,s = q n / C_g''   and   u_d = V_g0 - V_ch,d = u_s - V_DS,eff.
  [24]'s psi_s, psi_d are V_ch,s, V_ch,d here (notation map in spec D.6).

Checks run by __main__:
 Step 7a (Q_g)
  1. Eq. (19) as printed against the same expression written in u_s, u_d (identical algebra [CALC]). The printed
     form subtracts terms of size V_g0^2, so it loses digits once q n / C_g'' is small (subthreshold).
  2. Eq. (18) integrated along the channel: u(x) from the substitution printed under Eq. (18) (current continuity,
     (V_g0 - V_ch + phi_t) dV_ch/dx = constant), inverted by Brent's method [R. P. Brent, 1973] (scipy brentq);
     the x integral by adaptive quadrature [R. Piessens et al., QUADPACK, 1983] (scipy quad).
  3. The substitution as printed divides by V_DS; it maps the channel onto length L only with V_ch,d - V_ch,s
     (= V_DS,eff here) in its place [spec P13]. Eq. (19) itself follows from the corrected substitution.
  4. V_DS = 0: Eq. (19) must reduce to Q_g = q W L n, the gate charge balancing the channel electrons [CALC].
 Step 7b (Q_d, Q_s)
  5. Spec D.6 form (in V_g0, psi_s, psi_d) against the same expression written in u_s, u_d (identical algebra,
     checked symbolically with sympy [CALC]); the u form has only positive terms, so no digits are lost.
  6. Q_d and Q_s against the Ward-Dutton integrals along the channel, Q_d = -W C_g'' int (x/L) u dx and
     Q_s = -W C_g'' int (1 - x/L) u dx, same u(x) and quadrature as check 2; Q_s also tests Q_s = -(Q_g + Q_d).
  7. Limits: V_DS = 0 gives Q_d = Q_s = -Q_g/2; with phi_t -> 0 and u_d = 0 (long-channel pinch-off) the
     partition is the Ward-Dutton 40/60: Q_d = -0.4 Q_g, Q_s = -0.6 Q_g [Ward and Dutton 1978].
"""
import inspect
import json
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
import current as cur
from gate import Q, PHIT

QG_SCALE = cur.CG * cur.L * cur.W        # F, C_g'' L W


def model_params(path="params_device1.json"):
    """Step 5e parameters that the intrinsic model uses (access-region values u0, beta, vsat_acc, ... are left out)."""
    p = json.load(open(path))
    allp = {**p["fixed"], **p["fitted [CALIBRATE]"]}
    keep = inspect.signature(cur.bias_state).parameters
    return {k: v for k, v in allp.items() if k in keep}


def channel_potentials(V_GS, V_DS, **kw):
    """Channel ends at one intrinsic bias, as a dict (all in V): Vg0, Vs = V_ch,s, Vd = V_ch,d,
    us = V_g0 - V_ch,s, ud = V_g0 - V_ch,d, Vdeff; plus 'state', the full bias_state() dict."""
    s = cur.bias_state(V_GS, V_DS, **kw)
    us = s["Vq"]                                           # q n / C_g'', taken directly (no V_g0 - E_F subtraction)
    return {"Vg0": s["Vg0"], "Vs": s["EF"], "Vd": s["EF"] + s["Vdeff"], "us": us, "ud": us - s["Vdeff"],
            "Vdeff": s["Vdeff"], "state": s}


# ---------------------------------------------------------------- Step 7a: gate charge
def qg_printed(Vg0, Vs, Vd):
    """Gate charge in C, Eq. (19) exactly as printed [P]."""
    Vm = 0.5 * (Vs + Vd)
    bracket = Vg0**2 + (Vd**2 + Vs**2 + Vd * Vs) / 3 - Vg0 * (Vd + Vs - PHIT) - PHIT * Vm
    return QG_SCALE / (Vg0 - Vm + PHIT) * bracket


def qg(us, ud, phit=PHIT):
    """Gate charge in C, Eq. (19) [P] written in u_s and u_d:
    bracket = (u_s^2 + u_s u_d + u_d^2)/3 + phi_t (u_s + u_d)/2,  V_g0 - V_m + phi_t = (u_s + u_d)/2 + phi_t.
    Same algebra as printed [CALC]; nothing large is subtracted, so it keeps full precision in subthreshold.
    This is the form the model uses."""
    um = 0.5 * (us + ud)                                   # V, V_g0 - V_m
    return QG_SCALE * ((us**2 + us * ud + ud**2) / 3 + phit * um) / (um + phit)


# ---------------------------------------------------------------- Step 7b: drain and source charges
def qd_d6(Vg0, Vs, Vd):
    """Drain charge in C, [24] Eq. (13) as corrected in spec D.6 (psi_s = V_ch,s, psi_d = V_ch,d, V_go = V_g0,
    V_th = phi_t), in the variables as written there."""
    ps, pd, p, V0 = Vs, Vd, PHIT, Vg0
    Vm = 0.5 * (ps + pd)
    br = (12 * pd**3 + 8 * ps**3 + ps**2 * (16 * pd - 5 * (5 * p + 8 * V0))
          + 2 * ps * (12 * pd**2 - 5 * pd * (5 * p + 8 * V0) + 10 * (p + V0) * (p + 4 * V0))
          - 15 * pd**2 * (3 * p + 4 * V0) - 60 * V0 * (p + V0)**2 + 20 * pd * (p + V0) * (2 * p + 5 * V0))
    return QG_SCALE * br / (120 * (V0 - Vm + p)**2)


def qd(us, ud, phit=PHIT):
    """Drain charge in C, the spec D.6 form written in u_s and u_d (same algebra, checked with sympy [CALC]):
    Q_d = -C_g'' L W [8 u_s^3 + 16 u_s^2 u_d + 24 u_s u_d^2 + 12 u_d^3 + phi_t (25 u_s^2 + 50 u_s u_d + 45 u_d^2)
          + phi_t^2 (20 u_s + 40 u_d)] / [30 (u_s + u_d + 2 phi_t)^2]
    Every term is positive, so Q_d < 0 at every bias and no digits are lost. This is the form the model uses."""
    num = (8 * us**3 + 16 * us**2 * ud + 24 * us * ud**2 + 12 * ud**3
           + phit * (25 * us**2 + 50 * us * ud + 45 * ud**2) + phit**2 * (20 * us + 40 * ud))
    return -QG_SCALE * num / (30 * (us + ud + 2 * phit)**2)


def qs(us, ud, phit=PHIT):
    """Source charge in C, Q_s = -(Q_g + Q_d) [P, Sec. II-C]."""
    return -(qg(us, ud, phit) + qd(us, ud, phit))


# ---------------------------------------------------------------- integrals along the channel (checks)
def x_of_u(u, us, ud, denom):
    """Position x (m) where V_g0 - V_ch = u, from the substitution under Eq. (18) integrated from the source:
    x = L * integral (u' + phi_t) du' from u to u_s / [((u_s + u_d)/2 + phi_t) * denom].
    Printed: denom = V_DS.  Corrected [spec P13]: denom = V_ch,d - V_ch,s = u_s - u_d."""
    F = 0.5 * (us**2 - u**2) + PHIT * (us - u)             # V^2, the integral of a straight line, exact
    return cur.L * F / ((0.5 * (us + ud) + PHIT) * denom)


def charges_along_x(us, ud):
    """(Q_g, Q_d, Q_s) in C from the integrals along the channel: u(x) from the corrected x(u) by Brent's method,
    each x integral by QUADPACK. Q_g = W C int u dx (Eq. (18)); Ward-Dutton: Q_d = -W C int (x/L) u dx,
    Q_s = -W C int (1 - x/L) u dx."""
    def u_at(x):
        return brentq(lambda u: x_of_u(u, us, ud, us - ud) - x, ud, us, xtol=1e-18, rtol=1e-15)
    opts = dict(epsabs=0.0, epsrel=1e-13, limit=200)
    L, WC = cur.L, cur.W * cur.CG
    q_g = WC * quad(u_at, 0.0, L, **opts)[0]
    q_d = -WC * quad(lambda x: x / L * u_at(x), 0.0, L, **opts)[0]
    q_s = -WC * quad(lambda x: (1 - x / L) * u_at(x), 0.0, L, **opts)[0]
    return q_g, q_d, q_s


def qg_along_x(us, ud):
    """Gate charge in C from Eq. (18) itself (Step 7a check)."""
    return charges_along_x(us, ud)[0]


if __name__ == "__main__":
    kw = model_params()
    print(f"Device 1: W = {cur.W * 1e6:.0f} um, L = {cur.L * 1e6:.1f} um, C_g'' = {cur.CG:.4e} F/m^2, "
          f"C_g'' W L = {QG_SCALE * 1e15:.2f} fF, V_OFF = {kw['v_off']:+.3f} V; Step 5e parameters")

    print()
    print("========== Step 7a: gate charge Q_g ==========")
    print("Check 1: Eq. (19) as printed against the u_s, u_d form (V_DS = 1 V)")
    print("   V_GS     n (cm^-2)    Q_g u-form (fC)    Q_g printed (fC)   relative difference")
    for V in [-3.0, -2.0, -1.5, -1.0, 0.0, 3.0]:
        c = channel_potentials(V, 1.0, **kw)
        a, b = qg_printed(c["Vg0"], c["Vs"], c["Vd"]), qg(c["us"], c["ud"])
        print(f"  {V:+5.2f}    {c['state']['n']:10.3e}    {b * 1e15:13.6g}      {a * 1e15:13.6g}        {abs(a / b - 1):9.1e}")

    print()
    print("Check 2: Eq. (19) against Eq. (18) integrated along the channel")
    print("   V_GS   V_DS (V)   V_DS,eff (V)    Q_g Eq. (19) (fC)   Q_g along x (fC)   relative difference")
    for V, VD in [(-1.5, 1.0), (-1.0, 0.5), (0.0, 0.1), (0.0, 10.0), (1.0, 2.0), (3.0, 10.0)]:
        c = channel_potentials(V, VD, **kw)
        a, b = qg(c["us"], c["ud"]), qg_along_x(c["us"], c["ud"])
        print(f"  {V:+5.2f}   {VD:6.2f}     {c['Vdeff']:9.4f}      {a * 1e15:15.6f}    {b * 1e15:15.6f}       {abs(a / b - 1):9.1e}")

    print()
    print("Check 3: where the substitution under Eq. (18) puts the drain end, x(u_d) / L")
    print("   V_GS   V_DS (V)   V_DS,eff (V)   printed (V_DS)   corrected (V_ch,d - V_ch,s)")
    for V, VD in [(0.0, 0.1), (0.0, 10.0), (3.0, 10.0)]:
        c = channel_potentials(V, VD, **kw)
        r_p = x_of_u(c["ud"], c["us"], c["ud"], VD) / cur.L
        r_c = x_of_u(c["ud"], c["us"], c["ud"], c["us"] - c["ud"]) / cur.L
        print(f"  {V:+5.2f}   {VD:6.2f}     {c['Vdeff']:9.4f}        {r_p:9.4f}          {r_c:9.4f}")

    print()
    print("Check 4: V_DS = 0, Eq. (19) against q W L n")
    for V in [-3.0, -1.5, 0.0, 3.0]:
        c = channel_potentials(V, 0.0, **kw)
        q_g = qg(c["us"], c["ud"])
        qwln = Q * cur.W * cur.L * c["state"]["n"] * 1e4
        print(f"  V_GS = {V:+5.2f} V   Q_g = {q_g * 1e15:12.6g} fC   q W L n = {qwln * 1e15:12.6g} fC   "
              f"relative difference = {abs(q_g / qwln - 1):.1e}")

    print()
    print("Q_g (fC), intrinsic bias")
    print("   V_GS    V_DS = 0      1 V       10 V     Q_g(10 V)/Q_g(0)")
    for V in [-1.5, -1.0, -0.5, 0.0, 1.0, 2.0, 3.0]:
        q = []
        for VD in [0.0, 1.0, 10.0]:
            c = channel_potentials(V, VD, **kw)
            q.append(qg(c["us"], c["ud"]) * 1e15)
        print(f"  {V:+5.2f}  {q[0]:9.4g}  {q[1]:9.4g}  {q[2]:9.4g}        {q[2] / q[0]:.3f}")

    print()
    print("========== Step 7b: drain and source charges Q_d, Q_s ==========")
    print("Check 5: spec D.6 form against the u_s, u_d form (V_DS = 1 V)")
    print("   V_GS     n (cm^-2)    Q_d u-form (fC)    Q_d D.6 form (fC)   relative difference")
    for V in [-3.0, -2.0, -1.5, -1.0, 0.0, 3.0]:
        c = channel_potentials(V, 1.0, **kw)
        a, b = qd_d6(c["Vg0"], c["Vs"], c["Vd"]), qd(c["us"], c["ud"])
        print(f"  {V:+5.2f}    {c['state']['n']:10.3e}    {b * 1e15:13.6g}     {a * 1e15:13.6g}         {abs(a / b - 1):9.1e}")

    print()
    print("Check 6: Q_d and Q_s against the Ward-Dutton integrals along the channel")
    print("   V_GS   V_DS (V)    Q_d (fC)    along x (fC)   rel. diff.    Q_s (fC)    along x (fC)   rel. diff.")
    for V, VD in [(-1.5, 1.0), (-1.0, 0.5), (0.0, 0.1), (0.0, 10.0), (1.0, 2.0), (3.0, 10.0)]:
        c = channel_potentials(V, VD, **kw)
        _, xd, xs = charges_along_x(c["us"], c["ud"])
        d, s = qd(c["us"], c["ud"]), qs(c["us"], c["ud"])
        print(f"  {V:+5.2f}   {VD:6.2f}   {d * 1e15:10.4f}   {xd * 1e15:10.4f}    {abs(d / xd - 1):8.1e}   "
              f"{s * 1e15:10.4f}   {xs * 1e15:10.4f}    {abs(s / xs - 1):8.1e}")

    print()
    print("Check 7: limits of the partition, Q_d/Q_g and Q_s/Q_g")
    for V in [-1.5, 0.0, 3.0]:
        c = channel_potentials(V, 0.0, **kw)
        g = qg(c["us"], c["ud"])
        print(f"  V_DS = 0, V_GS = {V:+5.2f} V:   Q_d/Q_g = {qd(c['us'], c['ud']) / g:+.6f}   "
              f"Q_s/Q_g = {qs(c['us'], c['ud']) / g:+.6f}   (Ward-Dutton: -0.5, -0.5)")
    g0 = qg(1.0, 0.0, phit=0.0)
    print(f"  phi_t = 0, u_s = 1 V, u_d = 0:  Q_d/Q_g = {qd(1.0, 0.0, phit=0.0) / g0:+.6f}   "
          f"Q_s/Q_g = {qs(1.0, 0.0, phit=0.0) / g0:+.6f}   (Ward-Dutton 40/60: -0.4, -0.6)")

    print()
    print("Q_g, Q_d, Q_s (fC) at V_DS = 10 V, intrinsic bias")
    print("   V_GS     Q_g        Q_d        Q_s      Q_d/Q_g")
    for V in [-1.5, -1.0, -0.5, 0.0, 1.0, 2.0, 3.0]:
        c = channel_potentials(V, 10.0, **kw)
        g, d, s = qg(c["us"], c["ud"]), qd(c["us"], c["ud"]), qs(c["us"], c["ud"])
        print(f"  {V:+5.2f}  {g * 1e15:9.4g}  {d * 1e15:9.4g}  {s * 1e15:9.4g}    {d / g:+.3f}")