"""
va_twin.py - Step 10a: a Python twin of the Verilog-A model va/graded_hemt.va, line for line, in SI units.
It exists to test the Verilog-A algorithm before it reaches Spectre:
  1. the charge solve (safeguarded Newton + one final Newton step) against Brent's method on the same equation;
  2. the twin's drain current against the validated Python model (fit.py / fit2.py), both devices;
  3. the explicit inverse of the access resistance against the form in access.py.
Run:  python3 va_twin.py            (reads params_device1.json and params_device2.json; each device in its own process)

What differs from the Python model, and why
  F_1/2: the twin (and the Verilog-A) use the paper's closed form, Eq. (2) with the 2/sqrt(pi) of Eq. (1) [P, Step 1],
         because the exact integral cannot be evaluated cheaply inside a circuit simulator. The Python model uses the
         exact integral. Effect on n: at most 0.33 % (checked 17 Nov) [CALC].
  Access resistance: access.py solves I = ids(V_GS - I R_S(I), V_DS - I (R_S + R_D)) with
         R(I) = R0 / [1 - (I/I_sat)^beta]^(1/beta)   [T, Eqs. 5.1-5.3].
         Its exact inverse is I(V) = (V/R0) / [1 + (|V|/(R0 I_sat))^beta]^(1/beta) [CALC]: the same element, written as
         a current controlled by its voltage, which is how the Verilog-A puts it between the outer and inner nodes.
  Negative V_DS: the Python model returns 0; the twin and the Verilog-A swap source and drain [ASSUME: standard
         source-drain interchange; needed for harmonic balance, where V_DS can swing below 0].
"""
import json
import os
import numpy as np
from scipy.optimize import brentq

# ---- constants (SI) [CODATA]
Q = 1.602176634e-19
KB = 1.380649e-23
HBAR = 1.054571817e-34
M0 = 9.1093837e-31
EPS0 = 8.8541878128e-12
A_AY, B_AY = 9.6, 2.13            # Eq. (2) constants [P]
GAM32 = np.sqrt(np.pi) / 2        # Gamma(3/2)
LN2SP = np.log(2 / np.sqrt(np.pi))


def card(dev):
    """Model card for one device, from device1.py / device2.py [P, TCAD] and params_deviceN.json [CALIBRATE]."""
    if dev == 1:
        import device1 as d
        p = json.load(open("params_device1.json"))
        gamma = None
    else:
        import device2 as d
        p = json.load(open("params_device2.json"))
        gamma = p["fitted [CALIBRATE]"]["gamma"]
    f, x = p["fitted [CALIBRATE]"], p["fixed"]
    return dict(
        W=d.W, L=d.L, EPSR=d.EPS_R, DB=d.D_ALGAN, TG=d.T_GRADE, MEFF=0.217, TEMP=300.0,
        VOFF=x["v_off"], GAMMA=(2 / 3 if gamma is None else gamma),
        MULF=x["mu_lf"], MUA=x["mu_a"], MUB=x["mu_b"], THSAT=x["theta_sat"],
        VSAT0=f["vsat0"], TH1=x["th1"], TH2=f["th2"], NMAX=x["n_max"], DELTA=f["delta"],
        LSG=d.L_SG, LGD=d.L_GD, NACC=d.N_ACC * 1e4, U0=f["u0"], VSATACC=x["vsat_acc"], BETA=f["beta"])


class Twin:
    def __init__(self, c):
        self.c = c
        self.phit = KB * c["TEMP"] / Q
        x = c["MEFF"] * M0 * KB * c["TEMP"] / (np.pi * HBAR**2)       # m^-2
        self.lnNC = np.log(0.25 * (2 * x) ** 1.5)                      # ln N_C, N_C in m^-3 [P, below Eq. (1)]
        self.cg = c["EPSR"] * EPS0 / (2 * c["DB"] + c["TG"])           # F/m^2, printed form [P]
        self.u = self.cg / Q                                           # m^-2 per V
        self.m = 2 / (3 * c["GAMMA"])                                  # [DEVIATION, Step 8c]
        self.teff = c["DB"] + c["TG"]
        qacc = Q * c["NACC"]
        self.r0s = c["LSG"] / (c["W"] * qacc * c["U0"])
        self.r0d = c["LGD"] / (c["W"] * qacc * c["U0"])
        self.isat = c["W"] * qacc * c["VSATACC"]

    # ---------- F_1/2 by Eq. (2), in logs so that eta from -1000 to +1000 never overflows
    @staticmethod
    def lnF_dlnF(eta):
        s = np.sqrt((eta - B_AY) ** 2 + A_AY)
        Bk = B_AY + eta + s
        dB = 1 + (eta - B_AY) / s
        hi = 1.5 * 2**1.5 * Bk ** -1.5
        dhi = -1.5 * hi * dB / Bk
        if eta < 0:
            e = np.exp(eta)
            den = hi * e + 1 / GAM32
            return LN2SP + eta - np.log(den), -(dhi * e - 1 / GAM32) / den
        lo = np.exp(-eta) / GAM32
        den = hi + lo
        return LN2SP - np.log(den), -(dhi - lo) / den

    def sheet(self, eta):
        """S = (N_C F_1/2)^(2/3) in m^-2 and dS/deta."""
        lnF, dl = self.lnF_dlnF(eta)
        S = np.exp(2 / 3 * (self.lnNC + lnF))
        return S, 2 / 3 * S * dl

    def eta_solve(self, Vg0, maxit=200, tol=1e-13):
        """Root of f(eta) = S(eta) - u (V_g0 - m phi_t eta): safeguarded Newton, then one final Newton step
        (the final step makes the simulator's derivatives exact whatever path the loop took)."""
        c = self.u * self.m * self.phit
        hi = Vg0 / (self.m * self.phit)
        S, _ = self.sheet(hi)
        lo = hi - S / c - 1.0
        eta, it, dn = hi, 0, 1.0
        while it < maxit and abs(dn) > tol * (1 + abs(eta)):
            S, dS = self.sheet(eta)
            f = S - self.u * Vg0 + c * eta
            if f > 0:
                hi = eta
            else:
                lo = eta
            new = eta - f / (dS + c)
            if not (lo < new < hi):
                new = 0.5 * (lo + hi)
            dn, eta, it = new - eta, new, it + 1
        S, dS = self.sheet(eta)
        eta = eta - (S - self.u * Vg0 + c * eta) / (dS + c)
        return eta, it

    def state(self, Vgs, Vds):
        """Intrinsic I (A) and charges (C) for V_DS >= 0 [P, Eqs. (17)-(19); current.py, terminal_charges.py]."""
        c = self.c
        Vg0 = Vgs - c["VOFF"]
        eta, _ = self.eta_solve(Vg0)
        lnF, _ = self.lnF_dlnF(eta)
        n = np.exp(max(2 / 3 * (self.lnNC + lnF), -400.0))              # m^-2; floor only below V_g0 ~ -20 V
        Vq = n / self.u                                                 # V, q n / C_g''
        vsat = c["VSAT0"] / (c["TH1"] + (n / c["NMAX"]) ** c["TH2"])
        Eys = Vq / self.teff
        mus = c["MULF"] / (1 + c["MUA"] * Eys + c["MUB"] * Eys**2)
        Vdsat = Vq / (1 + Vq * mus / (2 * vsat * c["L"]))
        dl = c["DELTA"]
        if Vds <= Vdsat:                                                # V_DS/(1+(V_DS/V_DSAT)^delta)^(1/delta),
            Vdeff = Vds / (1 + (Vds / Vdsat) ** dl) ** (1 / dl)          # written with the ratio <= 1 on each side
        else:                                                           # (same function; no overflow or 0/0)
            Vdeff = Vdsat / (1 + (Vdsat / Vds) ** dl) ** (1 / dl)
        um = Vq - 0.5 * Vdeff                                           # V_g0 - V_m
        Ey = um / self.teff
        mu = c["MULF"] / (1 + c["MUA"] * Ey + c["MUB"] * Ey**2)
        mueff = mu / np.sqrt(1 + (c["THSAT"] * Vdeff) ** 2)
        I = c["W"] / c["L"] * mueff * self.cg * (um + self.phit) * Vdeff
        us, ud, pt = Vq, Vq - Vdeff, self.phit
        k = self.cg * c["W"] * c["L"]
        Qg = k * ((us**2 + us * ud + ud**2) / 3 + pt * 0.5 * (us + ud)) / (0.5 * (us + ud) + pt)
        Qd = -k * (8 * us**3 + 16 * us**2 * ud + 24 * us * ud**2 + 12 * ud**3
                   + pt * (25 * us**2 + 50 * us * ud + 45 * ud**2) + pt**2 * (20 * us + 40 * ud)) / (30 * (us + ud + 2 * pt) ** 2)
        return I, Qg, Qd, -(Qg + Qd), n

    def ids(self, Vgs, Vds):
        if Vds >= 0:
            return self.state(Vgs, Vds)[0]
        return -self.state(Vgs - Vds, -Vds)[0]

    def r_of_i(self, I, r0):
        return r0 / (1 - (I / self.isat) ** self.c["BETA"]) ** (1 / self.c["BETA"])

    def i_of_v(self, V, r0):
        """Explicit inverse of V = I R(I): what the Verilog-A contributes between outer and inner nodes."""
        return (V / r0) / (1 + (abs(V) / (r0 * self.isat)) ** self.c["BETA"]) ** (1 / self.c["BETA"])

    def id_ext(self, Vgs, Vds):
        """External drain current (A), solved like access.py (V_DS > 0)."""
        def f(I):
            vs, vd = I * self.r_of_i(I, self.r0s), I * self.r_of_i(I, self.r0d)
            return I - self.ids(Vgs - vs, max(Vds - vs - vd, 0.0))
        hi = min(Vds / (self.r0s + self.r0d), self.isat * (1 - 1e-9))
        return brentq(f, 0.0, hi, xtol=1e-16, rtol=1e-13)


if __name__ == "__main__":
    import subprocess
    import sys
    if os.environ.get("GAN_DEVICE") is None:          # one process per device: device.py picks the device at import
        for dev in ("1", "2"):
            subprocess.run([sys.executable, __file__], env={**os.environ, "GAN_DEVICE": dev}, check=True)
        sys.exit(0)
    for dev in (int(os.environ["GAN_DEVICE"]),):
        tw = Twin(card(dev))
        print(f"===== Device {dev}: C_g'' = {tw.cg:.4e} F/m^2, V_OFF = {tw.c['VOFF']:+.3f} V, m = {tw.m:.4f}, "
              f"R_S0 = {tw.r0s:.3f} ohm, R_D0 = {tw.r0d:.3f} ohm, I_sat = {tw.isat * 1e3:.1f} mA")

        # 1. charge solve against Brent's method, V_GS from -20 to +6 V
        worst, itmax = 0.0, 0
        for V in np.arange(-20.0, 6.0 + 1e-9, 0.05):
            Vg0 = V - tw.c["VOFF"]
            eta, it = tw.eta_solve(Vg0)
            g = lambda e: tw.sheet(e)[0] - tw.u * (Vg0 - tw.m * tw.phit * e)
            ref = brentq(g, -2000.0, Vg0 / (tw.m * tw.phit) + 1e-9, xtol=1e-15, rtol=1e-15)
            worst = max(worst, abs(eta - ref) / (1 + abs(ref)))
            itmax = max(itmax, it)
        print(f"  1. eta solve vs Brent, V_GS -20..6 V: worst relative error {worst:.1e}, at most {itmax} iterations")

        # 2. drain current against the validated Python model
        if dev == 1:
            import fit as F
            p = list(json.load(open("params_device1.json"))["fitted [CALIBRATE]"].values())
            ref = lambda vg, vd: F.model_id(p, vg, vd)
            grid = [(vg, 10.0) for vg in np.arange(-1.6, 3.01, 0.2)] + [(vg, vd) for vg in (-1, 0, 1, 2, 3) for vd in (0.5, 2.0, 5.0)]
        else:
            import fit2 as F
            p = list(json.load(open("params_device2.json"))["fitted [CALIBRATE]"].values())
            ref = lambda vg, vd: F.model_id(p, vg, vd)
            grid = [(vg, 5.0) for vg in np.arange(-1.5, 1.01, 0.1)] + [(vg, vd) for vg in (-1, -0.5, 0, 1) for vd in (0.5, 2.0, 5.0)]
        rows = [(vg, vd, tw.id_ext(vg, vd) * 1e3, ref(vg, vd)) for vg, vd in grid]
        on = [r for r in rows if r[3] > 0.05 * max(x[3] for x in rows)]
        dev_pct = max(abs(r[2] / r[3] - 1) for r in on) * 100
        print(f"  2. I_D twin vs Python model, {len(rows)} biases: worst {dev_pct:.2f} % where I_D > 5 % of max;"
              f" worst absolute {max(abs(r[2] - r[3]) for r in rows) * 1e3:.1f} uA")
        for r in rows[:: max(1, len(rows) // 6)]:
            print(f"     V_GS {r[0]:+5.2f} V_DS {r[1]:4.1f}:  twin {r[2]:9.4f} mA   Python {r[3]:9.4f} mA")

        # 3. explicit inverse of the access resistance
        e = max(abs(tw.i_of_v(I * tw.r_of_i(I, tw.r0d), tw.r0d) / I - 1) for I in np.linspace(1e-6, 0.999 * tw.isat, 400))
        print(f"  3. I(V(I)) round trip over 0..I_sat: worst relative error {e:.1e}")

        # 4. source-drain symmetry at V_DS = 0 (needed for the swap)
        vg = 0.0
        g1 = (tw.ids(vg, 1e-6) - tw.ids(vg, -1e-6)) / 2e-6
        print(f"  4. source-drain swap at V_GS = 0: channel conductance at V_DS = 0 {g1 * 1e3:.3f} mS, continuous;"
              f" I(+1 mV) = {tw.ids(vg, 1e-3) * 1e6:.3f} uA, I(-1 mV) = {tw.ids(vg, -1e-3) * 1e6:.3f} uA")
