"""
rf_twin.py - Step 10c: small-signal two-port of the Verilog-A model (graded_hemt.va), computed from va_twin.py.
Same circuit as the Verilog-A: R_G between g and gi; access elements d-di and si-s (their differential conductance
at the bias point); the channel current and the terminal charges Q_g, Q_d, Q_s on gi, di, si; constant fringing
capacitors gi-si, gi-di, di-si. Common source, source grounded.
The nodal admittance matrix is reduced to the (g, d) two-port, then [standard two-port relations]:
  H21 = Y21 / Y11
  U   = |Y21 - Y12|^2 / (4 (Re Y11 Re Y22 - Re Y12 Re Y21))       [S. J. Mason, IRE Trans. CT-1, 20-25, 1954]
  f_T, f_MAX: either where |H21| = 1 and U = 1 on the model curve ("cross"), or, as the paper does for Fig. 11(b),
  extrapolated at -20 dB/decade from one frequency f0: f_T = f0 |H21(f0)|, f_MAX = f0 sqrt(U(f0)) ("extrap").
  The paper's own Fig. 11(b) lines equal the extrapolation of its Fig. 11(a) lines (21.7 / 44.3 GHz for any
  f0 = 1-5 GHz), while its U line crosses 0 dB near 60 GHz [CALC, Step 10c].
It is used to calibrate R_G and the fringing capacitances; Spectre then confirms with the same values.
"""
import numpy as np
from va_twin import Twin, card

N = {"g": 0, "d": 1, "gi": 2, "di": 3, "si": 4}       # node s is ground


class RF:
    def __init__(self, dev, rg=0.0, cgs=0.0, cgd=0.0, cds=0.0, over=None):
        self.tw = Twin({**card(dev), **(over or {})})
        self.rg, self.cgs, self.cgd, self.cds = rg, cgs, cgd, cds

    def op(self, vgs, vds):
        """DC operating point: drain current and the intrinsic node voltages."""
        tw = self.tw
        I = tw.id_ext(vgs, vds)
        vsi = I * tw.r_of_i(I, tw.r0s) if I > 0 else 0.0
        vdi = vds - (I * tw.r_of_i(I, tw.r0d) if I > 0 else 0.0)
        return I, vgs, vdi, vsi

    def small_signal(self, vgs, vds, h=1e-4):
        """Conductance matrix G and capacitance matrix C on nodes (g, d, gi, di, si) at the bias point."""
        tw = self.tw
        I, vgi, vdi, vsi = self.op(vgs, vds)
        f = lambda a, b: np.array(tw.state(a, b)[:3])            # I, Q_g, Q_d at intrinsic V_GS, V_DS
        x = (vgi - vsi, vdi - vsi)
        dA = (f(x[0] + h, x[1]) - f(x[0] - h, x[1])) / (2 * h)   # d/dV_GS,i
        dB = (f(x[0], x[1] + h) - f(x[0], x[1] - h)) / (2 * h)   # d/dV_DS,i
        # derivative w.r.t. node voltages (gi, di, si): V_GS,i = gi - si, V_DS,i = di - si
        dn = lambda k: {"gi": dA[k], "di": dB[k], "si": -(dA[k] + dB[k])}
        G = np.zeros((5, 5)); C = np.zeros((5, 5))
        for col, val in dn(0).items():                            # channel current, di -> si
            G[N["di"], N[col]] += val; G[N["si"], N[col]] -= val
        for row, k in (("gi", 1), ("di", 2)):                     # charges on gi and di; si takes minus the sum
            for col, val in dn(k).items():
                C[N[row], N[col]] += val; C[N["si"], N[col]] -= val
        gacc = lambda V, r0: (tw.i_of_v(V + 1e-6, r0) - tw.i_of_v(V - 1e-6, r0)) / 2e-6
        for a, b, g in (("d", "di", gacc(vds - vdi, tw.r0d)), ("si", None, gacc(vsi, tw.r0s))):
            self._stamp(G, a, b, g)
        if self.rg > 0:
            self._stamp(G, "g", "gi", 1 / self.rg)
        else:
            self._stamp(G, "g", "gi", 1e6)                        # short
        for a, b, c in (("gi", "si", self.cgs), ("gi", "di", self.cgd), ("di", "si", self.cds)):
            self._stamp(C, a, b, c)
        return G, C, I

    @staticmethod
    def _stamp(M, a, b, y):
        M[N[a], N[a]] += y
        if b is not None:
            M[N[b], N[b]] += y; M[N[a], N[b]] -= y; M[N[b], N[a]] -= y

    def y2(self, freqs, vgs, vds):
        """Two-port Y (g, d) at each frequency; returns (Y array [nf, 2, 2], I_D)."""
        G, C, I = self.small_signal(vgs, vds)
        out = []
        for fr in np.atleast_1d(freqs):
            Y = G + 2j * np.pi * fr * C
            A, B, Cc, D = Y[:2, :2], Y[:2, 2:], Y[2:, :2], Y[2:, 2:]
            out.append(A - B @ np.linalg.solve(D, Cc))
        return np.array(out), I

    @staticmethod
    def gains(Y):
        y11, y12, y21, y22 = Y[:, 0, 0], Y[:, 0, 1], Y[:, 1, 0], Y[:, 1, 1]
        h21 = np.abs(y21 / y11)
        U = np.abs(y21 - y12) ** 2 / (4 * (y11.real * y22.real - y12.real * y21.real))
        return h21, U

    def ft_fmax(self, vgs, vds, mode="cross", f0=2e9):
        if mode == "extrap":
            Y, I = self.y2([f0], vgs, vds)
            h, U = self.gains(Y)
            return float(f0 * h[0]), float(f0 * np.sqrt(U[0])), I
        f = np.logspace(8, 12.5, 400)
        Y, I = self.y2(f, vgs, vds)
        h, U = self.gains(Y)
        cross = lambda g: float(np.exp(np.interp(0.0, -np.log(g), np.log(f)))) if g.min() < 1 < g.max() else np.nan  # g falls with f
        return cross(h), cross(U), I
