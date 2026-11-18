"""
fig7c.py - Step 12a: Fig. 7(c) of Nazir et al. 2025, simulation study of theta_SAT2 on G_M and V_SAT.
Run:  python3 fig7c.py      -> fig7c_device2.png
Device 2 Step 8c model (params_device2.json, no DIBL: the paper's model has none), with the printed study values
theta_SAT1 = 0.3, V_SAT0 = 1e5 m/s, n_max = 5e16 m^-2 [P, Fig. 7(c)] and theta_SAT2 = 0.1 ... 0.5 [P].
Bias V_DS = 5 V [ASSUME: not printed; the Fig. 6(b) bias]. G_M per mm (as Figs. 6(b), 6(c); printed item 11).
V_OFF = -0.63 V [TCAD n8, low V_DS; ASSUME]: with the Fig. 6(b) value (-1.01 V) every curve turns on 0.44-0.45 V
earlier than the paper's (at G_M = 150 mS/mm); with -0.63 V the peak sits at -0.18 V (paper -0.20 V). The paper's study
therefore used V_OFF ~ -0.6 V, not the Fig. 6(b) value [CALC, Step 12a].
V_SAT = V_SAT0 / (theta_SAT1 + (n/n_max)^theta_SAT2) [P], with n the channel density at the source end of the
intrinsic device. Model: the Verilog-A twin (va_twin.py) with DIBL off.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from va_twin import Twin, card

VDS, W_MM = 5.0, 0.075
TH2 = (0.1, 0.2, 0.3, 0.4, 0.5)
COL = {0.1: "#a00000", 0.2: "#22456a", 0.3: "#b68901", 0.4: "#27834b", 0.5: "#454567"}
V = np.arange(-3.0, 1.0 + 1e-9, 0.02)
L = lambda f: np.loadtxt(f, delimiter=",", skiprows=3)


VOFF = -0.63


def study(th2, voff=VOFF):
    tw = Twin({**card(2), "SIGMA": 0.0, "TH1": 0.3, "VSAT0": 1e5, "NMAX": 5e16, "TH2": th2, "VOFF": voff})
    I = np.array([tw.id_ext(v, VDS) for v in V])
    gm = np.gradient(I, V) / W_MM * 1e3                                   # mS/mm
    vs = []
    for v, i in zip(V, I):
        vsi = i * tw.r_of_i(i, tw.r0s) if i > 0 else 0.0
        vdi = VDS - (i * tw.r_of_i(i, tw.r0d) if i > 0 else 0.0)
        n = tw.state(v - vsi, vdi - vsi)[4]
        vs.append(1e5 / (0.3 + (n / 5e16) ** th2) / 1e3)                  # km/s
    return gm, np.array(vs)


if __name__ == "__main__":
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax2 = ax.twinx()
    print("theta_SAT2   G_M peak (mS/mm) at V_GS     paper            V_SAT at -3 V / +1 V (km/s)   paper")
    for th in TH2:
        gm, vs = study(th)
        pg, pv = L(f"paper_data/fig7c_paper_model_GM_th2_{th:.1f}.csv"), L(f"paper_data/fig7c_paper_model_VSAT_th2_{th:.1f}.csv")
        k = int(np.argmax(gm)); kp = int(np.argmax(pg[:, 1]))
        g1, _ = study(th, -1.01)
        half = lambda g, v: np.interp(150, g[:np.argmax(g)], v[:np.argmax(g)])
        rms = np.sqrt(np.mean((np.interp(pg[:, 0], V, gm) - pg[:, 1]) ** 2))
        print(f"      G_M = 150 mS/mm at {half(gm, V):+.3f} V (paper {half(pg[:, 1], pg[:, 0]):+.3f}; with V_OFF -1.01: "
              f"{half(g1, V):+.3f}); rms vs paper line {rms:.1f} mS/mm")
        print(f"   {th:.1f}       {gm[k]:6.1f} at {V[k]:+.2f}       {pg[kp, 1]:6.1f} at {pg[kp, 0]:+.2f}      "
              f"{vs[0]:5.0f} / {vs[-1]:4.0f}               {pv[0, 1]:5.0f} / {pv[-1, 1]:4.0f}")
        ax.plot(*pg.T, "-", c=COL[th], lw=1, alpha=0.5)
        ax.plot(V, gm, "--", c=COL[th], lw=2, label=f"theta_SAT2 = {th}")
        ax2.plot(*pv.T, "-", c=COL[th], lw=1, alpha=0.5)
        ax2.plot(V, vs, ":", c=COL[th], lw=2)
    ax.set(xlabel="V_GS (V)", ylabel="G_M (mS/mm)", xlim=(-3, 1), ylim=(0, 500))
    ax2.set(ylabel="V_SAT (km/s)", ylim=(0, 500))
    ax.set_title("Fig. 7(c), Device 2, V_DS = 5 V, V_OFF = -0.63 V: paper (thin), this work (G_M dashed, V_SAT dotted)", loc="left", fontsize=9)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.grid(c="#e1e0d9")
    plt.tight_layout()
    plt.savefig("fig7c_device2.png", dpi=150)
    print("saved fig7c_device2.png")
