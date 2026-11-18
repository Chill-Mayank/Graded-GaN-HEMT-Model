"""
fig10.py - Step 10d: Device 1 large signal at 10 GHz against Fig. 10 of Nazir et al. 2025, from Spectre hb.
Run:  python3 fig10.py      (reads va/hb_device1.raw and va/hb2_device1.raw, Spectre nutascii) -> fig10_device1.png
Circuit: va/hb_device1.scs, va/hb2_device1.scs (va_netlists.py). Nothing is fitted here: R_G and C_GD,Fr are the
Fig. 11(a) values, the bias is V_DS = 15 V, I_DS = 25 mA [P].
Definitions [standard; the paper does not define them]:
  P_DEL = power delivered into the gate network at f0 = -0.5 Re(V_p1 I_P1*)  (the port's source impedance and the
          series inductor are lossless apart from the port resistance, so this is the power entering the gate)
  P_OUT = |V_p2|^2 / (2 R_L) at f0;  Gain = P_OUT / P_DEL;  PAE = (P_OUT - P_DEL) / P_DC, P_DC = V_DS I_D(dc)
  two tones (per tone): P_FUND at 10 GHz, IM3 at 2 f1 - f2 = 9.99 GHz, OIP3 = P_FUND + (P_FUND - IM3)/2,
  OIP3/P_DC in dB = OIP3 (dBm) - P_DC (dBm)
"""
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

VDS = 15.0
dbm = lambda w: 10 * np.log10(np.asarray(w) * 1e3)


def hb_points(fn):
    """One dict of spectra per sweep point (the complex plot that carries the harmonic indices)."""
    txt = open(fn).read()
    pts = []
    for blk in re.split(r"(?m)^Plotname:", txt)[1:]:
        head_all = blk.split("Values:", 1)[0]
        if "complex" not in head_all or "harm_index" not in head_all:
            continue
        nv = int(re.search(r"No\. Variables:\s*(\d+)", head_all).group(1))
        head = re.split(r"(?m)^Variables:", head_all, maxsplit=1)[1]
        names = [ln.split()[1] for ln in head.strip().splitlines() if ln.split()][:nv]
        tok = re.split(r"(?m)^(?:Title|Date):", blk.split("Values:", 1)[1])[0].split()
        rows, i = [], 0
        while i + nv < len(tok) + 1 and i < len(tok):
            i += 1
            row = []
            for _ in range(nv):
                a, b = tok[i].split(",")
                row.append(complex(float(a), float(b)))
                i += 1
            rows.append(row)
        v = np.array(rows)
        pts.append({n: v[:, j] for j, n in enumerate(names)})
    return pts


def at(p, q, fr):
    return p[q][int(np.argmin(np.abs(p["freq"].real - fr)))]


def one_tone(fn="va/hb_device1.raw", rl=45.0):
    out = []
    for p in hb_points(fn):
        pdel = -0.5 * np.real(at(p, "p1", 10e9) * np.conj(at(p, "P1:p", 10e9)))
        pout = abs(at(p, "p2", 10e9)) ** 2 / (2 * rl)
        pdc = VDS * abs(at(p, "VD:p", 0.0).real)
        out.append((dbm(pdel), dbm(pout), dbm(pout) - dbm(pdel), 100 * (pout - pdel) / pdc))
    return np.array(sorted(out))


def two_tone(fn="va/hb2_device1.raw", rl=50.0):
    out = []
    for p in hb_points(fn):
        pdel = -0.5 * np.real(at(p, "p1", 10e9) * np.conj(at(p, "P1:p", 10e9)))
        pf = dbm(abs(at(p, "p2", 10e9)) ** 2 / (2 * rl))
        im3 = dbm(abs(at(p, "p2", 9.99e9)) ** 2 / (2 * rl))
        pdc = dbm(VDS * abs(at(p, "VD:p", 0.0).real))
        oip3 = pf + (pf - im3) / 2
        out.append((dbm(pdel), pf, im3, oip3 - pdc))
    return np.array(sorted(out))


if __name__ == "__main__":
    L = lambda f: np.loadtxt(f, delimiter=",", skiprows=3, ndmin=2)
    A, B = one_tone(), two_tone()
    rms = lambda ours, m: np.sqrt(np.nanmean((np.interp(m[:, 0], ours[:, 0], ours[:, 1], left=np.nan, right=np.nan) - m[:, 1]) ** 2))
    print(f"Fig. 10(a): {len(A)} power points, P_DEL {A[0, 0]:.1f} .. {A[-1, 0]:.1f} dBm")
    print("  rms against measured        this work   paper's model")
    for k, q in ((1, "Pout_dBm"), (2, "Gain_dB"), (3, "PAE_pct")):
        m, pm = L(f"paper_data/fig10a_measured_{q}.csv"), L(f"paper_data/fig10a_paper_model_{q}.csv")
        print(f"  {q:<12}               {rms(A[:, [0, k]], m):8.2f}   {rms(pm, m):10.2f}")
    print(f"Fig. 10(b): {len(B)} power points, P_DEL {B[0, 0]:.1f} .. {B[-1, 0]:.1f} dBm per tone")
    for k, q in ((1, "Pfund_dBm"), (2, "IM3_dBm"), (3, "OIP3_over_PDC_dB")):
        m, pm = L(f"paper_data/fig10b_measured_{q}.csv"), L(f"paper_data/fig10b_paper_model_{q}.csv")
        print(f"  {q:<18}         {rms(B[:, [0, k]], m):8.2f}   {rms(pm, m):10.2f}")
    if len(B) > 2:
        s = np.polyfit(B[:, 0], B[:, 2], 1)[0]
        print(f"  IM3 slope {s:.2f} dB/dB (3 expected for a third-order product)")

    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.4))
    a2 = ax[0].twinx()
    for k, q, c, axx in ((1, "Pout_dBm", "#2a78d6", ax[0]), (2, "Gain_dB", "#a94b17", ax[0]), (3, "PAE_pct", "#222222", a2)):
        m, pm = L(f"paper_data/fig10a_measured_{q}.csv"), L(f"paper_data/fig10a_paper_model_{q}.csv")
        axx.plot(m[:, 0], m[:, 1], "o", mfc="none", c=c, ms=5, label=f"{q.split('_')[0]} measured [P]")
        axx.plot(*pm.T, "-", c=c, lw=1.2, alpha=0.6, label=f"{q.split('_')[0]} paper's model [P]")
        axx.plot(A[:, 0], A[:, k], "--", c=c, lw=2, label=f"{q.split('_')[0]} this work")
    ax[0].set(xlabel="P_DEL (dBm)", ylabel="P_OUT (dBm), Gain (dB)", xlim=(-7, 17), ylim=(0, 30))
    a2.set(ylabel="PAE (%)", ylim=(0, 50))
    ax[0].set_title("(a) 10 GHz, V_DS = 15 V, I_DS = 25 mA; not fitted", loc="left", fontsize=10)
    ax[0].legend(frameon=False, fontsize=7, loc="upper left"); a2.legend(frameon=False, fontsize=7, loc="lower right")
    b2 = ax[1].twinx()
    for k, q, c, axx in ((1, "Pfund_dBm", "#a94b17", ax[1]), (2, "IM3_dBm", "#2a78d6", ax[1]), (3, "OIP3_over_PDC_dB", "#1d5b3a", b2)):
        m, pm = L(f"paper_data/fig10b_measured_{q}.csv"), L(f"paper_data/fig10b_paper_model_{q}.csv")
        lab = q.split("_dB")[0]
        axx.plot(m[:, 0], m[:, 1], "o", mfc="none", c=c, ms=5, label=f"{lab} measured [P]")
        axx.plot(*pm.T, "-", c=c, lw=1.2, alpha=0.6, label=f"{lab} paper's model [P]")
        axx.plot(B[:, 0], B[:, k], "--", c=c, lw=2, label=f"{lab} this work")
    ax[1].set(xlabel="P_DEL per tone (dBm)", ylabel="P_f0, IM3 (dBm)", xlim=(-2, 7), ylim=(-70, 30))
    b2.set(ylabel="OIP3/P_DC (dB)", ylim=(-60, 20))
    ax[1].set_title("(b) two tones 10 / 10.01 GHz; not fitted", loc="left", fontsize=10)
    ax[1].legend(frameon=False, fontsize=7, loc="lower right"); b2.legend(frameon=False, fontsize=7, loc="center right")
    for a in ax:
        a.grid(c="#e1e0d9")
    plt.tight_layout()
    plt.savefig("fig10_device1.png", dpi=150)
    print("saved fig10_device1.png")
