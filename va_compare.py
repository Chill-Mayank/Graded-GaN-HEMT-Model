"""
va_compare.py - Step 10a: compare Spectre's DC results for va/graded_hemt.va with the Python model and with va_twin.py.
Run:  python3 va_compare.py        (reads va/dc_device1.raw and va/dc_device2.raw, Spectre nutascii output)
Pass criteria [ASSUME, set before the run]:
  against va_twin.py (same equations, Eq. (2) F_1/2)   : I_D within 0.01 %  -> the Verilog-A is the twin, coded right
  against the Python model (exact F_1/2)               : I_D within 0.5 %, G_M within 1 % of its peak
"""
import os
import re
import subprocess
import sys
import numpy as np


def read_nutascii(fn):
    """All plots of a nutmeg ASCII raw file as {plotname: (names, values[npoints, nvars])}."""
    txt = open(fn, errors="replace").read()
    plots = {}
    for block in re.split(r"(?m)^Plotname:", txt)[1:]:
        name = block.splitlines()[0].strip()
        nv = int(re.search(r"No\. Variables:\s*(\d+)", block).group(1))
        head = re.split(r"(?m)^Variables:", block.split("Values:", 1)[0], maxsplit=1)[1]   # variable 0 may share the line
        names = [ln.split()[1] for ln in head.strip().splitlines() if ln.split()][:nv]
        vals = block.split("Values:", 1)[1]
        vals = re.split(r"(?m)^(?:Title|Date):", vals)[0]
        tok = np.array([float(t) for t in vals.replace(",", " ").split()])
        tok = tok[: (len(tok) // (nv + 1)) * (nv + 1)].reshape(-1, nv + 1)[:, 1:]
        plots[name] = (names, tok)
    return plots


def pick(plots, analysis):
    """(x, I_D) of one analysis: x is the swept variable, I_D = -(current into VD's + terminal)."""
    for name, (names, v) in plots.items():
        if re.search(rf"\b{analysis}\b", name):
            k = [i for i, s in enumerate(names) if s.lower().endswith(":p")][0]
            return v[:, 0], -v[:, k]
    raise KeyError(f"analysis {analysis} not found; plots are: {list(plots)}")


def check(dev):
    import va_twin
    from va_netlists import SWEEPS
    if dev == 1:
        import fit as F
        import json
        p = list(json.load(open("params_device1.json"))["fitted [CALIBRATE]"].values())
    else:
        import fit2 as F
        import json
        p = list(json.load(open("params_device2.json"))["fitted [CALIBRATE]"].values())
    tw = va_twin.Twin(va_twin.card(dev))
    s = SWEEPS[dev]
    plots = read_nutascii(f"va/dc_device{dev}.raw")

    print(f"===== Device {dev}: Spectre (graded_hemt.va) against va_twin.py and the Python model")
    vg, i_sp = pick(plots, "transfer")
    if np.median(i_sp[i_sp != 0]) < 0:
        i_sp = -i_sp
        print("  (current sign flipped: this Spectre reports VD:p as the current out of the + terminal)")
    i_py = np.array([F.model_id(p, v, s["vds"]) for v in vg]) * 1e-3
    i_tw = np.array([tw.id_ext(v, s["vds"]) for v in vg])
    on = i_py > 0.05 * i_py.max()
    gm_sp, gm_py = np.gradient(i_sp, vg), np.gradient(i_py, vg)
    print(f"  transfer, V_DS = {s['vds']} V, {len(vg)} points, V_GS {vg[0]:+.2f} .. {vg[-1]:+.2f} V")
    print(f"    I_D vs twin  : worst {100 * np.max(np.abs(i_sp[on] / i_tw[on] - 1)):.4f} %   (criterion 0.01 %)")
    print(f"    I_D vs Python: worst {100 * np.max(np.abs(i_sp[on] / i_py[on] - 1)):.3f} %    (criterion 0.5 %)")
    print(f"    G_M vs Python: worst {100 * np.max(np.abs(gm_sp - gm_py)) / gm_py.max():.3f} % of the peak"
          f" {gm_py.max() * 1e3:.2f} mS   (criterion 1 %)")
    worst_out = []
    for i, vgo in enumerate(s["vgs_out"]):
        vd, io = pick(plots, f"out_k{i}")
        io = io if np.median(io[io != 0]) > 0 else -io
        ip = np.array([F.model_id(p, vgo, v) for v in vd]) * 1e-3
        m = ip > 0.05 * max(ip.max(), 1e-12)
        worst_out.append(100 * np.max(np.abs(io[m] / ip[m] - 1)) if m.any() else 0.0)
    print(f"  output curves, V_GS = {', '.join(f'{v:+g}' for v in s['vgs_out'])} V: I_D vs Python worst "
          f"{max(worst_out):.3f} %   (criterion 0.5 %)")
    for v in (s["vg"][0], -1.0, 0.0, s["vg"][1]):
        k = int(np.argmin(abs(vg - v)))
        print(f"    V_GS {vg[k]:+5.2f} V: Spectre {i_sp[k] * 1e3:10.5f} mA   twin {i_tw[k] * 1e3:10.5f} mA"
              f"   Python {i_py[k] * 1e3:10.5f} mA")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        check(int(sys.argv[1]))
    else:
        for dev in ("1", "2"):
            if os.path.exists(f"va/dc_device{dev}.raw"):
                subprocess.run([sys.executable, __file__, dev], env={**os.environ, "GAN_DEVICE": dev}, check=True)
            else:
                print(f"va/dc_device{dev}.raw not found")
