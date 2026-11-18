"""
va_netlists.py - Step 10a: write the Spectre DC netlists that check va/graded_hemt.va against the Python model.
Run:  python3 va_netlists.py      -> va/dc_device1.scs, va/dc_device2.scs, va/sp_device1.scs
Every model parameter comes from device1.py / device2.py and params_device1.json / params_device2.json (va_twin.card),
so Spectre runs exactly the fitted model. Then, on the Spectre machine, inside va/:
  spectre dc_device1.scs -format nutascii -raw dc_device1.raw > dc_device1.log 2>&1
  spectre dc_device2.scs -format nutascii -raw dc_device2.raw > dc_device2.log 2>&1
and back here:  python3 va_compare.py
"""
import datetime
import os
from va_twin import card

# sweeps: transfer at the paper's bias [P, Figs. 6(a), 6(b)]; output curves up to the large-signal bias (Fig. 10: 15 V)
SWEEPS = {1: dict(vds=10.0, vg=(-3.0, 3.0), vgs_out=(-1.0, 0.0, 1.0, 2.0, 3.0), vd=(0.0, 15.0)),
          2: dict(vds=5.0, vg=(-3.0, 1.0), vgs_out=(-1.0, -0.5, 0.0, 0.5, 1.0), vd=(0.0, 7.0))}
NAMES = ["w", "l", "epsr", "db", "tg", "meff", "tdev", "voff", "gam", "mulf", "mua", "mub", "thsat", "vsat0", "th1",
         "th2", "nmax", "delta", "lsg", "lgd", "nacc", "u0", "vsatacc", "beta"]
KEYS = ["W", "L", "EPSR", "DB", "TG", "MEFF", "TEMP", "VOFF", "GAMMA", "MULF", "MUA", "MUB", "THSAT", "VSAT0", "TH1",
        "TH2", "NMAX", "DELTA", "LSG", "LGD", "NACC", "U0", "VSATACC", "BETA"]


def instance_params(c):
    pairs = [f"{n}={c[k]:.15g}" for n, k in zip(NAMES, KEYS)]
    lines, line = [], "    "
    for p in pairs:
        if len(line) + len(p) > 100:
            lines.append(line + "\\")
            line = "    "
        line += p + " "
    lines.append(line.rstrip())
    return "\n".join(lines)


def netlist(dev):
    c, s = card(dev), SWEEPS[dev]
    out = [f"// dc_device{dev}.scs - Step 10a: DC check of graded_hemt.va against the Python model, Device {dev}",
           f"// written by va_netlists.py on {datetime.date.today()} from params_device{dev}.json; do not edit by hand",
           "simulator lang=spectre",
           "global 0",
           'ahdl_include "graded_hemt.va"',
           f"parameters vgs=0 vds={s['vds']}",
           "VG (g 0) vsource dc=vgs",
           "VD (d 0) vsource dc=vds",
           "DUT (d g 0) graded_hemt \\",
           instance_params(c),
           "save VD:p",
           f"transfer dc param=vgs start={s['vg'][0]} stop={s['vg'][1]} step=0.01"]
    for i, vg in enumerate(s["vgs_out"]):
        out += [f"set_k{i} alter param=vgs value={vg}",
                f"out_k{i} dc param=vds start={s['vd'][0]} stop={s['vd'][1]} step=0.05"]
    return "\n".join(out) + "\n"


# Step 10c: S-parameters of Device 1 at the Fig. 11(a) bias, extrinsic values from fig11.py [CALIBRATE]
SP1 = dict(vgs=0.0, vds=10.0, rg=17.5, cgsfr=0.0, cgdfr=16.7e-15, cdsfr=0.0)


def sp_netlist(dev=1, s=SP1):
    c = card(dev)
    return "\n".join([
        f"// sp_device{dev}.scs - Step 10c: S-parameters of graded_hemt.va, Device {dev}, V_GS = {s['vgs']} V, V_DS = {s['vds']} V",
        f"// written by va_netlists.py on {datetime.date.today()}; ideal bias tees (1 H, 1 F) keep the 50-ohm ports off the DC path",
        "simulator lang=spectre",
        "global 0",
        'ahdl_include "graded_hemt.va"',
        f"VG (gb 0) vsource dc={s['vgs']}",
        "LG (gb g) inductor l=1",
        "CG (p1 g) capacitor c=1",
        "P1 (p1 0) port r=50 num=1",
        f"VD (db 0) vsource dc={s['vds']}",
        "LD (db d) inductor l=1",
        "CD (p2 d) capacitor c=1",
        "P2 (p2 0) port r=50 num=2",
        "DUT (d g 0) graded_hemt \\",
        instance_params(c) + " \\",
        f"    rg={s['rg']} cgsfr={s['cgsfr']:.6g} cgdfr={s['cgdfr']:.6g} cdsfr={s['cdsfr']:.6g}",
        f'sp1 sp ports=[P1 P2] start=1G stop=60G dec=20 file="sp_device{dev}.s2p" datafmt=touchstone',
    ]) + "\n"


# Step 10d: Fig. 10, Device 1 large signal at 10 GHz [P, Sec. III]: V_DS = 15 V, I_DS = 25 mA (V_GS = -0.33532 V from
# va_twin, [CALC]); Fig. 10(a) Z_S = 25 + 26j, Z_L = 45 + 96j ohm; Fig. 10(b) two tones 10 GHz and 10.01 GHz,
# Z_S = 20 + 26j, Z_L = 50 + 140j ohm. Each impedance is a resistor (the port) in series with an inductor giving the
# printed reactance at 10 GHz; at the harmonics it is R + j n X [ASSUME: the paper does not give harmonic terminations].
HB = dict(vgs=-0.33532, vds=15.0, f0=10e9, df=10e6)


def hb_netlist(two_tone=False, s=SP1, h=HB):
    c = card(1)
    rs, xs, rl, xl = (20, 26, 50, 140) if two_tone else (25, 26, 45, 96)
    w = 2 * 3.141592653589793 * h["f0"]
    src = (f'P1 (p1 0) port r={rs} num=1 type=sine freq={h["f0"]:.6g} dbm=prf fundname="f1" '
           f'freq2={h["f0"] + h["df"]:.6g} dbm2=prf fundname2="f2"' if two_tone else
           f'P1 (p1 0) port r={rs} num=1 type=sine freq={h["f0"]:.6g} dbm=prf fundname="f1"')
    # hb takes no start/stop/step (Spectre 20.1 warning SFE-106): the power sweep wraps it [Spectre RF manual]
    ana = ('swp sweep param=prf start=-6 stop=10 step=1 {\n'
           '  hb1 hb funds=["f1" "f2"] maxharms=[5 5] errpreset=conservative annotate=status\n}'
           if two_tone else
           'swp sweep param=prf start=-10 stop=22 step=1 {\n'
           '  hb1 hb funds=["f1"] maxharms=[9] errpreset=conservative annotate=status\n}')
    name = "hb2_device1" if two_tone else "hb_device1"
    return "\n".join([
        f"// {name}.scs - Step 10d: Fig. {'10(b) two-tone' if two_tone else '10(a) power sweep'}, Device 1, harmonic balance",
        f"// written by va_netlists.py on {datetime.date.today()}; Z_S = {rs}+{xs}j, Z_L = {rl}+{xl}j ohm at 10 GHz [P]",
        "simulator lang=spectre",
        "global 0",
        'ahdl_include "graded_hemt.va"',
        "parameters prf=0",
        f"VG (gb 0) vsource dc={h['vgs']}",
        "LG (gb g) inductor l=1",
        "CG (gs g) capacitor c=1",
        f"LS (p1 gs) inductor l={xs / w:.6g}",
        src,
        f"VD (db 0) vsource dc={h['vds']}",
        "LD (db d) inductor l=1",
        "CD (dl d) capacitor c=1",
        f"LL (p2 dl) inductor l={xl / w:.6g}",
        f"P2 (p2 0) port r={rl} num=2",
        "DUT (d g 0) graded_hemt \\",
        instance_params(c) + " \\",
        f"    rg={s['rg']} cgsfr={s['cgsfr']:.6g} cgdfr={s['cgdfr']:.6g} cdsfr={s['cdsfr']:.6g}",
        "save p1 p2 gs g d P1:p VD:p",
        ana,
    ]) + "\n"


if __name__ == "__main__":
    open("va/sp_device1.scs", "w").write(sp_netlist(1))
    print("wrote va/sp_device1.scs")
    open("va/hb_device1.scs", "w").write(hb_netlist(False))
    open("va/hb2_device1.scs", "w").write(hb_netlist(True))
    print("wrote va/hb_device1.scs, va/hb2_device1.scs")
    os.makedirs("va", exist_ok=True)
    for dev in (1, 2):
        fn = f"va/dc_device{dev}.scs"
        open(fn, "w").write(netlist(dev))
        print(f"wrote {fn}")
