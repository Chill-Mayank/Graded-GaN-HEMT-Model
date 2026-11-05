"""
fit.py - compact model step 5b: fit the Device 1 parameters to the measured points of the paper's Fig. 6(a) and Fig. 7(a).
The paper calls these "fitting parameters" and gives no values for Device 1 [P, below Eq. (17)].
Only the parameters the data can pin down are fitted (decided 2026-11-05: option A, then delta added).
Saves the result to params_device1.json, so later steps use the same numbers.
"""
import json
import time
import numpy as np
from scipy.optimize import least_squares
import current as cur
from access import id_ext, r_access, U0_ACC, N_ACC, L_SG, L_GD


def load(name):
    """Read one paper_data CSV (2 comment lines + 1 header line) into a numpy array."""
    return np.loadtxt(f"paper_data/{name}.csv", delimiter=",", skiprows=3)


# ---- the measured points we fit to [P, Fig. 6(a) and Fig. 7(a), measured = CW data of [11]]
id6 = load("fig6a_measured_ID")
id6 = id6[id6[:, 0] > -1.4]                                  # points on the zero line carry no information
gm6 = load("fig6a_measured_gm")
gm6 = gm6[(gm6[:, 0] > -1.4) & (gm6[:, 0] < 2.95)]           # the 3.0 V point disagrees with the I_D points [CALC]
out7 = load("fig7a_measured")
out7 = out7[out7[:, 1] > 0]                                  # V_DS = 0 gives I_D = 0 for any parameters

# ---- fitted: name in ids()/id_ext(), start, lower and upper bound, typical size, unit, meaning
FREE = [("u0",    U0_ACC,         100e-4, 3000e-4, 0.05, "m^2/Vs", "U0, access mobility [T Eq. 5.2]"),
        ("vsat0", cur.VSAT0,      2e4,    1e6,     1e5,  "m/s",    "V_SAT0 [P]"),
        ("th2",   cur.THETA_SAT2, 0.01,   2.0,     0.3,  "-",      "theta_SAT2 [P]"),
        ("delta", cur.DELTA_DS,   1.0,    10.0,    1.0,  "-",      "delta in V_DS,eff [P, no value]")]
NAMES = [f[0] for f in FREE]

# ---- fixed: the data cannot pin these down (seven-parameter test, Step 5b, 2026-11-04)
FIXED = {"vsat_acc": 1.8e5,     # m/s, access saturation velocity [TCAD, GaN.par saturation velocity]
         "beta": 2.0,           # access-resistance exponent [ASSUME]
         "mu_a": 0.0,           # m/V, the seven-parameter fit sent it to 0
         "theta_sat": 0.0}      # 1/V, the seven-parameter fit sent it to 0
W_ID, W_GM = 1.0, 2.0           # residual weights: per mA of I_D and per mS of g_M [ASSUME: g_M counts double]


def model_id(p, V_GS, V_DS):
    """Model drain current in mA for the fitted values p (in the order of NAMES)."""
    return id_ext(V_GS, V_DS, **FIXED, **dict(zip(NAMES, p))) * 1e3


def model_gm(p, V_GS, h=0.01):
    """Model g_M in mS at V_DS = 10 V, by a central difference of +-10 mV."""
    return (model_id(p, V_GS + h, 10.0) - model_id(p, V_GS - h, 10.0)) / (2 * h)


def residuals(p):
    """Model minus measurement for every point, weighted; least_squares makes the sum of their squares small."""
    r6 = [model_id(p, v, 10.0) - i for v, i in id6]
    r7 = [model_id(p, v, vd) - i for v, vd, i in out7]
    rg = [model_gm(p, v) - g for v, g in gm6]
    return np.r_[W_ID * np.array(r6), W_ID * np.array(r7), W_GM * np.array(rg)]


def rms(p):
    """Root-mean-square error of each data set: Fig. 6(a) I_D (mA), Fig. 7(a) I_D (mA), Fig. 6(a) g_M (mS)."""
    r = residuals(p)
    n6, n7 = len(id6), len(out7)
    return (np.sqrt(np.mean((r[:n6] / W_ID) ** 2)),
            np.sqrt(np.mean((r[n6:n6 + n7] / W_ID) ** 2)),
            np.sqrt(np.mean((r[n6 + n7:] / W_GM) ** 2)))


def run_fit(p0):
    """Least-squares fit from the start values p0; returns scipy's result."""
    return least_squares(residuals, p0,
                         bounds=([f[2] for f in FREE], [f[3] for f in FREE]),
                         x_scale=[f[4] for f in FREE], max_nfev=300)


if __name__ == "__main__":
    p0 = [f[1] for f in FREE]
    print(f"Fitting to {len(id6)} I_D points (Fig. 6(a)), {len(out7)} output points (Fig. 7(a)), {len(gm6)} g_M points")
    print("start:  rms I_D = {:.2f} mA, output = {:.2f} mA, g_M = {:.2f} mS".format(*rms(p0)))

    t0 = time.time()
    sol = run_fit(p0)
    print("fitted: rms I_D = {:.2f} mA, output = {:.2f} mA, g_M = {:.2f} mS".format(*rms(sol.x)),
          f"  ({time.time() - t0:.0f} s)")

    sol2 = run_fit([0.03, 3e5, 0.1, 6.0])            # a very different start, to check the answer does not change

    print()
    print(f"{'parameter':36s}{'start':>11s}{'fitted':>12s}{'other start':>14s}  unit")
    for (name, start, lo, hi, scale, unit, meaning), v, v2 in zip(FREE, sol.x, sol2.x):
        note = "  <- at its bound" if min(v - lo, hi - v) < 1e-3 * scale else ""
        print(f"{meaning:36s}{start:11.4g}{v:12.4g}{v2:14.4g}  {unit}{note}")
    Rs0 = r_access(0.0, L_SG, u0=sol.x[0])                 # ohm, access resistance at zero current
    Rd0 = r_access(0.0, L_GD, u0=sol.x[0])
    print(f"U0 = {sol.x[0] * 1e4:.0f} cm^2/Vs, so R_S0 = {Rs0:.2f} ohm ({Rs0 * cur.W * 1e3:.2f} ohm.mm), "
          f"R_D0 = {Rd0:.2f} ohm ({Rd0 * cur.W * 1e3:.2f} ohm.mm)")

    result = {
        "source": "fit.py, Step 5b: least squares to the measured points of Nazir et al. 2025 Fig. 6(a) and Fig. 7(a)",
        "fitted [CALIBRATE]": dict(zip(NAMES, map(float, sol.x))),
        "fixed": {**FIXED, "v_off": cur.V_OFF_TCADP, "mu_lf": cur.MU_LF, "th1": cur.THETA_SAT1, "n_max": cur.N_MAX,
                  "mu_b": cur.MU_B, "n_acc_cm2": N_ACC, "l_sg": L_SG, "l_gd": L_GD},
        "rms": dict(zip(["fig6a_ID_mA", "fig7a_ID_mA", "fig6a_gM_mS"], map(float, rms(sol.x)))),
    }
    with open("params_device1.json", "w") as f:
        json.dump(result, f, indent=2)
    print()
    print("saved params_device1.json")