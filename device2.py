"""
device2.py - Device 2 (short-channel PolFET, paper Fig. 5(b)): numbers for the compact model.
Step 8a, 2026-11-17. Decision (user, 17 Nov): "everything same as paper" -> the paper's printed forms where the paper
prints them (C_g'' = eps/(2 d_AlGaN + t)); our TCAD where the paper gives no value (V_OFF, mobility, access density).
Sources: [P] Nazir et al., IEEE TED 72(9), 2025; [TCAD] our Device 2 deck, runs n5-n8.
"""
Q = 1.602176634e-19           # C
EPS0_SI = 8.8541878128e-12    # F/m

NAME = "Device 2"
# ---- structure [P, Sec. III and Fig. 5(b)]
W = 75e-6                     # m, total width: two fingers of 37.5 um [P, Sec. III]
L = 50e-9                     # m, gate length [P, Sec. III and Fig. 5(b)]
D_ALGAN = 9e-9                # m, Al0.27GaN barrier [P, Sec. II, the Fig. 1(a)/4(b) stack; Fig. 5(b)]
T_GRADE = 6e-9                # m, graded AlGaN, x = 0.15 -> 0 [P, Fig. 5(b)]
EPS_R = 9.5                   # [ASSUME, as Device 1 and the Fig. 4(b) fit]
L_SG = 0.6e-6                 # m, SiN between source and gate [P, Fig. 5(b)]
L_GD = 0.6e-6                 # m, SiN between gate and drain [P, Fig. 5(b)]

# ---- gate capacitance: the printed form [P, below Eq. (8)]; TCAD n8 supports eps/(d + t) (printed item 13, open)
CG = EPS_R * EPS0_SI / (2 * D_ALGAN + T_GRADE) * 1e-4    # F/cm^2
U = CG / Q                                               # cm^-2 per V
CG_DT = EPS_R * EPS0_SI / (D_ALGAN + T_GRADE) * 1e-4     # F/cm^2, eps/(d + t), for the sensitivity check only

# ---- V_OFF [TCAD]: the paper prints no value for Device 2, and its Eq. (8) gives +0.12 V (normally-off), which the
# measured Fig. 6(b) rules out. -0.63 V: fit of the paper's charge equation to the n8 sheet densities (V_DS = 0.05 V);
# -0.38 V: shift of the ELR threshold from V_DS = 0.05 V (n6) to 5 V (n7), because the model has no DIBL and the
# fitted data are at V_DS = 5 V. Option A, decided 2026-11-17.
V_OFF_N8 = -0.63              # V [TCAD n8]
DIBL_5V = -0.38               # V [TCAD n6 -> n7, ELR]
V_OFF = V_OFF_N8 + DIBL_5V    # V, -1.01

# ---- transport and access [TCAD]
MU_LF = 1380e-4               # m^2/Vs, n-weighted eMobility under the gate, n8 at V_GS = 0 (Masetti, undoped GaN)
N_ACC = 3.6e12                # cm^-2, access sheet density, n8 V_GS = 0 cuts at y = 0.80 / 1.45 um (3.63 / 3.58e12)
U0_ACC = 1380e-4              # m^2/Vs, access mobility start value [TCAD, same cuts]
