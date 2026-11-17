"""plot_fig6b.py - Step 8b figure: Device 2 transfer and G_M at V_DS = 5 V, measured vs the paper's model vs this work.
Run: GAN_DEVICE=2 python3 plot_fig6b.py   (needs params_device2.json from fit2.py) -> fig6b_device2.png"""
import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
import fit2 as f2
p = list(json.load(open("params_device2.json"))["fitted [CALIBRATE]"].values())
V = np.linspace(-1.6, 1.0, 105)
I = [f2.model_id(p, v, 5.0) for v in V]; G = [f2.model_gm(p, v) for v in V]
mi, mg = f2.load("fig6b_measured_ID"), f2.load("fig6b_measured_gm")
li, lg = f2.load("fig6b_paper_model_ID"), f2.load("fig6b_paper_model_gm")
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].plot(mi[:, 0], mi[:, 1], "o", mfc="none", c="#52514e", label="measured [P]")
ax[0].plot(li[:, 0], li[:, 1], "-", c="#2a78d6", label="paper's model [P]")
ax[0].plot(V, I, "--", c="#eb6834", label="this work (A)")
ax[0].set(xlabel="V_GS (V)", ylabel="I_D (mA)", xlim=(-3, 1), title="Device 2, V_DS = 5 V, W = 75 um")
ax[1].plot(mg[:, 0], mg[:, 1] * 0.075, "o", mfc="none", c="#52514e")
ax[1].plot(lg[:, 0], lg[:, 1] * 0.075, "-", c="#2a78d6")
ax[1].plot(V, G, "--", c="#eb6834")
ax[1].set(xlabel="V_GS (V)", ylabel="G_M (mS)", xlim=(-3, 1), title="G_M (printed axis x 0.075 mm)")
for a in ax: a.grid(c="#e1e0d9"); a.spines[["top", "right"]].set_visible(False)
ax[0].legend(frameon=False)
plt.tight_layout(); plt.savefig("fig6b_device2.png", dpi=150)
