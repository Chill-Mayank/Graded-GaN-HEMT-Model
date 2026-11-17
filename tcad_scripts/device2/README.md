# Device 2 TCAD decks (Sentaurus T-2022.03-SP2, terminal flow)

Short-channel graded HEMT of the paper's Fig. 5(b). One change per run; each deck's header says what changed and why.

| File | What it does | Key result |
|---|---|---|
| `sde_dvs.cmd` | structure, mole fraction, regrown n+ contacts, mesh (Sentaurus Structure Editor) | |
| `n4_full.cmd` + `sdevice_n4.par` | equilibrium; relaxed Al0.04GaN buffer, layers strained to GaN (default a0) | n_s under the gate 3.18e12 cm^-2 |
| `n5_full.cmd` + `sdevice_n5.par` | as n4, layers strained to the buffer (a0 = 3.1859 A, Vegard) | n_s 2.35e12 cm^-2; junction charge -1.82e12 cm^-2 |
| `n6_full.cmd` | n5 + gate sweep -3 to +1 V at V_DS = 0.05 V | V_th (ELR) -0.84 V; 1 mA/mm at -0.86 V |
| `n7_full.cmd` | n6 with V_DS = 5 V (the bias of Fig. 6(b)) | 1 mA/mm at -1.464 V (measured -1.457 V) |
| `n8_full.cmd` | n6 + TDR snapshots every 0.1 V from -1.4 to +0.4 V | n_s(V_GS) under the gate -> C_g'', V_OFF |
| `n8_ns.py` | reads the n8 snapshots (h5py): n_s under the gate and vertical profiles | used for docs/data/dev2_n8_ns.csv |

Run from the folder that holds the decks (`pwd` first), in the background so a dropped SSH session does not stop it:

```bash
nohup sdevice n7_full.cmd > n7.log 2>&1 &
```

A 5 V gate sweep takes about 3 h of wall time (12.8 h CPU, 4.5 GB) on our server; a 0.05 V sweep about 1 h.
