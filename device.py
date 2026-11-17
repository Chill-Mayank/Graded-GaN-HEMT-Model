"""device.py - Step 8a: picks the device the model runs for. GAN_DEVICE=2 python3 fit2.py -> Device 2;
anything else (or unset) -> Device 1, so every earlier script runs exactly as before."""
import os
if os.environ.get("GAN_DEVICE", "1") == "2":
    import device2 as dev
else:
    import device1 as dev
