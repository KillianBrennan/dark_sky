"""Step 4 driver: azimuth-resolved artificial sky brightness on the 250 m grid, with
terrain blocking, for both bands. Saves kernel-unit radiance (scale applied later)."""
import json
import sys
import time

import numpy as np

from . import config as C
from . import lightdome as L
from .viirs import sources

SE_BINS = np.arange(45 // C.AZ_BIN, 225 // C.AZ_BIN + 1)        # bins centred 45..225 deg
H_IDX = {h: i for i, h in enumerate(C.ELEVATIONS)}


def want_masks():
    n_az = 360 // C.AZ_BIN
    wv = np.zeros((len(C.ELEVATIONS), n_az), bool)
    wv[H_IDX[90], 0] = True
    wv[H_IDX[25], SE_BINS] = True
    wv[H_IDX[45], np.arange(0, n_az, 45 // C.AZ_BIN)] = True
    wn = np.zeros_like(wv)
    for h in (20, 25, 30):
        wn[H_IDX[h], SE_BINS] = True
    return {"V": wv, "NIR": wn}


def run(tag="central", F=C.GARSTANG_F, leak=C.SHADOW_LEAK, bands=("V", "NIR")):
    cal = json.loads((C.OUT / "calibration.json").read_text())
    K = L.kernel(aod=cal["aod"], F=F)
    d = np.load(C.INTERIM / "terrain_250m.npz")
    near, far = sources()
    out = {}
    for ib, band in enumerate(C.BANDS):
        if band not in bands:
            continue
        t = time.time()
        B, Bb = L.site_brightness(K[ib], (d["x"], d["y"], d["z"]), d["beta"].astype(np.float32),
                                  d["beta_d_km"], near, far, want_masks()[band], use_terrain=True, leak=leak)
        out[f"B_{band}"] = B.astype(np.float32)
        out[f"Bblk_{band}"] = Bb.astype(np.float32)
        print(f"[{tag}] {band}: {time.time() - t:.0f}s", flush=True)
    np.savez_compressed(C.INTERIM / f"model_{tag}.npz", **out, aod=cal["aod"], F=F, leak=leak)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "central"
    if which == "central":
        run()
    elif which == "F030":
        run("F030", F=0.30, bands=("V", "NIR"))
    elif which == "leak005":
        run("leak005", leak=0.05, bands=("NIR",))
    elif which == "leak035":
        run("leak035", leak=0.35, bands=("NIR",))
