"""Step 3 driver: horizon-derived fields on the 100 m grid, and horizons + shadow
angles on the 250 m light-dome grid."""
import time

import numpy as np

from . import config as C
from .grid import make_grid
from .terrain import BETA_D_KM, horizons, open_fraction

FINE_RES = 100
DOME_RES = 250


def run():
    for res, with_beta in ((DOME_RES, True), (FINE_RES, False)):
        t = time.time()
        g = make_grid(res)
        hor, hdist, beta = horizons(g["x"], g["y"], g["z"], with_beta=with_beta)
        # 5-degree bins: max horizon within each bin (centred on 0, 5, 10 ...)
        h5 = np.max(np.roll(hor, 2, axis=1).reshape(len(hor), 72, 5), axis=2)
        out = dict(rows=g["rows"], cols=g["cols"], x=g["x"], y=g["y"], z=g["z"], drive_min=g["drive_min"],
                   shape=np.array(g["shape"]), transform=np.array(g["transform"])[:6],
                   horizon5=h5.astype(np.float32), open15=open_fraction(hor).astype(np.float32))
        if with_beta:
            out |= dict(beta=beta.astype(np.float16), beta_d_km=BETA_D_KM)
        np.savez(C.INTERIM / f"terrain_{res}m.npz", **out)
        print(f"{res} m grid: {len(g['x'])} px, {time.time() - t:.0f}s", flush=True)


if __name__ == "__main__":
    run()
