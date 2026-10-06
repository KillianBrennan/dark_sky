"""Run the full pipeline:  .venv/bin/python -m darksky [step ...]

Steps (default: all, in order):
  fetch      isochrones (+ KML), Copernicus DEM, WorldCover, Lorenz tiles
  viirs      clip VIIRS (requires config.VIIRS_FILE, see viirs.py)
  terrain    horizons on the 100 m and 250 m grids
  calibrate  fit scale + AOD against the Lorenz atlas
  model      light-dome model (central + sensitivity runs)
  layers     fields GeoTIFF, RGBA COG layers, legends, sensitivity summary
"""
import sys

STEPS = ["fetch", "viirs", "terrain", "calibrate", "model", "layers"]


def main(steps):
    for s in steps:
        print(f"== {s}", flush=True)
        if s == "fetch":
            from . import fetch, isochrone
            isochrone.fetch()
            isochrone.to_kml()
            fetch.fetch_dem(), fetch.fetch_worldcover()
            fetch.fetch_lorenz(2022), fetch.fetch_lorenz(2025)
        elif s == "viirs":
            from . import viirs
            viirs.clip()
        elif s == "terrain":
            from . import run_terrain
            run_terrain.run()
        elif s == "calibrate":
            from . import calibrate
            calibrate.run()
        elif s == "model":
            from . import run_model
            run_model.run()
            run_model.run("F030", F=0.30)
            run_model.run("leak005", leak=0.05, bands=("NIR",))
            run_model.run("leak035", leak=0.35, bands=("NIR",))
        elif s == "layers":
            from . import build_layers
            build_layers.run()


if __name__ == "__main__":
    main(sys.argv[1:] or STEPS)
