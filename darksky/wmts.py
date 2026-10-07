"""Static WMTS (capabilities XML + PNG tiles) for the all-sky layers, so map.geo.admin.ch
shows a title, description and colour legend under the layer's info button (COG layers
cannot carry a legend in the viewer).

  .venv/bin/python -m darksky.wmts

Layout (served from GitHub raw):
  outputs/wmts/WMTSCapabilities.xml
  outputs/wmts/<layer>/<TileMatrix>/<TileRow>/<TileCol>.png
Tile matrix set "lv95_darksky": EPSG:2056, origin at the domain's top-left corner,
256 px tiles at 800/400/200/100/50/25 m per pixel. Levels finer than the 100 m data are
pixel-replicated (crisp when zoomed in); coarser levels take the darkest opaque sub-pixel
for the masked layer (isolated good sites stay visible) and the central sub-pixel otherwise.
Map link: layers=WMTS|<capabilities url>|<layer identifier>
"""
from xml.sax.saxutils import escape

import numpy as np
import rasterio
from PIL import Image

from . import config as C
from .fetch import domain_bounds, domain_lv95

RAW = "https://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/"
RES = [800.0, 400.0, 200.0, 100.0, 50.0, 25.0]
TILE = 256
TMS = "lv95_darksky"

LAYERS = [
    dict(id="darksky_allsky30", src="layer_3_allsky30.tif", legend="legend_3_allsky30.png", reduce="darkest",
         title="Dark sky: sky above 30°, candidate sites",
         abstract=("Modelled mean V-band sky brightness over the sky above 30° elevation (mag/arcsec²), including "
                   "local light sources and direct glare from lamps in view. Shown only where terrain stays below 10° "
                   "over at least 75 % of the horizon, at least 200 m from any building, less than 50 % forest, and "
                   "artificial light at most 2× the natural sky. Darker blue = darker sky. Model output (VIIRS 2025, "
                   "calibrated to the Lorenz 2025 atlas), not a measurement. github.com/KillianBrennan/dark_sky")),
    dict(id="darksky_allsky30_unmasked", src="layer_3b_allsky30_unmasked.tif", legend="legend_3b_allsky30_unmasked.png",
         reduce="centre", title="Dark sky: sky above 30°, everywhere",
         abstract=("Modelled mean V-band sky brightness over the sky above 30° elevation (mag/arcsec²), including "
                   "local light sources and direct glare, for every pixel (no site masks). The lightest class is "
                   "brighter than the candidate-site threshold (artificial > 2× natural). Darker blue = darker sky. "
                   "Model output, not a measurement. github.com/KillianBrennan/dark_sky")),
]


def _level(rgba, factor, mode):
    """Resample the (H, W, 4) 100 m image by `factor` (>1 coarser, <1 finer)."""
    if factor < 1:
        k = int(round(1 / factor))
        return np.repeat(np.repeat(rgba, k, axis=0), k, axis=1)
    if factor == 1:
        return rgba
    k = int(factor)
    h, w = rgba.shape[0] // k * k, rgba.shape[1] // k * k
    blk = rgba[:h, :w].reshape(h // k, k, w // k, k, 4).transpose(0, 2, 1, 3, 4).reshape(h // k, w // k, k * k, 4)
    if mode == "darkest":
        lum = blk[..., :3].astype(np.int32).sum(-1)
        lum = np.where(blk[..., 3] > 0, lum, 10_000)
        i = lum.argmin(-1)
    else:
        i = np.full(blk.shape[:2], (k // 2) * k + k // 2)
    return np.take_along_axis(blk, i[..., None, None], axis=2)[:, :, 0, :]


def build():
    out = C.OUT / "wmts"
    minx, miny, maxx, maxy = domain_lv95()
    matrices = []
    for z, res in enumerate(RES):
        n_w = int(np.ceil((maxx - minx) / res / TILE))
        n_h = int(np.ceil((maxy - miny) / res / TILE))
        matrices.append((z, res, n_w, n_h))
    for L in LAYERS:
        with rasterio.open(C.OUT / L["src"]) as s:
            assert abs(s.transform.a - 100.0) < 1e-6 and abs(s.transform.c - minx) < 1e-6
            rgba = s.read().transpose(1, 2, 0)
        n_files = 0
        for z, res, n_w, n_h in matrices:
            img = _level(rgba, res / 100.0, L["reduce"])
            canvas = np.zeros((n_h * TILE, n_w * TILE, 4), np.uint8)
            canvas[:img.shape[0], :img.shape[1]] = img
            for r in range(n_h):
                for c in range(n_w):
                    d = out / L["id"] / str(z) / str(r)
                    d.mkdir(parents=True, exist_ok=True)
                    Image.fromarray(canvas[r * TILE:(r + 1) * TILE, c * TILE:(c + 1) * TILE], "RGBA").save(
                        d / f"{c}.png", optimize=True)
                    n_files += 1
        print(f"{L['id']}: {n_files} tiles", flush=True)
    lon0, lat0, lon1, lat1 = domain_bounds()
    tm = "".join(
        f"""
      <TileMatrix><ows:Identifier>{z}</ows:Identifier><ScaleDenominator>{res / 0.00028:.6f}</ScaleDenominator>
        <TopLeftCorner>{minx:.1f} {maxy:.1f}</TopLeftCorner><TileWidth>{TILE}</TileWidth><TileHeight>{TILE}</TileHeight>
        <MatrixWidth>{n_w}</MatrixWidth><MatrixHeight>{n_h}</MatrixHeight></TileMatrix>""" for z, res, n_w, n_h in matrices)
    layers = "".join(
        f"""
    <Layer>
      <ows:Title>{escape(L['title'])}</ows:Title>
      <ows:Abstract>{escape(L['abstract'])}</ows:Abstract>
      <ows:WGS84BoundingBox><ows:LowerCorner>{lon0:.5f} {lat0:.5f}</ows:LowerCorner><ows:UpperCorner>{lon1:.5f} {lat1:.5f}</ows:UpperCorner></ows:WGS84BoundingBox>
      <ows:Identifier>{L['id']}</ows:Identifier>
      <Style isDefault="true"><ows:Identifier>default</ows:Identifier>
        <LegendURL format="image/png" xlink:href="{RAW}{L['legend']}"/></Style>
      <Format>image/png</Format>
      <TileMatrixSetLink><TileMatrixSet>{TMS}</TileMatrixSet></TileMatrixSetLink>
      <ResourceURL format="image/png" resourceType="tile" template="{RAW}wmts/{L['id']}/{{TileMatrix}}/{{TileRow}}/{{TileCol}}.png"/>
    </Layer>""" for L in LAYERS)
    xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<Capabilities xmlns="http://www.opengis.net/wmts/1.0" xmlns:ows="http://www.opengis.net/ows/1.1"
  xmlns:xlink="http://www.w3.org/1999/xlink" version="1.0.0">
  <ows:ServiceIdentification>
    <ows:Title>Dark sky layers around Bern</ows:Title>
    <ows:ServiceType>OGC WMTS</ows:ServiceType><ows:ServiceTypeVersion>1.0.0</ows:ServiceTypeVersion>
  </ows:ServiceIdentification>
  <ows:ServiceProvider><ows:ProviderName>dark_sky (GitHub)</ows:ProviderName>
    <ows:ProviderSite xlink:href="https://github.com/KillianBrennan/dark_sky"/></ows:ServiceProvider>
  <Contents>{layers}
    <TileMatrixSet>
      <ows:Identifier>{TMS}</ows:Identifier>
      <ows:SupportedCRS>urn:ogc:def:crs:EPSG::2056</ows:SupportedCRS>{tm}
    </TileMatrixSet>
  </Contents>
</Capabilities>
"""
    (out / "WMTSCapabilities.xml").write_text(xml)
    return out / "WMTSCapabilities.xml"


if __name__ == "__main__":
    print(build())
