"""Static WMTS (capabilities XML + PNG tiles) for the all-sky layers, so map.geo.admin.ch
shows a title, description and colour legend under the layer's info button (COG layers
cannot carry a legend in the viewer).

  .venv/bin/python -m darksky.wmts

Layout (served from GitHub raw):
  outputs/wmts/WMTSCapabilities.xml
  outputs/wmts/<layer>/<TileMatrix>/<TileRow>/<TileCol>.png
Tile matrix set "lv95_darksky": EPSG:2056, origin at the domain's top-left corner,
256 px tiles at the viewer's own LV95 zoom resolutions 650/500/250/100/50/20/10 m, so no
resampling (= no blur) happens at those zooms. Finer levels replicate the 100 m pixels
(nearest neighbour); coarser levels take the darkest opaque source pixel in each output
pixel's footprint for the masked layer (isolated good sites stay visible) and the nearest
source pixel otherwise.
Map link: layers=WMTS|<capabilities url>|<layer identifier>
"""
from xml.sax.saxutils import escape

import numpy as np
import rasterio
from PIL import Image

from . import config as C
from .fetch import domain_bounds, domain_lv95

RAW = "https://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/"
# exactly the viewer's LV95 zoom resolutions (zoom 0..6), so tiles are drawn 1:1 without
# resampling; finer zooms upscale the 10 m tiles, in which a 100 m pixel is a 10x10 block
RES = [650.0, 500.0, 250.0, 100.0, 50.0, 20.0, 10.0]
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


def _render(rgba, res, n_w, n_h, mode, src_res=100.0):
    """Resample the (H, W, 4) 100 m image onto a canvas of n_h x n_w tiles at `res`."""
    H, W = rgba.shape[:2]
    oh, ow = n_h * TILE, n_w * TILE
    out = np.zeros((oh, ow, 4), np.uint8)
    if mode == "darkest" and res > src_res:
        # footprint of each output pixel in source pixels: [edge_k, edge_k+1)
        er = np.minimum(np.floor(np.arange(oh + 1) * res / src_res).astype(int), H)
        ec = np.minimum(np.floor(np.arange(ow + 1) * res / src_res).astype(int), W)
        rows = np.nonzero(er[:-1] < H)[0]
        cols = np.nonzero(ec[:-1] < W)[0]
        key = np.where(rgba[..., 3] > 0, 765 - rgba[..., :3].astype(np.int32).sum(-1), -1)
        k1 = np.maximum.reduceat(key, er[rows], axis=0)
        k2 = np.maximum.reduceat(k1, ec[cols], axis=1)
        # map darkness keys back to colours (few distinct colours)
        pal = {}
        flat_k, flat_c = key.ravel(), rgba.reshape(-1, 4)
        for k in np.unique(flat_k):
            pal[k] = flat_c[np.argmax(flat_k == k)] if k >= 0 else np.zeros(4, np.uint8)
        block = np.zeros(k2.shape + (4,), np.uint8)
        for k, c in pal.items():
            block[k2 == k] = c
        out[rows[0]:rows[-1] + 1, cols[0]:cols[-1] + 1] = block
        return out
    ri = np.floor((np.arange(oh) + 0.5) * res / src_res).astype(int)
    ci = np.floor((np.arange(ow) + 0.5) * res / src_res).astype(int)
    vr, vc = ri < H, ci < W
    out[np.ix_(vr, vc)] = rgba[np.ix_(ri[vr], ci[vc])]
    return out


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
            canvas = _render(rgba, res, n_w, n_h, L["reduce"])
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
