# Dark-sky layers around Bern

Shaded raster layers of modelled night-sky quality for a 140 × 140 km square centred on Bern
(LV95 2531000–2671000 E, 1130000–1270000 N, 100 m pixels), one layer per observing target.
They import directly into [map.geo.admin.ch](https://map.geo.admin.ch).

**These are model outputs, not measurements.** Read [Caveats](#caveats) before trusting a
specific hill.

## Importing into map.geo.admin.ch

1. Open map.geo.admin.ch → menu → **Tools / Werkzeuge** → **Import** → **Local file**.
2. Pick one of the `outputs/layer_*.tif` files (Cloud-Optimised GeoTIFF, RGBA, EPSG:2056).
   The import runs in your browser; the file is not uploaded anywhere.
3. Set the layer opacity in the layer list (≈60–70 % over the topographic basemap works well).
4. Optional: import `outputs/isochrones_bern.kml` to overlay the 15/30/45/60-min night-time
   drive-time contours from Bern.

Colour keys are in `outputs/legend_*.png`. All brightness layers use one blue ramp:
**darker blue = darker sky (better)**.

| File | Target | Shows |
|---|---|---|
| `layer_1_airglow_SE.tif` | 1. OH airglow gravity waves | Artificial NIR (Gen 3, ~850 nm) skyglow in the darkest 30°-wide window at 20–30° elevation within azimuth 45–225° (E through S), in magnitudes relative to the best 1 % of the square. Grey hatching marks pixels where terrain blocks every E–S window. |
| `layer_2_zenith.tif` | 2. Hα nebulae, 3. deep sky | Zenith sky brightness, V band, LED-corrected (mag/arcsec²). |
| `layer_3_allsky30.tif` | 3. Milky Way / binocular | Mean V sky brightness above 30° elevation (mag/arcsec²). |
| `layer_4_meteors.tif` | 4. Meteors | Relative visual meteor rate: open-sky fraction above 15° × r^(NELM − NELM_pristine), r = 2.5. |
| `layer_5_above_inversion.tif` | Autumn/winter bonus | Probability that the site lies above the Oct–Feb Mittelland stratus top (orange; transparent = below). |
| `layer_ref_lorenz2025_zenith.tif` | Reference | Lorenz 2025 atlas zenith brightness, an independent model with no terrain and no LED correction. |

`outputs/darksky_fields_lv95.tif` holds every field as float32 bands (named) for analysis in
QGIS/Python, including the S- and E-sector 25° V brightness, the best airglow azimuth, the
terrain-blocked fraction, open-sky fraction, elevation and drive-time bin.

## Results (central run)

**Calibration vs Lorenz 2025** (terrain off, all 313,600 pixels at 250 m): best AOD 0.12,
RMS 0.042 mag, 95th-percentile |residual| 0.09 mag, r = 0.995. Mean residual is within
±0.005 mag in every brightness class from 18 to 22 mag/arcsec² (`outputs/calibration.json`).
This is agreement with another VIIRS-driven model, **not** validation against the real sky.

**Terrain effect.** With terrain on, the zenith comes out a median 0.12 mag darker than the
terrain-free calibration (2–98 %: 0.03–0.27 mag). At 25° elevation in the S/E sector,
terrain removes a median 28 % of the artificial V light (2–98 %: 16–39 %), averaged over
all sources. For one city behind a ridge the per-source fraction is higher; for the zenith
it is far lower. This agrees in direction with the 50–70 % / "far less at zenith" guide in
the brief, but it is the outcome of the geometry, not an imposed split.

**Sensitivity** (`outputs/sensitivity.json`):

| Change | Airglow rank ρ | Top-5 % airglow area kept | Zenith |
|---|---|---|---|
| F = 0.30 (2× unseen horizontal emission) | 0.9997 | 97 % | 0.48 mag brighter, rank ρ 0.999 |
| shadow leak 5 % instead of 15 % | 0.9997 | 95 % | – |
| shadow leak 35 % | 0.9993 | 94 % | – |
| LED factor 1.2 / 2.0 instead of 1.5 | – | – | +0.14 / −0.19 mag, uniform (no rank change) |

The *ranking* is robust to the model's free parameters. The *absolute* magnitudes are not
(± ≈0.3 mag from the LED and horizontal-emission priors alone), and none of these runs
probes the structural problems listed under Caveats.

**Where it points** (darkest 2 % of the airglow layer within the 60-min isochrone;
coordinates are the darkest pixel of each patch, names approximate):

| Area (approx.) | Lat, Lon | Elev. | Best window centre | Zenith (mag) | 25° S / E (mag) | P(above fog) |
|---|---|---|---|---|---|---|
| Hills between Thun and the Emmental (Süderen–Buchholterberg) | 46.817, 7.781 | 1065 m | 102° (E) | 21.49 | 20.65 / 20.72 | 0.53 |
| Lower Diemtigtal | 46.608, 7.563 | 1252 m | 212° (SSW) | 21.59 | 20.78 / 20.73 | 0.84 |
| Gurnigel pass | 46.729, 7.454 | 1445 m | 198° (S) | 21.44 | 20.69 / 20.46 | 0.98 |
| Gantrisch / Selibüel ridge | 46.732, 7.401 | 1587 m | 188° (S) | 21.46 | 20.69 / 20.55 | 1.00 |
| Kiental | 46.597, 7.693 | 1355 m | 88° (E) | 21.52 | 20.68 / 20.73 | 0.94 |
| Engstligental (Frutigen–Adelboden) | 46.545, 7.597 | 1118 m | 78° (E) | 21.60 | 20.70 / 20.74 | 0.63 |

The Gurnigel/Gantrisch ridge is the only area that is dark, open to the S *and* almost
certainly above autumn stratus, at roughly 45–60 min from Bern. That is the "high
terrain toward the Mittelland, open toward S/E" optimum the brief asked about. The
differences between these areas (≈0.1–0.2 mag) are smaller than the model's absolute
uncertainty, so treat them as a shortlist for ground truth, not a ranking.

## Method

1. **Domain.** A 140 × 140 km LV95 square centred on Bern. Night-time drive-time isochrones
   (Valhalla on OSM, Saturday 23:00 departure → free-flow speeds) are exported as KML only;
   they do not mask anything.
2. **Terrain / horizons.** Copernicus GLO-30 DSM, reprojected to LV95 at 30 m, 3×3 median
   filtered (removes single-pixel masts, e.g. the Chasseral tower). For every pixel the DEM
   is ray-cast in 1° azimuth steps to 60 km (step max(30 m, s/150)), with Earth curvature
   and standard refraction (k = 4/3), eye height 1.6 m. This runs on a 100 m grid for the
   horizon-dependent fields and on a 250 m grid for the light-dome model.
3. **Light sources.** EOG VIIRS DNB VNL v2.2, 2025 annual `average_masked` (nW cm⁻² sr⁻¹),
   treated as an *emission source term* (radiance × pixel area ∝ nadir intensity), not as
   sky brightness. Sources lie within 250 km and are binned to 1 km LV95 cells, with 5 km
   blocks beyond 15 km. Each cell is softened as an area source (RMS radius cell/√6).
4. **Light-dome model.** Single-scattering line-of-sight integral (Garstang 1986 /
   Cinzano et al. 2000 type), tabulated once and looked up per source–pixel pair:
   - Rayleigh (0.0116 km⁻¹ at 550 nm, H = 8 km, λ⁻⁴) plus aerosol (fitted AOD, H = 1.5 km,
     Ångström 1.3, Henyey–Greenstein g = 0.65), with extinction on both legs;
   - the Garstang emission function normalised to the nadir value VIIRS sees,
     I(ψ)/I(0) = cos ψ + 0.554 F ψ⁴ / (2G(1−F)), F = G = 0.15;
   - Earth curvature for the line of sight and the source position, the source vertical
     tilted by D/R, and the observer–source height difference;
   - two bands: V (550 nm) and NIR (850 nm, the Gen 3 GaAs response region).
5. **Partial terrain blocking (two-component).** The kernel is binned by the elevation
   angle ε at which each scattering point sits as seen *from the source*. For every pixel,
   azimuth (5°) and source distance, the ray-caster stores the terrain *shadow angle* β:
   how high above its own horizon a ground source must emit to clear the terrain towards
   the observer. Scattering points with ε < β (the near-surface, blockable component) lose
   their direct illumination; points above (the high-altitude component) are untouched.
   SHADOW_LEAK = 15 % of the shadowed light is kept to stand in for multiple scattering.
   The low/high split is therefore set by geometry per source, not by a fixed fraction;
   see Results for the effective split it produces.
6. **Calibration.** With terrain switched off, as in the atlas, the modelled zenith
   brightness was fitted to the Lorenz 2025 atlas on every 250 m pixel. Free parameters:
   overall scale and aerosol optical depth (grid search over AOD, least squares in
   magnitudes).
7. **Bias corrections that cannot be calibrated** (the atlas is VIIRS-derived and shares
   them) are applied as stated priors:
   - **LED / blue light.** VIIRS DNB (500–900 nm) under-reads blue-rich LED. V-band
     artificial brightness is multiplied by 1.5 (range 1.2–2.0). For NIR the factor is 0.85
     (range 0.6–1.0), because white LEDs emit little beyond 700 nm while the HPS sodium
     they replaced emits strongly at 819 nm.
   - **Horizontal emission.** VIIRS looks at nadir, so the F ψ⁴ term is invisible to it.
     F is tested at 0.30, about twice Garstang's value.
8. **Scores.** These are separate per target, not blended (see the layer table).
   Inversion: P(above stratus top) = Φ((z − 1050 m)/200 m), an assumed distribution of
   Oct–Feb Swiss Plateau stratus-top heights. MeteoSwiss climatology was not ingested.

## Caveats

- **The input is VIIRS, so everything inherits its biases.**
  - *Spectral:* the DNB is blind below ~500 nm and under-reads blue-rich white LED. The
    V-band layers carry an assumed ×1.5 correction (range 1.2–2.0, ≈ +0.14 / −0.19 mag),
    applied uniformly. In reality the LED fraction varies by municipality, which would
    shift relative rankings; that is not modelled.
  - *Geometric:* VIIRS looks at nadir. Façades, billboards, floodlit sports grounds and
    ski-area lighting emit near-horizontally and are largely missed, and they build
    exactly the low-elevation domes that matter most for the airglow programme. The
    F = 0.30 run is a rough stand-in for this (+0.5 mag at zenith); local sources such as
    a floodlit pitch 3 km away, or a lit pass road, are not represented at all.
- **The calibration cannot fix the above.** The Lorenz atlas is itself VIIRS-driven and
  terrain-free, so the 0.04 mag fit shows the propagation model reproduces it, not that
  either matches the real sky. The LED factor, F and the shadow leak are priors.
- **Single scattering.** The kernel is single-scattering with a fixed 15 % shadow leak for
  multiple scattering. Aerosol properties are one fitted AOD (0.12) everywhere; real
  nights vary (AOD 0.03–0.3), which changes dome size and the zenith/horizon ratio.
- **NIR channel.** The NIR (850 nm) layer uses the same source map with an assumed ×0.85
  LED factor; the actual 800–900 nm emission of Swiss street lighting is unknown. It is
  reported *relative* to the best 1 % of pixels because the natural NIR sky (dominated by
  the OH emission you are observing) is not modelled.
- **Terrain.** Copernicus GLO-30 is a *surface* model: tree lines and buildings appear in
  the horizons, which is realistic for a pixel at a forest edge but means a 100 m pixel's
  horizon is not that of every point within it. 30 m resolution smooths sharp peaks (the
  Jungfrau sits ~0.2° low from Bern). Horizons are cast from the pixel centre at 1.6 m eye
  height. The shadow-angle table uses 28 distance nodes and 5° sectors, which leaves faint
  ring/texture artefacts (< 0.05 mag) around the brightest towns.
- **Inversion layer.** This is an *assumed* stratus-top distribution (normal, 1050 ± 200 m),
  not MeteoSwiss climatology. Real tops vary from ~700 m to > 1800 m between events and
  are lower in Jan–Feb than in Oct–Nov. Also, being above the fog *hides* the Mittelland
  lights beneath it, which this clear-sky model does not credit.
- **No access information.** Per your request the layers are not masked by roads, forest,
  land cover or drive time; only lakes are blanked. Check road access, winter closures
  (e.g. Gurnigel and the Gantrisch roads), parking and legality of tracks on the map.
- **Ground truth needed to confirm the ranking.**
  - *Instruments:* an SQM-L (narrow, ~20° FWHM) or, better, a calibrated all-sky camera
    (fisheye + photometric calibration) on several clear, moonless nights at a transect
    of 6–10 sites spanning the layer values: at least two of the shortlist above, two
    mid-range sites and one near Bern as an anchor.
  - *Targets 2–4:* measure the zenith *and* 25° elevation toward S and E.
  - *Airglow (target 1):* repeat the 25° S/E frames with the Gen 3 setup or a camera
    behind an ~800 nm long-pass filter. VIIRS-based modelling is weakest exactly there
    (LED spectrum, horizontal emitters).

## Reproducing

```bash
uv venv .venv --python 3.12
uv pip install --python .venv/bin/python numpy scipy rasterio geopandas shapely pyproj \
    matplotlib requests numba pandas pyarrow
# VIIRS needs a free EOG login: download
#   https://eogdata.mines.edu/nighttime_light/annual/v22/2025/VNL_npp_2025_global_vcmslcfg_v2_c202604011200.average_masked.dat.tif.gz
# into data/raw/viirs/ (or point config.VIIRS_FILE at any DNB radiance GeoTIFF)
.venv/bin/python -m darksky            # all steps, ~30 min on 14 cores
.venv/bin/python -m darksky layers     # just re-colour / re-export
```

Parameters live in `darksky/config.py`; colour ranges are in `darksky/build_layers.py`.

## Deviations from the original brief

- **Map layers instead of a candidate list** (per request). There is no candidate grid,
  OSM road/landcover filtering, top-15 table, per-site horizon plots or GeoJSON. Horizon
  profiles are computed for every pixel; the per-sector values are in the fields GeoTIFF.
- **Domain:** a ±70 km square around Bern instead of the isochrone polygon (per request).
  The isochrones are exported as KML.
- **Lorenz 2025 instead of 2022,** to match the 2025 VIIRS composite (2022 tiles are still
  downloaded by `fetch`).
- **Copernicus GLO-30 instead of swissALTI3D,** for one consistent DEM across the border
  out to 60 km. swissALTI3D (2 m tiles) was not used.
- **No MeteoSwiss stratus climatology;** an assumed distribution is used instead (see Caveats).
- **Overpass was unusable** (HTTP 406/504 for this area), so the earlier candidate version
  used the Geofabrik extract; that step is no longer needed.
