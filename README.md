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
4. Optional: import `outputs/isochrones_bern.kml` to overlay the 30/45/60-min night-time
   drive-time contours from Bern.

Colour keys are in `outputs/legend_*.png`. All brightness layers use one blue ramp:
**darker blue = darker sky (better)**.

| File | Target | Shows |
|---|---|---|
| `layer_1_airglow_SE.tif` | 1. OH airglow gravity waves | Artificial NIR (Gen 3, ~850 nm) skyglow in the darkest 30°-wide window at 20–30° elevation within azimuth 45–225° (E through S), in magnitudes relative to the best 1 % of the square. Grey hatching marks pixels where terrain blocks every E–S window. |
| `layer_2_zenith.tif` | 2. Hα nebulae, 3. deep sky | Zenith sky brightness, V band, LED-corrected (mag/arcsec²). |
| `layer_3_allsky30.tif` | **Main summary layer**; 3. Milky Way / binocular | Mean V sky brightness over the whole sky above 30° (mag/arcsec²), **including local sources** (settlements, jobs, buildings, larger roads at 100 m) **and direct glare** from lamps in view. **Shown only where terrain stays below 10° over ≥ 75 % of the horizon, the point is ≥ 200 m from any building, the 100 m pixel is < 50 % forest, and artificial light (incl. glare) is at most 2× the natural sky (all-sky ≥ 20.53)**; everything else is transparent. Stepped colour scale: lightest class = brighter than the threshold (unmasked layer only), then six equal 0.125 mag steps to the darkest sky; breaks 20.52 / 20.65 / 20.77 / 20.90 / 21.02 / 21.15. See [All-sky layer](#all-sky-layer-with-local-sources). |
| `layer_3b_allsky30_unmasked.tif` | Context for the main layer | Same values and 7-class scale as `layer_3_allsky30`, with no masks at all: every pixel of the square, lakes and rivers included. |
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

## Shareable map links

**[Sky above 30°, candidate sites + everywhere, with legends](https://map.geo.admin.ch/#/map?lang=en&center=2601000,1200000&z=2&bgLayer=ch.swisstopo.pixelkarte-grau&layers=WMTS%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/wmts/WMTSCapabilities.xml%7Cdarksky_allsky30_unmasked,f,0.65;WMTS%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/wmts/WMTSCapabilities.xml%7Cdarksky_allsky30,,0.8;KML%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/isochrones_bern.kml)** (recommended).
These two all-sky layers are served as a static WMTS (`outputs/wmts/`), so map.geo.admin.ch
shows their title, description and colour legend under the layer's ⓘ button in "Maps
displayed", like native layers. The candidate-site layer is on; the unmasked one is
loaded but hidden. To add them by hand: Advanced tools → Import → Online, paste
`https://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/wmts/WMTSCapabilities.xml`.

The COG links below load the same data as plain files; COG layers cannot show a legend
in the viewer.

Each link opens map.geo.admin.ch with the layer streamed from this repository (65 % opacity)
and the drive-time contours on top. Toggle or fade layers in the "Maps displayed" panel.

**[All layers in one map](https://map.geo.admin.ch/#/map?lang=en&center=2601000,1200000&z=2&bgLayer=ch.swisstopo.pixelkarte-grau&layers=COG%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/layer_1_airglow_SE.tif,,0.65;COG%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/layer_2_zenith.tif,f,0.65;COG%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/layer_3_allsky30.tif,f,0.65;COG%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/layer_3b_allsky30_unmasked.tif,f,0.65;COG%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/layer_4_meteors.tif,f,0.65;COG%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/layer_5_above_inversion.tif,f,0.65;COG%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/layer_ref_lorenz2025_zenith.tif,f,0.65;KML%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/isochrones_bern.kml)** (airglow shown, the others loaded but hidden)

| Layer | Link |
|---|---|
| OH airglow E–S (target 1) | [open](https://map.geo.admin.ch/#/map?lang=en&center=2601000,1200000&z=2&bgLayer=ch.swisstopo.pixelkarte-grau&layers=COG%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/layer_1_airglow_SE.tif,,0.65;KML%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/isochrones_bern.kml) |
| Zenith (targets 2/3) | [open](https://map.geo.admin.ch/#/map?lang=en&center=2601000,1200000&z=2&bgLayer=ch.swisstopo.pixelkarte-grau&layers=COG%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/layer_2_zenith.tif,,0.65;KML%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/isochrones_bern.kml) |
| Sky above 30° (target 3) | [open](https://map.geo.admin.ch/#/map?lang=en&center=2601000,1200000&z=2&bgLayer=ch.swisstopo.pixelkarte-grau&layers=COG%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/layer_3_allsky30.tif,,0.65;KML%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/isochrones_bern.kml) |
| Sky above 30°, unmasked | [open](https://map.geo.admin.ch/#/map?lang=en&center=2601000,1200000&z=2&bgLayer=ch.swisstopo.pixelkarte-grau&layers=COG%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/layer_3b_allsky30_unmasked.tif,,0.65;KML%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/isochrones_bern.kml) |
| Sky above 30°, masked + unmasked (toggle) | [open](https://map.geo.admin.ch/#/map?lang=en&center=2601000,1200000&z=2&bgLayer=ch.swisstopo.pixelkarte-grau&layers=COG%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/layer_3b_allsky30_unmasked.tif,f,0.65;COG%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/layer_3_allsky30.tif,,0.8;KML%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/isochrones_bern.kml) |
| Meteors (target 4) | [open](https://map.geo.admin.ch/#/map?lang=en&center=2601000,1200000&z=2&bgLayer=ch.swisstopo.pixelkarte-grau&layers=COG%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/layer_4_meteors.tif,,0.65;KML%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/isochrones_bern.kml) |
| Above autumn/winter stratus | [open](https://map.geo.admin.ch/#/map?lang=en&center=2601000,1200000&z=2&bgLayer=ch.swisstopo.pixelkarte-grau&layers=COG%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/layer_5_above_inversion.tif,,0.65;KML%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/isochrones_bern.kml) |
| Lorenz 2025 reference | [open](https://map.geo.admin.ch/#/map?lang=en&center=2601000,1200000&z=2&bgLayer=ch.swisstopo.pixelkarte-grau&layers=COG%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/layer_ref_lorenz2025_zenith.tif,,0.65;KML%7Chttps://raw.githubusercontent.com/KillianBrennan/dark_sky/main/outputs/isochrones_bern.kml) |

## Data sources and credits

- Night lights: VIIRS DNB VNL v2.2 annual 2025, Earth Observation Group, Payne Institute,
  Colorado School of Mines (CC BY 4.0). Elvidge et al. (2021), *Remote Sensing* 13, 922.
- Sky-brightness reference/calibration: World Atlas of the Artificial Night Sky Brightness,
  2025 edition, D. J. Lorenz (djlorenz.github.io/astronomy/lp), after Cinzano/Falchi.
- Terrain: Copernicus DEM GLO-30 © DLR e.V. 2010-2014 and © Airbus Defence and Space GmbH
  2014-2018, provided under COPERNICUS by the European Union and ESA.
- Water and forest masks: ESA WorldCover 10 m 2021 v200 (CC BY 4.0).
- Drive times: Valhalla (FOSSGIS public instance) on OpenStreetMap data © OpenStreetMap
  contributors (ODbL).

## All-sky layer with local sources

`layer_3_allsky30.tif` is built by `darksky/run_allsky.py`, separately from the other layers
(which still use VIIRS alone).

**Local light sources.** VIIRS has a ~0.5–1 km footprint and drops pixels below
~0.45 nW cm⁻² sr⁻¹. Inside Switzerland, BFS hectare statistics (residents STATPOP 2024, jobs
STATENT 2023, buildings 2024) and OSM roads say where light can be at 100 m:

- *Proxy.* Each hectare gets c·[residents, jobs, buildings, m of `lit=yes` road, m of
  untagged motorway/trunk/primary/secondary road]. The weights c ≥ 0 come from a
  non-negative least-squares fit of VIIRS over 5 km cells in Switzerland (log r = 0.94).
  Per unit, relative to one resident: job 0.99, building 3.3, lit road 0.05 per m,
  untagged major road 0.47 per m. Share of the total proxy: residents 38 %, jobs 19 %,
  buildings 25 %, untagged major roads 17 %, lit roads 1 %. The lit-road term is small
  because it is collinear with towns.
- *Downscaling.* The light of each Swiss VIIRS pixel is redistributed onto hectares within
  600 m (proxy × Gaussian, σ = 400 m), conserving the VIIRS total. 2.5 % of Swiss VIIRS
  light has no proxy nearby (roads through empty land, industry, ski lifts) and stays put.
- *Fill.* 54 k inhabited hectares with no lit VIIRS pixel within 600 m get their proxy
  intensity, capped so that no VIIRS pixel would exceed the detection floor. This adds
  1.8 % on top of Swiss VIIRS.
- *Abroad,* VIIRS pixels are used unchanged (softened as 385 m areas).

**Two-scale model.** Sources farther than 1.5 km go through the 250 m model (terrain
shadow from the ray-cast table). Hectare sources within 1.5 km are summed on the 100 m grid,
each with its own terrain shadow traced along the lamp → observer path in the DEM (lamp at
5 m). A smoothstep taper of ±300 m hands sources over, so each is counted once.
Recalibrated against Lorenz 2025 (terrain off, 250 m): AOD 0.12, RMS 0.047 mag.

**Sky sampling.** Zenith plus rings at 30°, 45° and 60° (12 azimuths each), weighted by
solid angle: bands 30–37.5°, 37.5–52.5°, 52.5–75°, 75–90° → 0.218, 0.369, 0.345, 0.068.
The previous version used only the zenith and a 45° ring. That underweighted the brighter
low part of the band and read about 0.13 mag too dark (`allsky30_minus_previous_layer`).

**Masks.**

- *Visibility:* the 5°-sector horizon maxima (Copernicus DSM, so trees count) must stay
  below 10° in at least 54 of 72 sectors.
- *Buildings:* at least 200 m from any OSM building outline (centroid distance minus
  outline radius, so slightly conservative). This uses 1.62 M buildings from the
  Switzerland, Franche-Comté, Alsace and Freiburg extracts, barns and sheds included.
- *Forest:* less than 50 % ESA WorldCover 2021 tree cover (10 m) within the 100 m pixel
  (`MAX_TREE_FRACTION`). The visibility test cannot catch closed forest, because the DSM
  puts the observer on top of the canopy.
- *Lakes:* blanked.

8.6 % of land pixels pass. 54 % fail visibility, 54 % fail the building distance and 42 %
fail the forest test, with much overlap. The forest mask alone removed 9 % of land that
passed the other two tests, i.e. about half of the previously shown pixels. What remains
is mostly open pasture, ridgelines and scattered single 100 m pixels.

**What the local sources change.** Within the shown pixels, very little:

| Variant (same sky sampling) | Δ mag vs central, 2 % / median / 98 % | Rank ρ |
|---|---|---|
| Raw VIIRS only | −0.025 / +0.002 / +0.049 | 0.9994 |
| No fill | +0.001 / +0.005 / +0.045 | 0.9997 |
| Fill × 2 | −0.022 / −0.002 / 0.000 | 0.9999 |

Sources within 1.5 km give a median 10 % of the artificial zenith light at shown pixels. Requiring 200 m from buildings already removes the places where
hectare-scale light dominates, so at a valid site the sky glow is still set by towns
kilometres away. Direct glare from lamps in view is a separate, larger local effect;
it is now included (see below).

**Direct glare.** A lamp in view scatters light inside the eye and lays a veil over the
retinal image that acts like extra background luminance. It is added to the sky
brightness using the CIE 146:2002 general disability-glare equation,
L_veil/E = 10/θ³ + (5/θ² + 0.1p/θ)[1 + (A/62.5)⁴] + 0.0025p, with age A = 40 and
pigmentation p = 0.5 (`CIE_AGE`, `CIE_PIGMENT`).

- *Illuminance at the eye:* E for each lamp within 6 km is computed in the sky model's own
  calibration, so the veil comes out in natural-sky units with no absolute photometry.
  It is scale × source weight × Garstang emission toward the eye (the near-horizontal
  ψ⁴ term VIIRS cannot see) × extinction / d².
- *Line of sight:* a lamp counts only if it is in direct view: per 100 m pixel, the
  running-maximum DSM angle along 72 rays at 11 distances (int8, 0.25° steps, 1.5 GB),
  with a ±0.25° visibility ramp.
- *Gaze:* the veil is averaged over the same gaze directions and weights as the sky
  brightness (rings at 30/45/60° and the zenith).
- *Size:* at shown pixels glare adds a median 10 % to the sky luminance (90th percentile
  41 %, 99th 98 %), i.e. −0.10 mag median, −0.37 mag at the 10th percentile and −0.74 mag
  at the 1st. The implied lamp output is a median ~300 cd per inhabited hectare towards
  the horizon (90th percentile ~1300 cd), plausible for a few streetlights plus façades
  and windows.

| Glare variant | Δ mag vs central, 2 % / median / 98 % | Rank ρ |
|---|---|---|
| Horizontal emission F = 0.30 (×2.4 lamp intensity) | −0.53 / −0.13 / 0.00 | 0.984 |
| Observer age 25 | 0.00 / +0.007 / +0.035 | 0.9999 |
| Observer age 60 | −0.155 / −0.033 / 0.00 | 0.998 |

The glare term is the least certain part of the model. It scales directly with the
near-horizontal lamp intensity, which VIIRS cannot measure, and OSM/BFS say nothing about
shielding. A fully shielded LED street lamp sends almost nothing toward a distant observer,
an old globe lamp a lot. The sky-only value is kept as `allsky30_sky_only_mag` and the
veil/sky ratio as `glare_veil_ratio` in the fields file.

**Brightness filter (masked layer).** Pixels are shown only if artificial light including
glare is at most `MAX_ARTIFICIAL_RATIO` = 2× the natural sky, i.e. all-sky ≥ 20.53 (the
pristine value for this sky average is 21.72). That is roughly a Bortle-4 sky in which the
Milky Way keeps clear structure. A stricter "artificial ≤ natural" (ratio 1, ≥ 20.97) leaves
only ~0.8 km² within 60 min of Bern: the Gantrisch ridge (~46.720 N, 7.446 E, 1630 m) and
the hills around Süderen–Chuderhüsi (~46.83–46.85 N, 7.76–7.79 E, 1100–1180 m). Both
already sit at about 21.0 with negligible glare. With ratio 2: 2.3 % of land is shown,
12.3 km² within 60 min and 2.3 km² within 45 min.

**Colour scale.** Shared by both layers. The lightest class is everything brighter than the
masked-layer threshold, so it appears only in the unmasked layer. The other six are equal
steps in magnitude (equal luminance ratios, i.e. logarithmic) from the threshold to the
darkest pixels (99.9th percentile): 0.125 mag ≈ 12 % in luminance per step, above the
model's ~0.05 mag relative resolution and the ±0.1 mag of an SQM. Steps never go below
`CLASS_MIN_WIDTH_MAG` = 0.1 mag; a stricter threshold would anchor the scale at the dark end
instead.

Data in `outputs/allsky30_fields_lv95.tif` (not in git; regenerate with
`python -m darksky allsky`): `allsky30_mag` (unmasked), `zenith_mag`, `local_share`,
`allsky30_sky_only_mag`, `glare_veil_ratio`, `visible_fraction_below10`, `building_dist_m`, `tree_fraction`,
`shown`, `allsky30_minus_previous_layer` (vs the original VIIRS-only layer).
Summaries: `outputs/allsky30_summary.json`, `allsky30_sensitivity.json`, `allsky30_glare_sensitivity.json`,
`calibration_allsky.json`, `sources_*.json`.

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
- **No access information.** The layers are not masked by roads, forest, land cover or
  drive time (only lakes are blanked; the all-sky layer also applies its visibility and
  building masks). Check road access, winter closures
  (e.g. Gurnigel and the Gantrisch roads), parking and legality of tracks on the map.
- **Local sources (all-sky layer).** The hectare proxy assumes every resident, job or
  building emits like the Swiss average; untagged motorways are treated as possibly lit,
  though most Swiss motorways are dark; OSM `lit` tags are incomplete. Abroad there is no
  hectare data. Glare assumes the Garstang near-horizontal emission for every lamp
  (no information on shielding) and the CIE eye-scatter function for a 40-year-old.
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
.venv/bin/python -m darksky            # all steps, ~40 min on 14 cores
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
