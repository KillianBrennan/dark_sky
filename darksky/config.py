"""Shared configuration. All tunable assumptions live here so they are explicit."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"
OUT = ROOT / "outputs"
for _p in (RAW, INTERIM, OUT):
    _p.mkdir(parents=True, exist_ok=True)

BERN = (46.9480, 7.4474)  # lat, lon
MAX_DRIVE_MIN = 60
# Saturday 23:00 departure -> Valhalla uses night / free-flow speeds where available
VALHALLA_URL = "https://valhalla1.openstreetmap.de"
DEPART = "2026-10-10T23:00"

# Analysis domain: square of +/- DOMAIN_HALF_KM around Bern (LV95), no drive-time mask.
# The isochrones are still computed and exported as an optional KML overlay.
DOMAIN_HALF_KM = 70

# Metric working CRS: Swiss LV95
CRS_M = "EPSG:2056"
CRS_LL = "EPSG:4326"

# ---- horizon ----
HORIZON_MAX_KM = 60
AZ_STEP_DEG = 1
DEM_RES_M = 30              # working DEM resolution (Copernicus GLO-30 native ~30 m)
K_REFRACTION = 4.0 / 3.0    # effective Earth radius factor
R_EARTH = 6_371_000.0

# ---- inversion (Mittelland stratus) ----
# Assumed distribution of Oct-Feb stratus/fog top heights over the Swiss Plateau,
# used when MeteoSwiss climatology is not ingested. Literature (e.g. Bendix 2002;
# Scherrer & Appenzeller 2014; MeteoSwiss Hochnebel climatology) places typical
# tops at 800-1200 m a.s.l., median ~1000 m, occasional tops to 1500-1800 m.
# Modelled as normal(mean, sd) truncated to >= 600 m.
INVERSION_TOP_MEAN_M = 1050
INVERSION_TOP_SD_M = 200

# ---- light dome model (single-scattering line-of-sight integral, Garstang-type) ----
BANDS = {"V": 550.0, "NIR": 850.0}      # nm; NIR ~ Gen3 GaAs response peak region
RAYLEIGH_BETA0_550 = 0.0116e-3          # 1/m at sea level, 550 nm
H_MOL = 8000.0                          # m
AOD_550 = 0.10                          # prior; refit in calibration
H_AER = 1500.0                          # m, aerosol scale height
ANGSTROM = 1.3
HG_G = 0.65                             # aerosol asymmetry parameter
GARSTANG_F = 0.15                       # fraction emitted near-horizontally (Garstang 1986)
GARSTANG_G = 0.15                       # ground reflectance
SOURCE_ALT = 450.0                      # m, reference source altitude for the kernel table
SHADOW_LEAK = 0.15                      # fraction of terrain-shadowed single-scatter light that still
                                        # arrives (multiple scattering fills shadows; not modelled)
ELEVATIONS = (15, 20, 25, 30, 45, 90)   # deg, output elevations
AZ_BIN = 5                              # deg, output azimuth bins

# ---- VIIRS bias handling ----
# VIIRS DNB (500-900 nm) under-reads blue-rich LED. Visual (V) skyglow factor applied
# on top of the Lorenz-calibrated model: central 1.5, range 1.2-2.0 (motivated by
# Sanchez de Miguel et al. 2021; Kyba et al. 2023 report ~7-10 %/yr ground-based
# brightening vs ~2 %/yr in VIIRS). For the NIR (Gen3, ~850 nm) channel LEDs emit
# little beyond 700 nm while HPS has strong 819 nm Na lines, so the factor is <1:
# central 0.85, range 0.6-1.0. These are priors, NOT fitted: the Lorenz atlas is
# itself VIIRS-derived and shares the same blindness, so calibration cannot see it.
LED_FACTOR = {"V": (1.5, 1.2, 2.0), "NIR": (0.85, 0.6, 1.0)}
NATURAL_ZENITH_MPSAS = 22.0             # Lorenz convention
V_EXTINCTION_K = 0.17                   # mag/airmass, for natural-sky brightening off zenith
VIIRS_MIN_RADIANCE = 0.0               # nW/cm2/sr; the *masked* VNL product already zeroes background
SOURCE_RADIUS_KM = 250
VIIRS_FILE = RAW / "viirs" / "VNL_npp_2025_global_vcmslcfg_v2_c202604011200.average_masked.dat.tif.gz"
LORENZ_YEAR = 2025                      # same year as the VIIRS composite -> like-for-like calibration

# ---- local sources from BFS hectare statistics (STATPOP / STATENT / building statistics) ----
BFS_FILES = {
    "pop": "volkszaehlung-bevoelkerungsstatistik_einwohner_2024_2056.csv",
    "jobs": "betriebszaehlungen-beschaeftigte_vollzeitaequivalente_2023_2056.csv",
    "bldg": "volkszaehlung-gebaeudestatistik_gebaeude_2024_2056.csv",
}
BFS_URL = "https://data.geo.admin.ch/{coll}/{item}/{file}"
DOWNSCALE_SIGMA_M = 400.0     # VIIRS footprint + geolocation: Gaussian redistribution scale
DOWNSCALE_RADIUS_M = 600.0
VIIRS_DETECTION_FLOOR = 0.45  # nW/cm2/sr, ~1st percentile of lit VNL pixels in the domain
MISSED_FILL_FACTOR = 1.0      # scale on the per-unit intensity given to VIIRS-dark inhabited hectares
LOCAL_RADIUS_M = 1500.0       # sources nearer than this are summed at 100 m with exact terrain
LOCAL_TAPER_M = 300.0         # half-width of the smooth hand-over between local and far terms

# ---- all-sky layer masks ----
VIS_MAX_HORIZON_DEG = 10.0    # "good visibility": terrain below this ...
VIS_MIN_FRACTION = 0.75       # ... over at least this fraction of the horizon (5-deg sector maxima)
BUILDING_MIN_DIST_M = 200.0   # minimum distance to any OSM building outline
