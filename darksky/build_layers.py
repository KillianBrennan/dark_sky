"""Step 6 driver: write data GeoTIFF, RGBA COG layers + legends, sensitivity summary."""
import json

import numpy as np

from . import config as C
from . import layers as Ly

# (file name, field, vmin, vmax, dark_is_high, legend title, ticks, tick labels, note)
LAYERS = [
    ("1_airglow_SE", "airglow_excess", 0.0, 3.0, False,
     "1  OH airglow, 20-30° elev., darkest 30° window E-S (NIR, Gen3)",
     [0, 0.5, 1, 1.5, 2, 2.5, 3], ["best", "+0.5", "+1", "+1.5", "+2", "+2.5", "≥ +3 mag"],
     "Artificial NIR skyglow relative to the best 1 % of the area. Grey hatch: terrain above 15° in every 30° E–S window."),
    ("2_zenith", "zenith_mag", 19.8, 21.7, True,
     "2/3  Zenith sky brightness, V (mag/arcsec²), LED-corrected",
     [19.8, 20.2, 20.6, 21.0, 21.4, 21.7], None, None),
    ("3_allsky30", "allsky30_mag", 19.4, 21.4, True,
     "3  Mean sky brightness above 30°, V (mag/arcsec²), LED-corrected",
     [19.4, 19.8, 20.2, 20.6, 21.0, 21.4], None, None),
    ("4_meteors", "meteor_rel", 0.2, 0.75, True,
     "4  Relative visual meteor rate (1 = pristine sky, full horizon)",
     [0.2, 0.3, 0.4, 0.5, 0.6, 0.7], None, None),
    ("ref_lorenz2025_zenith", "lorenz_zenith_mag", 19.8, 21.7, True,
     "Reference: Lorenz 2025 atlas zenith (mag/arcsec²), no terrain, no LED corr.",
     [19.8, 20.2, 20.6, 21.0, 21.4, 21.7], None, None),
]


def _rank(a):
    return np.argsort(np.argsort(a))


def sensitivity(central):
    out = {}
    m0 = np.isfinite(central["airglow_excess"])
    top0 = central["airglow_excess"] <= np.nanquantile(central["airglow_excess"], 0.05)
    for tag in ("F030", "leak005", "leak035"):
        f, *_ = Ly.compute(tag)
        a = f["airglow_excess"]
        m = m0 & np.isfinite(a)
        top = a <= np.nanquantile(a, 0.05)
        mz = np.isfinite(central["zenith_mag"])
        out[tag] = dict(
            spearman_airglow=float(np.corrcoef(_rank(central["airglow_excess"][m]), _rank(a[m]))[0, 1]),
            top5pct_overlap=float((top & top0).sum() / top0.sum()),
            spearman_zenith=float(np.corrcoef(_rank(central["zenith_mag"][mz]), _rank(f["zenith_mag"][mz]))[0, 1]),
            median_zenith_shift_mag=float(np.nanmedian(f["zenith_mag"] - central["zenith_mag"])),
        )
    # LED factor range: uniform multiplicative -> a brightness shift, no rank change
    zen_ratio = 10 ** ((C.NATURAL_ZENITH_MPSAS - central["zenith_mag_noLED"]) / 2.5) - 1
    led_c, led_lo, led_hi = C.LED_FACTOR["V"]
    out["LED_V_range"] = {
        f"{lab}": float(np.nanmedian(Ly.ratio_to_mpsas(fac * zen_ratio) - central["zenith_mag"]))
        for lab, fac in (("x1.2", led_lo), ("x2.0", led_hi))}
    return out


def run():
    f, g, tr, sh = Ly.compute("central")
    print(Ly.write_fields(f, g, tr, sh))
    for name, field, vmin, vmax, dark_hi, title, ticks, labels, note in LAYERS:
        blocked = (~np.isfinite(f[field]) & np.isfinite(f["zenith_mag"])) if field == "airglow_excess" else None
        p = Ly.write_rgba(name, f[field], g, tr, sh, vmin, vmax, dark_is_high=dark_hi, hatch_nan=blocked)
        Ly.legend(name, title, vmin, vmax, ticks, labels, note, dark_is_high=dark_hi)
        print(p)
    p = Ly.write_rgba("5_above_inversion", f["p_above_inversion"], g, tr, sh, 0, 1,
                      alpha_from=f["p_above_inversion"])
    print(p)
    sens = sensitivity(f)
    (C.OUT / "sensitivity.json").write_text(json.dumps(sens, indent=2))
    print(json.dumps(sens, indent=2))
    return f, g, tr, sh


if __name__ == "__main__":
    run()
