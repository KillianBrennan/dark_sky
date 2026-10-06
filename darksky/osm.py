"""OpenStreetMap buildings (for the 200 m proximity mask) and larger roads (light proxy).

Extracts (Geofabrik): Switzerland, plus Franche-Comte, Alsace and Regierungsbezirk Freiburg
so that the foreign edges of the square are masked too. Buildings: every way tagged
building=* except ruins/demolished/construction, stored as centroid + outline radius
(LV95). Multipolygon buildings (relations, a small minority) are not included.
Roads (Swiss extract only, used in the light-source proxy): motorway, trunk, primary,
secondary (+ links) and tertiary, with their lit=* tag, rasterised to metres per BFS
hectare in two classes: lit=yes, and untagged major roads (motorway..secondary).
lit=no roads are dropped.
"""
import math

import numpy as np
import pandas as pd
import requests
from pyproj import Transformer

from . import config as C
from .fetch import domain_lv95

EXTRACTS = {"switzerland": "europe/switzerland", "franche-comte": "europe/france/franche-comte",
            "alsace": "europe/france/alsace", "freiburg-regbez": "europe/germany/baden-wuerttemberg/freiburg-regbez"}
SKIP_BUILDING = {"ruins", "demolished", "construction", "no", "proposed"}
ROAD_LIT_CLASSES = {"motorway", "motorway_link", "trunk", "trunk_link", "primary", "primary_link",
                    "secondary", "secondary_link", "tertiary", "tertiary_link"}
ROAD_MAJOR = {"motorway", "motorway_link", "trunk", "trunk_link", "primary", "primary_link",
              "secondary", "secondary_link"}


def _pbf(name):
    d = C.RAW / "osm"
    d.mkdir(exist_ok=True)
    p = d / f"{name}.osm.pbf"
    if not p.exists():
        with requests.get(f"https://download.geofabrik.de/{EXTRACTS[name]}-latest.osm.pbf", stream=True, timeout=600) as r:
            r.raise_for_status()
            with open(p, "wb") as f:
                for chunk in r.iter_content(1 << 20):
                    f.write(chunk)
    return p


def _ll_bounds(margin_m):
    minx, miny, maxx, maxy = domain_lv95()
    t = Transformer.from_crs(C.CRS_M, C.CRS_LL, always_xy=True)
    lon, lat = t.transform([minx - margin_m, maxx + margin_m, minx - margin_m, maxx + margin_m],
                           [miny - margin_m, miny - margin_m, maxy + margin_m, maxy + margin_m])
    return min(lon), min(lat), max(lon), max(lat)


def buildings():
    out = C.INTERIM / "osm_buildings.parquet"
    if out.exists():
        return pd.read_parquet(out)
    import osmium
    lo0, la0, lo1, la1 = _ll_bounds(500)
    rows = []
    for name in EXTRACTS:
        fp = osmium.FileProcessor(str(_pbf(name))).with_locations().with_filter(osmium.filter.KeyFilter("building"))
        for o in fp:
            if not o.is_way() or o.tags.get("building") in SKIP_BUILDING:
                continue
            try:
                lon = np.array([n.lon for n in o.nodes])
                lat = np.array([n.lat for n in o.nodes])
            except osmium.InvalidLocationError:
                continue
            cx, cy = lon.mean(), lat.mean()
            if not (lo0 <= cx <= lo1 and la0 <= cy <= la1):
                continue
            k = 111320.0 * math.cos(math.radians(cy))
            r = float(np.sqrt(((lon - cx) * k) ** 2 + ((lat - cy) * 110540.0) ** 2).max())
            rows.append((o.id, cx, cy, r, o.tags.get("building")))
    df = pd.DataFrame(rows, columns=["osm_id", "lon", "lat", "radius_m", "building"]).drop_duplicates("osm_id")
    df["x"], df["y"] = Transformer.from_crs(C.CRS_LL, C.CRS_M, always_xy=True).transform(df.lon.values, df.lat.values)
    df.to_parquet(out)
    return df


def roads_per_hectare():
    """Metres of road per BFS hectare (SW-corner keyed like the BFS tables), Swiss extract."""
    out = C.INTERIM / "osm_roads_hectare.parquet"
    if out.exists():
        return pd.read_parquet(out)
    import osmium
    t = Transformer.from_crs(C.CRS_LL, C.CRS_M, always_xy=True)
    xs, ys, cls = [], [], []
    fp = osmium.FileProcessor(str(_pbf("switzerland"))).with_locations().with_filter(osmium.filter.KeyFilter("highway"))
    stats = {"lit_yes": 0, "lit_no": 0, "lit_untagged_major": 0, "lit_untagged_tertiary": 0}
    for o in fp:
        if not o.is_way():
            continue
        hw = o.tags.get("highway")
        if hw not in ROAD_LIT_CLASSES:
            continue
        lit = o.tags.get("lit")
        if lit in ("no", "disused"):
            stats["lit_no"] += 1
            continue
        if lit not in (None, "no"):
            c = 0                                        # lit=yes / automatic / limited / interval ...
            stats["lit_yes"] += 1
        elif hw in ROAD_MAJOR:
            c = 1
            stats["lit_untagged_major"] += 1
        else:
            stats["lit_untagged_tertiary"] += 1
            continue                                     # untagged tertiary: no information, skip
        try:
            lon = np.array([n.lon for n in o.nodes])
            lat = np.array([n.lat for n in o.nodes])
        except osmium.InvalidLocationError:
            continue
        x, y = t.transform(lon, lat)
        # sample every ~10 m along the line
        seg = np.hypot(np.diff(x), np.diff(y))
        for i in np.nonzero(seg > 0)[0]:
            n = max(int(seg[i] // 10.0), 1)
            f = (np.arange(n) + 0.5) / n
            xs.append(x[i] + f * (x[i + 1] - x[i]))
            ys.append(y[i] + f * (y[i + 1] - y[i]))
            cls.append(np.full(n, c, np.int8))
    x = np.concatenate(xs)
    y = np.concatenate(ys)
    c = np.concatenate(cls)
    seglen = 10.0  # each sample represents ~10 m (exact per segment: seg/n; error < 10 %)
    df = pd.DataFrame({"E_KOORD": (x // 100 * 100).astype(np.int64), "N_KOORD": (y // 100 * 100).astype(np.int64),
                       "road_lit": np.where(c == 0, seglen, 0.0), "road_major": np.where(c == 1, seglen, 0.0)})
    df = df.groupby(["E_KOORD", "N_KOORD"], as_index=False).sum()
    df.attrs["way_counts"] = stats
    df.to_parquet(out)
    print("road ways:", stats)
    return df


if __name__ == "__main__":
    import time
    t = time.time()
    r = roads_per_hectare()
    print(f"roads: {len(r)} hectares, lit {r.road_lit.sum() / 1e3:.0f} km, untagged major {r.road_major.sum() / 1e3:.0f} km, {time.time() - t:.0f}s", flush=True)
    t = time.time()
    b = buildings()
    print(f"buildings: {len(b)} in domain, median radius {b.radius_m.median():.1f} m, {time.time() - t:.0f}s")
