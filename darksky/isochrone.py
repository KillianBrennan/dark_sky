"""Step 1: 60-min night-time driving isochrone from Bern (Valhalla, OSM data)."""
import json
import requests
from . import config as C


def fetch(contours=(15, 30, 45, 60)):
    feats = []
    if (C.INTERIM / "isochrones.geojson").exists():
        return C.INTERIM / "isochrones.geojson"
    # Valhalla allows max 4 contours per request
    for chunk in (contours[:4], contours[4:]):
        if not chunk:
            continue
        q = {
            "locations": [{"lat": C.BERN[0], "lon": C.BERN[1]}],
            "costing": "auto",
            "contours": [{"time": t} for t in chunk],
            "polygons": True,
            "denoise": 0.3,
            "generalize": 50,
            "date_time": {"type": 1, "value": C.DEPART},
        }
        r = requests.post(f"{C.VALHALLA_URL}/isochrone", data=json.dumps(q, separators=(",", ":")), timeout=120)
        r.raise_for_status()
        feats += r.json()["features"]
    fc = {"type": "FeatureCollection", "features": sorted(feats, key=lambda f: f["properties"]["contour"])}
    out = C.INTERIM / "isochrones.geojson"
    out.write_text(json.dumps(fc))
    return out


def to_kml(levels=(30, 45, 60)):
    """Drive-time contours as KML lines (importable in map.geo.admin.ch next to the layers).
    45 min is not a Valhalla contour we request, so it is skipped unless present."""
    from xml.sax.saxutils import escape
    import geopandas as gpd
    g = gpd.read_file(C.INTERIM / "isochrones.geojson")
    colours = {30: "ff2678eb", 45: "ff2678eb", 60: "ff1a1a1a", 20: "ff2678eb", 40: "ff2678eb", 50: "ff2678eb"}
    marks = []
    for lev in levels:
        row = g[g.contour == lev]
        if row.empty:
            continue
        geom = row.geometry.iloc[0].boundary
        lines = getattr(geom, "geoms", [geom])
        for ln in lines:
            if ln.length < 0.02:          # drop tiny islands of the generalised polygon
                continue
            coords = " ".join(f"{x:.5f},{y:.5f}" for x, y in ln.coords)
            marks.append(f"""<Placemark><name>{escape(f"{lev} min drive from Bern (night)")}</name>
<Style><LineStyle><color>{colours.get(lev, "ff2678eb")}</color><width>{3 if lev == 60 else 2}</width></LineStyle></Style>
<LineString><tessellate>1</tessellate><coordinates>{coords}</coordinates></LineString></Placemark>""")
    kml = ('<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2"><Document>'
           "<name>Drive time from Bern (Valhalla, Sat 23:00)</name>" + "".join(marks) + "</Document></kml>")
    out = C.OUT / "isochrones_bern.kml"
    out.write_text(kml)
    return out


if __name__ == "__main__":
    print(fetch())
    print(to_kml())
