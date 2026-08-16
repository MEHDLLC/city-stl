#!/usr/bin/env python3
"""
city_stl.py — Turn any city into a 3D-printable STL with terrain + buildings.

Usage:
    python city_stl.py "Seattle, WA" --x-km 4 --y-km 3 --width-mm 150 -o seattle.stl
    python city_stl.py "10001"       --x-km 4 --y-km 3 --width-mm 150      # US ZIP code
    python city_stl.py "40.7484, -73.9857" --x-km 4 --y-km 3 --width-mm 150   # lat, lon

Options:
    --z-exag 2.0        vertical exaggeration for terrain (buildings scale realistically)
    --base-mm 3         thickness of solid base under lowest terrain point
    --default-levels 2  floors to assume when OSM has no height data
    --min-area 20       ignore building footprints smaller than this (m²)
    --dem-zoom 13       terrain tile zoom (12≈38m/px, 13≈19m/px, 14≈10m/px, 15≈5m/px)
    --no-buildings      terrain only
    --sea-level 0       elevations below this (m) are flattened to it (hides seafloor data)

Data: OpenStreetMap (buildings, via osmnx) and AWS Terrain Tiles (elevation, no key needed).
Install: pip install osmnx trimesh shapely numpy scipy pillow requests pyproj mapbox-earcut
"""

import argparse
import math
import re
import sys
import numpy as np
import requests
from io import BytesIO
from PIL import Image

import osmnx as ox
import trimesh
from shapely.geometry import Polygon, MultiPolygon, box
from shapely.geometry.polygon import orient
from pyproj import Transformer
import mapbox_earcut as earcut

TILE_URL = "https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png"


# ----------------------------------------------------------------------------- geo helpers
def utm_epsg(lat, lon):
    zone = int((lon + 180) / 6) + 1
    return 32600 + zone if lat >= 0 else 32700 + zone


def lonlat_to_tile(lon, lat, z):
    n = 2 ** z
    x = (lon + 180) / 360 * n
    lat_r = math.radians(lat)
    y = (1 - math.log(math.tan(lat_r) + 1 / math.cos(lat_r)) / math.pi) / 2 * n
    return x, y


def tile_to_lonlat(x, y, z):
    n = 2 ** z
    lon = x / n * 360 - 180
    lat = math.degrees(math.atan(math.sinh(math.pi * (1 - 2 * y / n))))
    return lon, lat


# ----------------------------------------------------------------------------- terrain
def fetch_dem(west, south, east, north, zoom):
    """Return (heights[m] as 2D array, lon_axis, lat_axis) covering the bbox."""
    x0, y1 = lonlat_to_tile(west, south, zoom)
    x1, y0 = lonlat_to_tile(east, north, zoom)
    tx0, tx1 = int(math.floor(x0)), int(math.floor(x1))
    ty0, ty1 = int(math.floor(y0)), int(math.floor(y1))
    ntx, nty = tx1 - tx0 + 1, ty1 - ty0 + 1
    print(f"  fetching {ntx * nty} terrain tiles at zoom {zoom}...")
    mosaic = np.zeros((nty * 256, ntx * 256), dtype=np.float32)
    for tx in range(tx0, tx1 + 1):
        for ty in range(ty0, ty1 + 1):
            r = requests.get(TILE_URL.format(z=zoom, x=tx, y=ty), timeout=30)
            r.raise_for_status()
            img = np.asarray(Image.open(BytesIO(r.content)).convert("RGB"), dtype=np.float32)
            h = img[..., 0] * 256 + img[..., 1] + img[..., 2] / 256 - 32768
            mosaic[(ty - ty0) * 256:(ty - ty0 + 1) * 256, (tx - tx0) * 256:(tx - tx0 + 1) * 256] = h

    # pixel -> lon/lat axes, then crop to bbox
    px_x = tx0 + np.arange(mosaic.shape[1] + 1) / 256
    px_y = ty0 + np.arange(mosaic.shape[0] + 1) / 256
    lons = np.array([tile_to_lonlat(px, 0, zoom)[0] for px in px_x[:-1] + 0.5 / 256])
    lats = np.array([tile_to_lonlat(0, py, zoom)[1] for py in px_y[:-1] + 0.5 / 256])
    cx = (lons >= west) & (lons <= east)
    cy = (lats >= south) & (lats <= north)
    dem = mosaic[np.ix_(cy, cx)]
    return dem, lons[cx], lats[cy]


def terrain_mesh(dem, xs, ys, base_z):
    """Watertight solid: heightmap top, flat bottom at base_z, four walls.
    dem[i, j] corresponds to (xs[j], ys[i]) in projected metres."""
    ny, nx = dem.shape
    X, Y = np.meshgrid(xs, ys)
    top = np.column_stack([X.ravel(), Y.ravel(), dem.ravel()])
    bot = np.column_stack([X.ravel(), Y.ravel(), np.full(dem.size, base_z)])
    verts = np.vstack([top, bot])
    idx = np.arange(nx * ny).reshape(ny, nx)
    off = nx * ny
    faces = []
    a, b, c, d = idx[:-1, :-1], idx[:-1, 1:], idx[1:, :-1], idx[1:, 1:]
    faces.append(np.column_stack([a.ravel(), c.ravel(), b.ravel()]))
    faces.append(np.column_stack([b.ravel(), c.ravel(), d.ravel()]))
    faces.append(np.column_stack([a.ravel(), b.ravel(), c.ravel()]) + off)
    faces.append(np.column_stack([b.ravel(), d.ravel(), c.ravel()]) + off)

    def wall(line, flip):
        f = []
        for i in range(len(line) - 1):
            p, q = line[i], line[i + 1]
            t1, t2 = (p, q, q + off), (p, q + off, p + off)
            if flip:
                t1, t2 = t1[::-1], t2[::-1]
            f += [t1, t2]
        return np.array(f)

    faces.append(wall(idx[0, :], True))          # south edge (y min)
    faces.append(wall(idx[-1, :], False))        # north edge
    faces.append(wall(idx[:, 0], False))         # west edge
    faces.append(wall(idx[:, -1], True))         # east edge
    m = trimesh.Trimesh(vertices=verts, faces=np.vstack(faces), process=True)
    m.fix_normals()
    return m


# ----------------------------------------------------------------------------- buildings
def building_height(row, default_levels):
    for key in ("height", "building:height"):
        v = row.get(key)
        if isinstance(v, str):
            try:
                return float(v.split()[0].replace("m", ""))
            except ValueError:
                pass
        elif isinstance(v, (int, float)) and not np.isnan(v):
            return float(v)
    lv = row.get("building:levels")
    try:
        lv = float(lv)
        if not np.isnan(lv):
            return lv * 3.2 + 1.0
    except (TypeError, ValueError):
        pass
    return default_levels * 3.2 + 1.0


def fetch_buildings(bbox_wgs, epsg, clip_box, default_levels, min_area):
    west, south, east, north = bbox_wgs
    print("  fetching OSM buildings...")
    try:
        gdf = ox.features_from_bbox((west, south, east, north), tags={"building": True})
    except Exception as e:
        print(f"  no buildings returned ({e})")
        return []
    gdf = gdf[gdf.geometry.type.isin(["Polygon", "MultiPolygon"])]
    if "building" in gdf.columns:
        gdf = gdf[gdf["building"].astype(str) != "no"]
    gdf = gdf.to_crs(epsg=epsg)
    geoms = gdf.geometry.intersection(clip_box)
    cols = {c: gdf[c].tolist() for c in ("height", "building:height", "building:levels") if c in gdf.columns}
    n = len(gdf)
    out = []
    for k in range(n):
        geom = geoms.iloc[k]
        if geom.is_empty:
            continue
        row = {c: v[k] for c, v in cols.items()}
        h = building_height(row, default_levels)
        polys = list(geom.geoms) if isinstance(geom, MultiPolygon) else [geom]
        for p in polys:
            if isinstance(p, Polygon) and p.area >= min_area:
                out.append((p, h))
    print(f"  {len(out)} building footprints kept")
    return out


def sample_dem(dem, xs, ys, px, py):
    j = np.clip(np.searchsorted(xs, px), 0, len(xs) - 1)
    i = np.clip(np.searchsorted(ys, py), 0, len(ys) - 1)
    return dem[i, j]


def _extrude(poly, z0, z1):
    """Fast prism from a shapely polygon (with holes). Returns (verts, faces)."""
    poly = orient(poly, -1.0)          # exterior CW, holes CCW -> outward normals
    rings = [np.asarray(poly.exterior.coords)[:-1]] + \
            [np.asarray(r.coords)[:-1] for r in poly.interiors]
    rings = [r for r in rings if len(r) >= 3]
    if not rings:
        return None
    pts = np.vstack(rings)
    ends = np.cumsum([len(r) for r in rings]).astype(np.uint32)
    tri = earcut.triangulate_float64(pts, ends).reshape(-1, 3)
    if len(tri) == 0:
        return None
    n = len(pts)
    top = np.column_stack([pts, np.full(n, z1)])
    bot = np.column_stack([pts, np.full(n, z0)])
    verts = np.vstack([top, bot])
    faces = [tri, tri[:, ::-1] + n]                     # roof (up), floor (down)
    off = 0
    for r in rings:
        m = len(r)
        i = np.arange(m) + off
        j = (np.arange(m) + 1) % m + off
        # wall quads: top_i, top_j, bot_j / top_i, bot_j, bot_i
        faces.append(np.column_stack([i, j, j + n]))
        faces.append(np.column_stack([i, j + n, i + n]))
        off += m
    return verts, np.vstack(faces)


def buildings_mesh(bldgs, dem, xs, ys, z_exag):
    """Batch all buildings into one mesh with numpy — orders of magnitude faster
    than one trimesh object per building."""
    V, F, base = [], [], 0
    for poly, h in bldgs:
        c = poly.centroid
        gz = sample_dem(dem, xs, ys, c.x, c.y) * z_exag - 2.0   # sink 2 m into terrain
        r = _extrude(poly.simplify(0.3), gz, gz + h)
        if r is None:
            continue
        v, f = r
        V.append(v)
        F.append(f + base)
        base += len(v)
    if not V:
        return []
    return [trimesh.Trimesh(vertices=np.vstack(V), faces=np.vstack(F), process=False)]


# ----------------------------------------------------------------------------- location
def resolve_location(text):
    """Accepts a place name, a US ZIP code, or 'lat, lon' decimal coordinates."""
    t = text.strip()
    m = re.fullmatch(r"\s*(-?\d+(?:\.\d+)?)\s*[, ]\s*(-?\d+(?:\.\d+)?)\s*", t)
    if m:
        lat, lon = float(m.group(1)), float(m.group(2))
        if -90 <= lat <= 90 and -180 <= lon <= 180:
            print(f"Using coordinates {lat}, {lon}")
            return lat, lon
    if re.fullmatch(r"\d{5}(-\d{4})?", t):
        t = f"{t}, USA"                       # bare ZIP -> help the geocoder
    print(f"Geocoding '{t}'...")
    return ox.geocode(t)


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("city", help="place name, US ZIP code, or 'lat, lon'")
    ap.add_argument("--x-km", type=float, required=True)
    ap.add_argument("--y-km", type=float, required=True)
    ap.add_argument("--width-mm", type=float, required=True, help="printed size along X")
    ap.add_argument("-o", "--output", default="city.stl")
    ap.add_argument("--z-exag", type=float, default=2.0)
    ap.add_argument("--base-mm", type=float, default=3.0)
    ap.add_argument("--default-levels", type=float, default=2)
    ap.add_argument("--min-area", type=float, default=20)
    ap.add_argument("--dem-zoom", type=int, default=13)
    ap.add_argument("--no-buildings", action="store_true")
    ap.add_argument("--sea-level", type=float, default=0.0,
                    help="clamp elevations below this (m) — flattens ocean/bathymetry")
    a = ap.parse_args()

    # 1. resolve location & bbox
    lat, lon = resolve_location(a.city)
    epsg = utm_epsg(lat, lon)
    to_utm = Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True)
    to_wgs = Transformer.from_crs(f"EPSG:{epsg}", "EPSG:4326", always_xy=True)
    cx, cy = to_utm.transform(lon, lat)
    hx, hy = a.x_km * 500, a.y_km * 500
    xmin, xmax, ymin, ymax = cx - hx, cx + hx, cy - hy, cy + hy
    west, south = to_wgs.transform(xmin, ymin)
    east, north = to_wgs.transform(xmax, ymax)
    print(f"  center {lat:.4f},{lon:.4f}  bbox {a.x_km}x{a.y_km} km  UTM EPSG:{epsg}")

    # 2. terrain
    dem_ll, lons, lats = fetch_dem(west, south, east, north, a.dem_zoom)
    # resample onto a regular metric grid
    res = (xmax - xmin) / max(dem_ll.shape[1] - 1, 1)
    xs = np.arange(xmin, xmax + res / 2, res)
    ys = np.arange(ymin, ymax + res / 2, res)
    X, Y = np.meshgrid(xs, ys)
    LON, LAT = to_wgs.transform(X, Y)
    li = np.clip(np.searchsorted(lons, LON), 0, len(lons) - 1)
    lj = np.clip(np.searchsorted(-lats, -LAT), 0, len(lats) - 1)  # lats descend
    dem = dem_ll[lj, li]
    dem = np.nan_to_num(dem, nan=0.0)
    dem = np.maximum(dem, a.sea_level)          # flatten water / bathymetry
    dem_z = dem * a.z_exag
    print(f"  terrain grid {dem.shape}, elev {dem.min():.0f}–{dem.max():.0f} m")

    scale = a.width_mm / (a.x_km * 1000)               # mm per metre
    base_z = dem_z.min() - a.base_mm / scale           # base thickness in real metres
    parts = [terrain_mesh(dem_z, xs, ys, base_z)]

    # 3. buildings
    if not a.no_buildings:
        clip = box(xmin, ymin, xmax, ymax)
        bldgs = fetch_buildings((west, south, east, north), epsg, clip, a.default_levels, a.min_area)
        parts += buildings_mesh(bldgs, dem, xs, ys, a.z_exag)

    # 4. merge, scale, export
    print("  merging & scaling...")
    mesh = trimesh.util.concatenate(parts)
    mesh.apply_translation([-xmin, -ymin, -base_z])
    mesh.apply_scale(scale)
    mesh.export(a.output)
    ext = mesh.bounds[1] - mesh.bounds[0]
    print(f"Wrote {a.output}: {ext[0]:.1f} x {ext[1]:.1f} x {ext[2]:.1f} mm, "
          f"{len(mesh.faces):,} triangles")
    if not mesh.is_watertight:
        print("  note: buildings overlap the terrain rather than being booleaned in — "
              "slicers handle this fine, but if yours complains, run a mesh repair "
              "(e.g. Meshmixer 'Make Solid' or PrusaSlicer 'Fix through Netfabb').")


if __name__ == "__main__":
    sys.exit(main())
