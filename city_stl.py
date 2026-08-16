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
    --landmarks FILE    JSON list of real 3D models to swap in for landmarks (see README)
    --tiles 2x2         split into COLSxROWS separate STLs that butt together (edges shared)

Buildings use OSM building:part pieces (setbacks, spires, towers) and roof:shape tags
(pyramidal/hipped/gabled -> pointed, dome/onion -> domed) when mapped.

Data: OpenStreetMap (buildings, via osmnx) and AWS Terrain Tiles (elevation, no key needed).
Install: pip install osmnx trimesh shapely numpy scipy pillow requests pyproj mapbox-earcut
"""

import argparse
import json
import math
import os
import re
import sys
import numpy as np
import pandas as pd
import requests
from io import BytesIO
from PIL import Image

import osmnx as ox
import trimesh
from shapely.geometry import Polygon, MultiPolygon, Point, box
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
def _num(v):
    """Parse an OSM numeric tag ('12', '12 m', '40 ft', 12.0). Returns float or None."""
    if v is None:
        return None
    if isinstance(v, (int, float)):
        return None if np.isnan(v) else float(v)
    t = str(v).strip().lower().replace(",", ".")
    try:
        if t.endswith("ft") or t.endswith("'"):
            return float(t.rstrip("ft'").strip()) * 0.3048
        return float(t.rstrip("m").strip())
    except ValueError:
        return None


ROOF_POINTED = {"pyramidal", "gabled", "hipped", "half-hipped", "mansard", "cone",
                "gambrel", "round", "saltbox"}
ROOF_DOMED = {"dome", "onion"}


def building_geometry(row, default_levels):
    """From OSM tags -> (min_h, wall_top_h, roof_shape, roof_h), all metres above ground."""
    total = _num(row.get("height")) or _num(row.get("building:height"))
    if total is None:
        lv = _num(row.get("building:levels"))
        total = (lv if lv is not None else default_levels) * 3.2 + 1.0
    min_h = _num(row.get("min_height")) or 0.0
    if min_h >= total:
        min_h = 0.0
    shape = str(row.get("roof:shape") or "flat").lower()
    roof_h = _num(row.get("roof:height"))
    if roof_h is None:
        rl = _num(row.get("roof:levels"))
        roof_h = rl * 3.0 if rl else None
    if shape in ROOF_POINTED or shape in ROOF_DOMED:
        if roof_h is None:
            roof_h = min(max(0.35 * (total - min_h), 2.0), 0.6 * (total - min_h))
        roof_h = min(roof_h, 0.9 * (total - min_h))
        return min_h, total - roof_h, shape, roof_h
    return min_h, total, "flat", 0.0


def fetch_buildings(bbox_wgs, epsg, clip_box, default_levels, min_area):
    """Returns list of (polygon, min_h, wall_top_h, roof_shape, roof_h, name)."""
    west, south, east, north = bbox_wgs
    print("  fetching OSM buildings & building parts...")
    try:
        gdf = ox.features_from_bbox((west, south, east, north),
                                    tags={"building": True, "building:part": True})
    except Exception as e:
        print(f"  no buildings returned ({e})")
        return []
    gdf = gdf[gdf.geometry.type.isin(["Polygon", "MultiPolygon"])].to_crs(epsg=epsg)

    def col(name):
        return gdf[name] if name in gdf.columns else pd.Series([None] * len(gdf), index=gdf.index)

    is_part = col("building:part").notna() & (col("building:part").astype(str) != "no")
    is_bldg = col("building").notna() & (col("building").astype(str) != "no") & ~is_part
    parts, outlines = gdf[is_part], gdf[is_bldg]

    # drop outlines that are (mostly) covered by parts — the parts carry the 3D detail
    if len(parts):
        parts_union = parts.geometry.union_all() if hasattr(parts.geometry, "union_all") \
            else parts.geometry.unary_union
        covered = outlines.geometry.intersection(parts_union).area / outlines.geometry.area
        outlines = outlines[covered.fillna(0) < 0.5]
    keep = pd.concat([outlines, parts])
    print(f"  {len(outlines)} outlines + {len(parts)} parts")

    geoms = keep.geometry.intersection(clip_box)
    tag_cols = ["height", "building:height", "building:levels", "min_height",
                "roof:shape", "roof:height", "roof:levels", "name"]
    cols = {c: keep[c].tolist() for c in tag_cols if c in keep.columns}
    out = []
    for k in range(len(keep)):
        geom = geoms.iloc[k]
        if geom.is_empty:
            continue
        row = {c: v[k] for c, v in cols.items()}
        min_h, top_h, shape, roof_h = building_geometry(row, default_levels)
        name = row.get("name")
        name = str(name) if isinstance(name, str) else None
        polys = list(geom.geoms) if isinstance(geom, MultiPolygon) else [geom]
        for p in polys:
            if isinstance(p, Polygon) and p.area >= min_area:
                out.append((p, min_h, top_h, shape, roof_h, name))
    print(f"  {len(out)} footprints kept")
    return out


def sample_dem(dem, xs, ys, px, py):
    j = np.clip(np.searchsorted(xs, px), 0, len(xs) - 1)
    i = np.clip(np.searchsorted(ys, py), 0, len(ys) - 1)
    return dem[i, j]


def _side_faces(upper, lower):
    """Quads between two rings of equal length (index arrays), consistent winding."""
    m = len(upper)
    i, j = np.arange(m), (np.arange(m) + 1) % m
    return np.vstack([np.column_stack([upper[i], upper[j], lower[j]]),
                      np.column_stack([upper[i], lower[j], lower[i]])])


def _extrude(poly, z0, z1, roof="flat", roof_h=0.0):
    """Prism (optionally with pointed/domed roof) from a shapely polygon. -> (verts, faces)."""
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
    V = [top, bot]
    F = [tri[:, ::-1] + n]                              # floor (facing down)
    off = 0
    for r in rings:                                     # walls
        m = len(r)
        idx = np.arange(m) + off
        F.append(_side_faces(idx, idx + n))
        off += m
    nv = 2 * n
    shaped = roof in ROOF_POINTED or roof in ROOF_DOMED
    if shaped and roof_h > 0 and len(rings) == 1:
        ext = rings[0]
        c = np.asarray(poly.centroid.coords)[0]
        m = len(ext)
        cur = np.arange(m)                              # indices of top ring
        if roof in ROOF_DOMED:
            K = 5
            for k in range(1, K):
                t = k * (math.pi / 2) / K
                ring = c + (ext - c) * math.cos(t)
                V.append(np.column_stack([ring, np.full(m, z1 + roof_h * math.sin(t))]))
                nxt = np.arange(m) + nv
                F.append(_side_faces(nxt, cur))
                cur, nv = nxt, nv + m
        V.append(np.array([[c[0], c[1], z1 + roof_h]]))
        apex = nv
        nv += 1
        i, j = cur, np.roll(cur, -1)
        F.append(np.column_stack([j, i, np.full(m, apex)]))
    else:
        F.append(tri)                                   # flat roof
    return np.vstack(V), np.vstack(F)


def buildings_mesh(bldgs, dem, xs, ys, z_exag):
    """Batch all buildings into one mesh with numpy."""
    V, F, base = [], [], 0
    for poly, min_h, top_h, shape, roof_h, _name in bldgs:
        c = poly.centroid
        gz = sample_dem(dem, xs, ys, c.x, c.y) * z_exag
        z0 = gz + min_h - (0.0 if min_h > 0 else 2.0)   # sink ground-floor parts 2 m
        r = _extrude(poly.simplify(0.3), z0, gz + top_h, shape, roof_h)
        if r is None:
            continue
        v, f = r
        V.append(v)
        F.append(f + base)
        base += len(v)
    if not V:
        return []
    return [trimesh.Trimesh(vertices=np.vstack(V), faces=np.vstack(F), process=False)]


# ----------------------------------------------------------------------------- landmarks
def load_landmarks(path):
    """landmarks.json: list of {name, lat, lon, file, height_m, radius_m, rotation_deg, rotate_x_deg}"""
    if not path or not os.path.exists(path):
        return []
    with open(path) as fh:
        items = json.load(fh)
    base = os.path.dirname(os.path.abspath(path))
    for it in items:
        it["file"] = os.path.join(base, it["file"]) if not os.path.isabs(it["file"]) else it["file"]
    return items


def apply_landmarks(items, bldgs, to_utm, xmin, xmax, ymin, ymax, dem, xs, ys, z_exag):
    """Remove footprints near each landmark and return placed landmark meshes."""
    meshes = []
    for it in items:
        x, y = to_utm.transform(it["lon"], it["lat"])
        if not (xmin <= x <= xmax and ymin <= y <= ymax):
            continue
        if not os.path.exists(it["file"]):
            print(f"  landmark '{it.get('name')}' skipped: {it['file']} not found")
            continue
        radius = float(it.get("radius_m", 60))
        nm = (it.get("match_name") or "").lower()
        before = len(bldgs)
        bldgs[:] = [b for b in bldgs
                    if not (b[0].centroid.distance(Point(x, y)) <= radius or
                            (nm and b[5] and nm in b[5].lower()))]
        m = trimesh.load(it["file"], force="mesh")
        if float(it.get("rotate_x_deg", 0)):
            m.apply_transform(trimesh.transformations.rotation_matrix(
                math.radians(float(it["rotate_x_deg"])), [1, 0, 0]))
        ext = m.bounds[1] - m.bounds[0]
        m.apply_scale(float(it["height_m"]) / ext[2])
        if float(it.get("rotation_deg", 0)):
            m.apply_transform(trimesh.transformations.rotation_matrix(
                math.radians(float(it["rotation_deg"])), [0, 0, 1]))
        lo, hi = m.bounds
        gz = sample_dem(dem, xs, ys, x, y) * z_exag
        m.apply_translation([x - (lo[0] + hi[0]) / 2, y - (lo[1] + hi[1]) / 2, gz - lo[2] - 1.0])
        meshes.append(m)
        print(f"  landmark '{it.get('name')}' placed ({before - len(bldgs)} footprints replaced)")
    return meshes


def clip_mesh_to_box(mesh, x0, x1, y0, y1):
    """Return the part of a mesh inside an XY box (None if nothing). Used to split
    landmark models across tile boundaries."""
    lo, hi = mesh.bounds
    if hi[0] <= x0 or lo[0] >= x1 or hi[1] <= y0 or lo[1] >= y1:
        return None
    if lo[0] >= x0 and hi[0] <= x1 and lo[1] >= y0 and hi[1] <= y1:
        return mesh
    cx, cy = (lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2
    m = mesh.copy()
    m.apply_translation([-cx, -cy, 0])          # slice in local coords (numerically safer)
    try:
        for normal, origin in (([1, 0, 0], [x0 - cx, 0, 0]), ([-1, 0, 0], [x1 - cx, 0, 0]),
                               ([0, 1, 0], [0, y0 - cy, 0]), ([0, -1, 0], [0, y1 - cy, 0])):
            m = trimesh.intersections.slice_mesh_plane(m, normal, origin, cap=True)
            if m is None or len(m.faces) == 0:
                return None
        m.apply_translation([cx, cy, 0])
        return m
    except Exception as e:
        print(f"  (landmark clip failed: {e}; placing whole model in one tile)")
        return mesh if (x0 <= cx < x1 and y0 <= cy < y1) else None


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
    ap.add_argument("--landmarks", default="landmarks/landmarks.json",
                    help="JSON list of landmark models to swap in for OSM footprints")
    ap.add_argument("--tiles", default="1x1",
                    help="split into COLSxROWS tiles, e.g. 3x2 (tiles share edges, print separately)")
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
    bldgs, landmark_meshes = [], []

    # 3. buildings & landmarks (fetched once for the whole area)
    if not a.no_buildings:
        clip = box(xmin, ymin, xmax, ymax)
        bldgs = fetch_buildings((west, south, east, north), epsg, clip, a.default_levels, a.min_area)
        landmark_meshes = apply_landmarks(load_landmarks(a.landmarks), bldgs, to_utm,
                                          xmin, xmax, ymin, ymax, dem, xs, ys, a.z_exag)

    # 4. split into tiles, scale, export
    try:
        cols, rows = (int(v) for v in a.tiles.lower().split("x"))
    except ValueError:
        sys.exit("--tiles must look like 2x2")
    cols, rows = max(cols, 1), max(rows, 1)
    xi = np.linspace(0, len(xs) - 1, cols + 1).round().astype(int)   # tile edges on grid lines
    yi = np.linspace(0, len(ys) - 1, rows + 1).round().astype(int)
    stem, ext = os.path.splitext(a.output)
    ext = ext or ".stl"
    written = []
    for r in range(rows):
        for c_ in range(cols):
            i0, i1, j0, j1 = xi[c_], xi[c_ + 1], yi[r], yi[r + 1]
            tx0, tx1, ty0, ty1 = xs[i0], xs[i1], ys[j0], ys[j1]
            tile_box = box(tx0, ty0, tx1, ty1)
            tparts = [terrain_mesh(dem_z[j0:j1 + 1, i0:i1 + 1], xs[i0:i1 + 1], ys[j0:j1 + 1], base_z)]
            if not a.no_buildings:
                tb = []
                for b in bldgs:
                    g = b[0].intersection(tile_box)
                    if g.is_empty:
                        continue
                    for pp in (g.geoms if isinstance(g, MultiPolygon) else [g]):
                        if isinstance(pp, Polygon) and pp.area > 1.0:
                            tb.append((pp,) + b[1:])
                tparts += buildings_mesh(tb, dem, xs, ys, a.z_exag)
                for lm in landmark_meshes:
                    piece = clip_mesh_to_box(lm, tx0, tx1, ty0, ty1)
                    if piece is not None:
                        tparts.append(piece)
            mesh = trimesh.util.concatenate(tparts)
            mesh.apply_translation([-tx0, -ty0, -base_z])
            mesh.apply_scale(scale)
            name = a.output if cols * rows == 1 else f"{stem}_r{r + 1}_c{c_ + 1}{ext}"
            mesh.export(name)
            e = mesh.bounds[1] - mesh.bounds[0]
            print(f"Wrote {name}: {e[0]:.1f} x {e[1]:.1f} x {e[2]:.1f} mm, {len(mesh.faces):,} triangles")
            written.append(name)
    if cols * rows > 1:
        print(f"\n{cols}x{rows} tiles — r1 is the SOUTH row, c1 the WEST column; "
              f"whole model {a.width_mm:.0f} x {a.width_mm * a.y_km / a.x_km:.0f} mm.")
        with open(f"{stem}_layout.txt", "w") as fh:
            fh.write("Tile layout (north at top):\n")
            for r in range(rows, 0, -1):
                fh.write("  " + "  ".join(f"r{r}_c{c_}" for c_ in range(1, cols + 1)) + "\n")
    print("  note: buildings overlap the terrain rather than being booleaned in — "
          "slicers handle this fine; run a mesh repair (Meshmixer/PrusaSlicer) if yours complains.")


if __name__ == "__main__":
    sys.exit(main())
