# City STL

Generate a 3D-printable STL of any city — terrain plus buildings — with no local install.

## Run it (GitHub Actions)

1. **Actions → Build city STL → Run workflow**.
2. Location can be a place name (`Chelsea, Manhattan, New York`), a US ZIP (`10001`),
   or coordinates (`40.7484, -73.9857`). Enter size in km and printed width in mm.
3. Download the STL from **Artifacts** on the finished run.

## Run it locally

```
pip install -r requirements.txt
python city_stl.py "Seattle, WA" --x-km 4 --y-km 3 --width-mm 150 -o seattle.stl
```

## Previews

Every run also produces PNG renders — one per tile plus, for tiled runs, an
`_assembled.png` showing all tiles together with the seams marked — so you can
check the result before downloading or slicing anything.

## Tiling (big prints)

Set **tiles** to e.g. `3x3` and the model is written as separate STLs that butt
together edge-to-edge (`name_r1_c1.stl` … r1 is the south row, c1 the west column;
a `_layout.txt` shows the grid). Printed width is the whole assembled model, so
10 km at 600 mm wide in 3x3 tiles gives nine 200 mm pieces at 1:16,700 — three
times the height detail of a single 200 mm print. Buildings and landmarks that
cross a tile edge are cut cleanly so each half prints on its own tile.

## Building detail

Buildings come from OpenStreetMap. Where mapped, the script uses `building:part`
pieces (setbacks, spires, towers) and `roof:shape` tags — pointed roofs for
pyramidal/hipped/gabled, rounded for dome/onion — so well-mapped landmarks look
like themselves rather than plain boxes.

## Scouting: which landmarks need a model?

Before a big build, run the workflow with **scout** ticked. Instead of building, it
lists everything in the area that won't look right as an extruded box — statues,
monuments, towers, cathedrals, stadiums, and the tallest named buildings — with
coordinates and heights pre-filled, as `NAME_landmark_suggestions.json` /
`.txt`. Copy the entries you want into `landmarks/landmarks.json`, upload matching
STLs, then run the real build.

## Landmarks (real 3D models)

Statues and sculptural buildings can't be reconstructed from footprints. Drop a
printable STL into `landmarks/` and describe it in `landmarks/landmarks.json`:

| field          | meaning                                                                 |
| -------------- | ----------------------------------------------------------------------- |
| `name`         | label for the log                                                       |
| `lat`, `lon`   | where to place it                                                       |
| `file`         | STL/OBJ filename inside `landmarks/`                                    |
| `height_m`     | real-world height; model is scaled uniformly to this                    |
| `radius_m`     | OSM footprints within this distance are removed (default 60)            |
| `match_name`   | also remove any OSM building whose name contains this text (optional)   |
| `rotation_deg` | spin about vertical axis (optional)                                     |
| `rotate_x_deg` | set to 90 for Y-up models that import lying down (optional)             |

Entries whose file is missing are skipped, so you can keep the list long and only
add models as you find them (Printables / Thingiverse — check the licence).

Data: OpenStreetMap buildings (via Overpass), AWS Terrain Tiles for elevation.
