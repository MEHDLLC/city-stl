# City STL

Generate a 3D-printable STL of any city — terrain plus buildings — with no local install.

## Run it (GitHub Actions)

1. Go to **Actions → Build city STL → Run workflow**.
2. Enter a city, size in km, and printed width in mm.
3. When the run finishes, download the STL from the **Artifacts** section on the run page.

## Run it locally

```
pip install -r requirements.txt
python city_stl.py "Seattle, WA" --x-km 4 --y-km 3 --width-mm 150 -o seattle.stl
```

Data: OpenStreetMap buildings (via Overpass), AWS Terrain Tiles for elevation.
