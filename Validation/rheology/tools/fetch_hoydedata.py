#!/usr/bin/env python3
# Architect: Hans Bihs
"""Download a terrain model (DTM) of one Kartverket laser project from hoydedata.no as a GeoTIFF.

  fetch_hoydedata.py <LAS project id> <out.tif> --bbox xmin ymin xmax ymax [--res 5] [--epsg 25832]
  fetch_hoydedata.py --list --bbox xmin ymin xmax ymax [--epsg 25832]

Uses the public ArcGIS ImageServer Prosjekt_DTM of hoydedata.no (exportImage with the project as the mosaic rule,
LAS_PROJECT_ID; --list prints the projects that cover the box with their year and flight date). Bilinear
resampling of the 0.25 - 1 m project DTM onto the cells of --res; no data -9999. Tiles of at most 4000 x 4000
pixels are merged. Needs only the Python standard library (and numpy for more than one tile).

Gjerdrum 2020 (EPSG:25832): pre-event 4241 "Gjerdrum Ullensaker Nannestad 5pkt 2020" (flown 19.04.2020),
post-event 5001 "Gjerdrum 5pkt 2021" (08.01.2021):
  python3 fetch_hoydedata.py 4241 gjerdrum_pre_2020_DTM_5m.tif  --bbox 611000 6658000 615500 6662500 --res 5
  python3 fetch_hoydedata.py 5001 gjerdrum_post_2021_DTM_5m.tif --bbox 611000 6658000 615500 6662500 --res 5
"""
import sys, json, argparse, struct, urllib.request, urllib.parse
from datetime import datetime, timezone

BASE = "https://hoydedata.no/arcgis/rest/services/Prosjekt_DTM/ImageServer"


def get(url, params):
    with urllib.request.urlopen(url + "?" + urllib.parse.urlencode(params), timeout=300) as r:
        return r.read()


def list_projects(bbox, epsg):
    g = json.dumps({"xmin": bbox[0], "ymin": bbox[1], "xmax": bbox[2], "ymax": bbox[3],
                    "spatialReference": {"wkid": epsg}})
    q = json.loads(get(BASE + "/query", {"f": "json", "geometry": g, "geometryType": "esriGeometryEnvelope",
                                         "inSR": epsg, "spatialRel": "esriSpatialRelIntersects",
                                         "where": "CATEGORY=1", "returnGeometry": "false",
                                         "outFields": "LAS_PROJECT_ID,LAS_PROJECT_NAME,AARSTALL,SISTEFLYDATO,OPPLOSNING,PUNKTTETTHET"}))
    seen = {}
    for f in q.get("features", []):
        a = f["attributes"]
        seen[a["LAS_PROJECT_ID"]] = a
    for pid, a in sorted(seen.items(), key=lambda kv: kv[1]["SISTEFLYDATO"] or 0):
        d = datetime.fromtimestamp((a["SISTEFLYDATO"] or 0) / 1000, timezone.utc).strftime("%Y-%m-%d")
        print("%6d  %s  flown %s  %.2f m  %s pts/m2" % (pid, a["LAS_PROJECT_NAME"], d, a["OPPLOSNING"] or 0,
                                                       a["PUNKTTETTHET"]))


def tile(pid, x0, y0, x1, y1, nx, ny, epsg):
    """float32 rows north to south"""
    b = get(BASE + "/exportImage", {"bbox": "%f,%f,%f,%f" % (x0, y0, x1, y1), "bboxSR": epsg, "imageSR": epsg,
                                    "size": "%d,%d" % (nx, ny), "format": "bsq", "pixelType": "F32",
                                    "noData": -9999, "interpolation": "RSP_BilinearInterpolation",
                                    "mosaicRule": json.dumps({"mosaicMethod": "esriMosaicNone",
                                                              "where": "LAS_PROJECT_ID=%d" % pid}),
                                    "f": "image"})
    n = nx * ny
    return list(struct.unpack("<%df" % n, b[:4 * n]))


def write_geotiff(path, z, nx, ny, x0, ymax, res, epsg):
    """minimal uncompressed single-strip float32 GeoTIFF (ModelPixelScale, ModelTiepoint, GeoKeys, GDAL_NODATA)"""
    data = struct.pack("<%df" % (nx * ny), *z)
    entries = []
    extra = b""
    off0 = 8

    def add(tag, typ, count, value_bytes):
        entries.append((tag, typ, count, value_bytes))
    nodata = b"-9999\x00"
    geokeys = [1, 1, 0, 3, 1024, 0, 1, 1, 1025, 0, 1, 1, 3072, 0, 1, epsg]
    tags = [(256, 4, 1, struct.pack("<I", nx)), (257, 4, 1, struct.pack("<I", ny)),
            (258, 3, 1, struct.pack("<H", 32)), (259, 3, 1, struct.pack("<H", 1)),
            (262, 3, 1, struct.pack("<H", 1)), (273, 4, 1, None), (277, 3, 1, struct.pack("<H", 1)),
            (278, 4, 1, struct.pack("<I", ny)), (279, 4, 1, struct.pack("<I", len(data))),
            (284, 3, 1, struct.pack("<H", 1)), (339, 3, 1, struct.pack("<H", 3)),
            (33550, 12, 3, struct.pack("<3d", res, res, 0.0)),
            (33922, 12, 6, struct.pack("<6d", 0, 0, 0, x0, ymax, 0.0)),
            (34735, 3, len(geokeys), struct.pack("<%dH" % len(geokeys), *geokeys)),
            (42113, 2, len(nodata), nodata)]
    nent = len(tags)
    ifd_size = 2 + 12 * nent + 4
    extra_off = off0 + ifd_size
    blobs = b""
    ent = b""
    data_off = None
    for tag, typ, cnt, vb in tags:
        if vb is None:
            ent += struct.pack("<HHI", tag, typ, cnt) + b"\x00\x00\x00\x00"
            continue
        if len(vb) <= 4:
            ent += struct.pack("<HHI", tag, typ, cnt) + vb.ljust(4, b"\x00")
        else:
            ent += struct.pack("<HHII", tag, typ, cnt, extra_off + len(blobs))
            blobs += vb + (b"\x00" if len(vb) % 2 else b"")
    data_off = extra_off + len(blobs)
    # patch StripOffsets
    i = [t[0] for t in tags].index(273)
    ent = ent[:12 * i] + struct.pack("<HHII", 273, 4, 1, data_off) + ent[12 * (i + 1):]
    with open(path, "wb") as f:
        f.write(b"II*\x00" + struct.pack("<I", off0))
        f.write(struct.pack("<H", nent) + ent + struct.pack("<I", 0))
        f.write(blobs)
        f.write(data)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("project", nargs="?", type=int)
    ap.add_argument("out", nargs="?")
    ap.add_argument("--bbox", nargs=4, type=float, required=True)
    ap.add_argument("--res", type=float, default=5.0)
    ap.add_argument("--epsg", type=int, default=25832)
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    if a.list or a.project is None:
        list_projects(a.bbox, a.epsg)
        return
    x0, y0, x1, y1 = a.bbox
    nx, ny = int(round((x1 - x0) / a.res)), int(round((y1 - y0) / a.res))
    T = 4000
    z = [0.0] * (nx * ny)
    for jt in range(0, ny, T):              # jt from the north
        for it in range(0, nx, T):
            tx, ty = min(T, nx - it), min(T, ny - jt)
            v = tile(a.project, x0 + it * a.res, y1 - (jt + ty) * a.res, x0 + (it + tx) * a.res, y1 - jt * a.res,
                     tx, ty, a.epsg)
            for r in range(ty):
                z[(jt + r) * nx + it:(jt + r) * nx + it + tx] = v[r * tx:(r + 1) * tx]
            print("tile %d %d: %d x %d" % (it, jt, tx, ty))
    write_geotiff(a.out, z, nx, ny, x0, y1, a.res, a.epsg)
    ok = [v for v in z if v > -9000]
    print("%s: %d x %d, %.0f %% with data, %.2f - %.2f m" % (a.out, nx, ny, 100.0 * len(ok) / len(z),
                                                          min(ok) if ok else 0, max(ok) if ok else 0))


if __name__ == "__main__":
    main()
