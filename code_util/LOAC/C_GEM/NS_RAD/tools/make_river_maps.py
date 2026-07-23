"""
Plan-view river-network maps: SWORD channel network overlaid on satellite imagery,
one panel per river, into docs/ns_rad_river_networks.pdf.

Imagery: Sentinel-2 cloudless (EOX, 10 m, cloud-free true-color mosaic). NOTE the user
asked for Landsat, but NASA GIBS Landsat (WELD) and USGS LandsatLook do NOT cover the
North Slope (70 N) -- they return blank. Sentinel-2 cloudless is the best available
cloud-free true-color basemap for the Arctic and is higher resolution than Landsat's
30 m, so it is used instead (flagged on the figure).

SWORD nodes (v17b) are coloured by channel width, so the wide seaward delta and the
narrow inland channels read directly; this is the same data behind the delta-mouth
widths (tools/extract_sword.py). The model's seaward boundary width is annotated.

Usage:  python tools/make_river_maps.py         (needs /tmp/sword/netcdf/na_sword_v17b.nc
                                                 and network access for the imagery)
"""
import io
import json
import urllib.request
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nsrad_style as _S  # shared NS-RAD publication/presentation style
_S.apply()
_S.install_autoscale(1.2)  # embed fonts + presentation-size bump
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import LineCollection  # noqa: F401  (kept for future reach lines)
from PIL import Image
from netCDF4 import Dataset

ROOT = Path(__file__).resolve().parent.parent
SWORD_PATH = Path("/tmp/sword/netcdf/na_sword_v17b.nc")
OUT = ROOT / "docs" / "ns_rad_river_networks.pdf"
SWJSON = ROOT / "docs" / "sword_widths.json"

LABEL = {"colville": "Colville", "kuparuk": "Kuparuk",
         "sagavanirktok": "Sagavanirktok", "canning": "Canning"}
# delta-region bounding boxes (lon0, lon1, lat0, lat1)
BBOX = {"colville": (-151.5, -150.0, 70.10, 70.60),
        "kuparuk": (-149.35, -148.45, 70.12, 70.52),
        "sagavanirktok": (-148.45, -147.65, 70.02, 70.42),
        "canning": (-146.35, -145.35, 69.92, 70.28)}
INK, MUTED = "#0b0b0b", "#8a8880"

WMS = ("https://tiles.maps.eox.at/wms?SERVICE=WMS&REQUEST=GetMap&VERSION=1.3.0&"
       "LAYERS=s2cloudless-2020&CRS=EPSG:4326&BBOX={lat0},{lon0},{lat1},{lon1}&"
       "WIDTH={w}&HEIGHT={h}&FORMAT=image/png")


def fetch_imagery(lon0, lon1, lat0, lat1):
    latm = 0.5 * (lat0 + lat1)
    km_x = (lon1 - lon0) * 111.0 * np.cos(np.radians(latm))
    km_y = (lat1 - lat0) * 111.0
    w = 900
    h = int(w * km_y / km_x)
    url = WMS.format(lat0=lat0, lon0=lon0, lat1=lat1, lon1=lon1, w=w, h=h)
    with urllib.request.urlopen(url, timeout=60) as r:
        img = Image.open(io.BytesIO(r.read())).convert("RGB")
    return np.asarray(img)


def load_nodes():
    d = Dataset(SWORD_PATH); n = d.groups["nodes"]
    return dict(x=n.variables["x"][:], y=n.variables["y"][:],
               w=n.variables["width"][:], lake=n.variables["lakeflag"][:],
               name=np.array([str(s) for s in n.variables["river_name"][:]]))


def main():
    OUT.parent.mkdir(exist_ok=True)
    nodes = load_nodes()
    sw = json.load(open(SWJSON)) if SWJSON.exists() else {}
    with PdfPages(OUT) as pdf:
        for site in ["colville", "kuparuk", "sagavanirktok", "canning"]:
            lon0, lon1, lat0, lat1 = BBOX[site]
            try:
                img = fetch_imagery(lon0, lon1, lat0, lat1)
            except Exception as e:
                print(f"  {site}: imagery fetch failed ({e}); skipping")
                continue
            latm = 0.5 * (lat0 + lat1)
            fig, ax = plt.subplots(figsize=(9, 9 * img.shape[0] / img.shape[1] + 0.6))
            ax.imshow(img, extent=[lon0, lon1, lat0, lat1], origin="upper")
            ax.set_aspect(1.0 / np.cos(np.radians(latm)))  # equirectangular at this latitude
            # SWORD channel network in the box (all channels, coloured by width)
            m = ((nodes["x"] > lon0) & (nodes["x"] < lon1) & (nodes["y"] > lat0)
                 & (nodes["y"] < lat1) & (nodes["lake"] == 0) & (nodes["w"] > 0))
            sc = ax.scatter(nodes["x"][m], nodes["y"][m], c=nodes["w"][m], s=7,
                            cmap="cool", vmin=0, vmax=600, edgecolor="none", alpha=0.9)
            cb = fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.02)
            cb.set_label("SWORD channel width (m)", fontsize=8)
            info = sw.get(site, {})
            dsum = info.get("delta_sum_m")
            ttl = f"{LABEL[site]} — SWORD channel network on Sentinel-2 cloudless"
            ax.set_title(ttl, loc="left", fontsize=12, weight="bold", color=INK)
            note = (f"model seaward width B_lb = {dsum} m (SWORD distributary sum"
                    f" of {info.get('chan_widths')})" if dsum else "")
            ax.set_xlabel(f"longitude   ·   {note}", fontsize=8, color=MUTED)
            ax.set_ylabel("latitude", fontsize=8, color=MUTED)
            ax.tick_params(labelsize=7)
            fig.text(0.012, 0.012, "Imagery: Sentinel-2 cloudless (EOX). Landsat (GIBS WELD / "
                     "USGS LandsatLook) does not cover 70°N. Network: SWORD v17b nodes.",
                     fontsize=6.5, color=MUTED)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
            print(f"  {site}: map done ({img.shape[1]}x{img.shape[0]} px, "
                  f"{m.sum()} SWORD nodes)")
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
