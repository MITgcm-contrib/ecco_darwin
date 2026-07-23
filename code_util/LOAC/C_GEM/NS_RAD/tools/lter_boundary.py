"""
Derive the Kuparuk riverine boundary chemistry from the Arctic LTER Streams record.

Source: EDI knb-lter-arc.10303 (Arctic LTER Streams Chemistry, Toolik Field Station,
1978-2019). NATURAL "Reference" reach only -- the other reaches are a long-term
phosphorus-fertilization EXPERIMENT and are excluded.

Spatial caveat (the important one): the LTER Kuparuk site is at the Dalton Highway
crossing, 68.6 N / 720 m elevation -- ~190 km straight-line (~266 river-km) from the
mouth, i.e. ~163 km upstream of this model's own upstream boundary. These are
HEADWATER concentrations, a first-order proxy for the delta boundary, not a
measurement of it. Alkalinity in particular (soft tundra headwater) is likely a
lower bound for the delta.

Units: LTER uM = umol/L = mmol/m^3 = the model's concentration units, so no
conversion for NO3/NH4/PO4/TOC. Alkalinity ueq/L ~ mmol/m^3. DIC is not measured;
it is scaled from ALK to keep a freshwater DIC:ALK ratio.

Only the Kuparuk has LTER data; the other three rivers keep the shared placeholder.
The values printed here are what is hard-coded into sites/kuparuk.py.

Usage:  python tools/lter_boundary.py     (re-derives from a cached CSV or downloads)
"""
import csv
import io
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "forcing" / "lter_streams_chem.csv"
ENTITY = ("https://pasta.lternet.edu/package/data/eml/knb-lter-arc/10303/7/"
          "13c2451ec91b6d5151190124dcf32a5b")

# LTER 'Type' -> (model species, old placeholder cub)
MAP = [
    ("Nitrate+Nitrite", "NO3", 7.68),
    ("Ammonium", "NH4", 3.9),
    ("Soluble Reactive Phosphorus", "PO4", 0.01),
    ("Dissolved Organic Carbon", "TOC", 1582),
    ("pH", "pH", 7.5),
    ("Alkalinity", "ALK", 1596),
]
DIC_OVER_ALK = 1717 / 1596   # placeholder ratio, to scale the unmeasured DIC


def load_rows():
    if not CACHE.exists():
        print(f"downloading LTER streams chemistry -> {CACHE.name}")
        urllib.request.urlretrieve(ENTITY, CACHE)
    text = CACHE.read_text(encoding="utf-8", errors="replace")
    return list(csv.DictReader(io.StringIO(text)))


def main():
    rows = [r for r in load_rows()
            if "upar" in r.get("River", "") and "Reference" in r.get("Reach", "")]
    print(f"Kuparuk reference-reach rows: {len(rows)}  (fertilized reaches excluded)\n")
    print(f"{'species':6s} {'LTER median':>12s}  {'placeholder':>11s}  note")
    out = {}
    for lter, tag, plc in MAP:
        v = []
        for r in rows:
            if r["Type"] == lter:
                try:
                    v.append(float(r["Value"]))
                except ValueError:
                    pass
        if not v:
            continue
        med = float(np.median(v))
        out[tag] = med
        print(f"{tag:6s} {med:12.2f}  {plc:11.1f}  n={len(v)}")
    dic = out["ALK"] * DIC_OVER_ALK
    out["DIC"] = dic
    print(f"{'DIC':6s} {dic:12.2f}  {1717:11.1f}  scaled from ALK (LTER has no DIC)")
    print("\nThese are the cub values hard-coded in sites/kuparuk.py "
          "(marine clb kept from _baseline).")


if __name__ == "__main__":
    main()
