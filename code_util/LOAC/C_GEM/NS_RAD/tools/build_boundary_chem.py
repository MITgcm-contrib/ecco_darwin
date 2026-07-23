"""
Derive riverine (upstream, cub) boundary chemistry from sufficiently-upstream data.

Rule (per project decision): use chemistry from a station far enough upstream to be
freshwater river water (NOT the brackish near-mouth/estuary grabs, which belong to
the SEAWARD/ocean boundary that is already set). Where no adequate upstream data
exists, keep the shipped placeholder.

Sources
  Water Quality Portal grabs (cached CSVs in forcing/: alk*, chem_*,
  ph*), filtered to each river's MAINSTEM station name and to stations well inland
  (freshwater). Colville "R nr Nuiqsut" (~32 km, at the delta head); Sagavanirktok
  "Sagwon"/"Franklin Bluffs" (~55-137 km). Kuparuk keeps its Arctic-LTER headwater
  values (sites/kuparuk.py, tools/lter_boundary.py). Canning has only n=2-3 upstream
  grabs (ambiguous, ~seawater alkalinity) -> placeholder kept.

Unit conversions to model mmol/m^3 (= umol/L):
  alkalinity  mg/L CaCO3   x 19.98   (2 eq per mmol CaCO3, /100.09 g/mol)
  DOC/TOC     mg/L         x 83.26   (/12.011 g C/mol)
  NO3, NH4    mg/L as N    x 71.39   (/14.007 g N/mol)
  PO4         mg/L as PO4  x 10.53   (/94.97 g/mol)
DIC is not measured; it is SOLVED from the observed ALK and pH via the model's
carbonate system (fun_module) -- same method as the Kuparuk/LTER boundary.

Resulting cub values (hard-coded in the site files, with per-river provenance):
  Colville       DIC 1302  ALK 1319  TOC 358  NO3 5.0  NH4 4.0  PO4 0.5  pH 7.80
  Kuparuk        DIC  285  ALK  274  TOC 296  NO3 3.5  NH4 0.4  PO4 0.05 pH 7.32  (LTER)
  Sagavanirktok  DIC 2017  ALK 2078  TOC 250* NO3 11   NH4 1.0  PO4 0.1  pH 8.00  (*TOC lit)
  Canning        placeholder (DIC 1717 ALK 1596 TOC 1582 ...) -- too sparse

Findings worth keeping: the placeholder alkalinity/DIC were NOT wildly wrong for the
carbonate-rich rivers (Colville, Sagavanirktok drain Brooks Range limestone). The
consistent, result-changing error is TOC/DOC -- placeholder 1582 vs actual 250-358,
~5x too high, which overstated respiration and thus CO2 outgassing.

This script re-derives the WQP-based numbers; the literature TOC for the Sag and the
Canning placeholder decision are documented, not computed. Usage: python tools/build_boundary_chem.py
"""
import csv
from pathlib import Path

import numpy as np

FORC = Path(__file__).resolve().parent.parent / "forcing"

# river -> (mainstem name prefix, max latitude to count as sufficiently upstream)
UPSTREAM = {
    "colville": ("COLVILLE R", 70.25),
    "sagavanirktok": ("SAGAVANIRKTOK R", 70.00),
}
SPECIES = {"alk": ("alk", 19.98), "doc": ("chem_doc", 83.26),
           "no3": ("chem_no3", 71.39), "nh4": ("chem_nh4", 71.39),
           "po4": ("chem_po4", 10.53), "ph": ("ph", 1.0)}


def stations(base):
    sf = FORC / (f"{base}_st.csv")
    st = {}
    for r in csv.DictReader(open(sf, encoding="utf-8", errors="replace")):
        try:
            st[r["MonitoringLocationIdentifier"]] = (
                r.get("MonitoringLocationName", "").upper(),
                float(r["LatitudeMeasure"]))
        except (ValueError, KeyError):
            pass
    return st


def median_upstream(base, conv, prefix, maxlat, ph_clip=False):
    st = stations(base)
    rows = list(csv.DictReader(open(FORC / f"{base}.csv", encoding="utf-8", errors="replace")))
    vals = []
    for r in rows:
        sid = r.get("MonitoringLocationIdentifier")
        v = r.get("ResultMeasureValue")
        if sid not in st or v in ("", "nan"):
            continue
        nm, lat = st[sid]
        if nm.startswith(prefix) and lat < maxlat:
            try:
                x = float(v)
                if ph_clip and not (5 < x < 9.5):
                    continue
                vals.append(x)
            except ValueError:
                pass
    return (np.median(vals) * conv, len(vals)) if vals else (None, 0)


def main():
    print(f"{'river':15s} {'ALK':>10s} {'TOC':>10s} {'NO3':>9s} {'NH4':>9s} {'PO4':>9s} {'pH':>7s}")
    for riv, (prefix, maxlat) in UPSTREAM.items():
        out = {}
        for tag, (base, conv) in SPECIES.items():
            out[tag] = median_upstream(base, conv, prefix, maxlat, ph_clip=(tag == "ph"))
        def s(t):
            return f"{out[t][0]:6.0f}(n{out[t][1]})" if out[t][0] is not None else "   --   "
        print(f"{riv:15s} {s('alk'):>10s} {s('doc'):>10s} {s('no3'):>9s} "
              f"{s('nh4'):>9s} {s('po4'):>9s} {out['ph'][0]:7.2f}")
    print("\nDIC solved from ALK+pH in the site files (fun_module carbonate system).")
    print("Kuparuk = LTER headwater; Canning = placeholder (too sparse).")


if __name__ == "__main__":
    main()
