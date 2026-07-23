"""
Build per-river tidal-constituent forcing from NOAA CO-OPS harmonic analyses.

Writes forcing/tidal_constituents.json. config.py loads the active
site's constituents and fun_module.Tide() reconstructs the sea-surface elevation as

    eta(t) = sum_i  A_i * cos( speed_i * (t/3600) - G_i )

with A_i in metres, speed_i in deg/hr, G_i the GMT phase lag in degrees. This
replaces the shipped single-sinusoid tide (AMPL=0.2 m, pfun=0.0787 cyc/hr).

Each river is matched to its nearest CO-OPS harmonic station on the Beaufort coast:

    Colville       -> 9496345  Cape Halkett
    Kuparuk        -> 9497645  Prudhoe Bay
    Sagavanirktok  -> 9497778  Cross Island (Dinkum Sands)
    Canning        -> 9498306  Point Thomson

Only SHORT-PERIOD tidal constituents are kept. The station analyses also report
SA and SSA (solar annual / semi-annual) with large amplitude (~0.38 ft) -- these are
the SEASONAL sea-level cycle, not tides, and are excluded so they do not conflate
with the model's discharge-driven seasonal signal.

Reality check: the Beaufort coast is microtidal. M2 amplitude is only ~0.06-0.07 m
and the full short-period range is ~0.4 m, so this is a modest, realistic refinement
of the old 0.2 m idealisation -- it adds spring-neap and diurnal-inequality structure
and per-river phase, not a large amplitude change.

Usage:  python tools/build_tides.py
"""
import json
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "forcing" / "tidal_constituents.json"
FT_TO_M = 0.3048

STATIONS = {
    "colville":      ("9496345", "Cape Halkett"),
    "kuparuk":       ("9497645", "Prudhoe Bay"),
    "sagavanirktok": ("9497778", "Cross Island (Dinkum Sands)"),
    "canning":       ("9498306", "Point Thomson"),
}

# short-period tidal constituents; SA/SSA/MM/MSF (seasonal, long-period) excluded
TIDAL = {"M2", "S2", "N2", "K2", "2N2", "NU2", "L2", "T2", "MU2", "LAM2",
         "K1", "O1", "P1", "Q1", "J1", "M1", "OO1", "RHO1", "2Q1",
         "M4", "MN4", "MS4", "M6", "M8", "S4", "S6"}
AMP_FLOOR_M = 0.005   # drop constituents below 5 mm -- noise at these stations


def fetch(sid):
    u = (f"https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/"
         f"stations/{sid}/harcon.json")
    return json.load(urllib.request.urlopen(u, timeout=60))["HarmonicConstituents"]


def main():
    out = {}
    for site, (sid, name) in STATIONS.items():
        cons = fetch(sid)
        kept = []
        for c in cons:
            if c["name"] not in TIDAL:
                continue
            amp = c["amplitude"] * FT_TO_M
            if amp < AMP_FLOOR_M:
                continue
            kept.append({"name": c["name"], "amp_m": round(amp, 5),
                         "phase_deg": round(c["phase_GMT"], 2),
                         "speed_deg_hr": round(c["speed"], 6)})
        kept.sort(key=lambda c: -c["amp_m"])
        rng = 2.0 * sum(c["amp_m"] for c in kept)
        out[site] = {"station_id": sid, "station_name": name,
                     "n_constituents": len(kept), "approx_range_m": round(rng, 3),
                     "constituents": kept}
        print(f"{site:15s} {name:26s} {len(kept):2d} constituents, "
              f"range ~{rng:.2f} m (M2 amp {kept[0]['amp_m']:.3f} m)")

    OUT.write_text(json.dumps(out, indent=1))
    print(f"\nwrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
