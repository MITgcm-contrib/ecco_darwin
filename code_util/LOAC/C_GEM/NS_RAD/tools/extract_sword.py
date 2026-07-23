"""
Recover raw SWORD v17b channel widths for the four North Slope rivers and write
docs/sword_widths.json — the delta-mouth distributary SUM (the conveyance/salt-exchange
width) plus a raw width-vs-distance profile for the geometry document overlay.

WHY the sum. SWORD traces individual channels. A delta discharges through several
distributaries in PARALLEL that never rejoin, so the estuarine mouth width is the SUM of
the seaward-most distinct channels, not the single main thread (which SWORD keeps narrow,
and which for Colville/Kuparuk it stops tracing ~30-50 km inland — the delta below is
unnamed 'NODATA' channels). We therefore work by REACH (each reach_id = one channel
segment) in a seaward band, and sum one width per distinct reach.

Input:  the extracted na_sword_v17b.nc (SWORD_PATH below).
Output: docs/sword_widths.json  { river: {delta_sum_m, n_chan, chan_widths,
                                          dist_km:[...], width_raw:[...]} }

Re-run only when the SWORD file changes; the site configs then read delta_sum_m.
"""
import json
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

ROOT = Path(__file__).resolve().parent.parent
SWORD_PATH = Path("/tmp/sword/netcdf/na_sword_v17b.nc")
OUT = ROOT / "docs" / "sword_widths.json"

# coastal mouth point + lon half-window for each river's delta region
MOUTH = {"colville": (-150.7, 70.44, 0.65),
         "kuparuk": (-148.9, 70.42, 0.45),
         "sagavanirktok": (-148.05, 70.31, 0.40),
         "canning": (-145.85, 70.13, 0.40)}
BAND_DEG = 0.03      # seaward latitude band (~3 km) over which distributaries are summed
DOMAIN_KM = 27.2     # model domain length


def main():
    d = Dataset(SWORD_PATH)
    n = d.groups["nodes"]
    x = n.variables["x"][:]; y = n.variables["y"][:]
    w = n.variables["width"][:]; lake = n.variables["lakeflag"][:]
    rid = n.variables["reach_id"][:]

    nch = np.asarray(n.variables["n_chan_mod"][:], float)
    so = np.asarray(n.variables["stream_order"][:])
    trib = np.asarray(n.variables["trib_flag"][:])
    w = np.asarray(w, float); x = np.asarray(x, float); y = np.asarray(y, float)

    out = {}
    for name, (mx, my, hw) in MOUTH.items():
        box = ((x > mx - hw) & (x < mx + hw) & (y > 70.0) & (y < my + 0.15)
               & (lake == 0) & (w > 0))
        xi, yi, wi, ri = x[box], y[box], w[box], rid[box]

        # --- delta-mouth conveyance: SUM of one RAW width per distinct reach in the
        # seaward band. Raw (not per-channel) because at the delta the distributaries
        # convey the flow in parallel -- that sum is B_lb. ------------------------------
        ymax = yi.max()
        band = yi > (ymax - BAND_DEG)
        chan = [float(wi[(ri == r) & band].mean()) for r in np.unique(ri[band])]
        delta_sum = float(np.sum(chan))

        # --- width PROFILE for the doc overlay: the PER-CHANNEL width (raw / n_chan_mod),
        # main stem only (top stream orders, non-tributary). Raw SWORD width is the total
        # across up to ~33 braids, which for a single conveyance channel must be divided
        # by n_chan_mod; and tributaries/side channels are dropped. Binned by distance
        # from the mouth with a robust median + IQR. This is what the flare's prismatic
        # B_ub should match; the delta_sum above is the separate mouth widening. ---------
        nci, soi, ti = nch[box], so[box], trib[box]
        keep = (soi >= soi.max() - 1) & (ti == 0)
        pc = wi[keep] / np.clip(nci[keep], 1, None)
        dist = np.sqrt(((xi[keep] - mx) * np.cos(np.radians(my)) * 111.0) ** 2
                       + ((yi[keep] - my) * 111.0) ** 2)
        dk, wmed, wlo, whi = [], [], [], []
        for b in np.arange(0, DOMAIN_KM, 1.5):
            m = (dist >= b) & (dist < b + 1.5)
            if m.sum() < 2:
                continue
            dk.append(round(float(b + 0.75), 2))
            wmed.append(round(float(np.median(pc[m])), 1))
            wlo.append(round(float(np.percentile(pc[m], 25)), 1))
            whi.append(round(float(np.percentile(pc[m], 75)), 1))
        # every individual per-channel node in the domain, for the full scatter overlay
        inbox = dist <= DOMAIN_KM
        scat = sorted(zip([round(float(v), 2) for v in dist[inbox]],
                          [round(float(v), 1) for v in pc[inbox]]))

        out[name] = {"delta_sum_m": round(delta_sum),
                     "n_chan": len(chan),
                     "chan_widths": [round(c) for c in sorted(chan, reverse=True)],
                     "dist_km": dk, "width_perchan_med": wmed,
                     "width_perchan_p25": wlo, "width_perchan_p75": whi,
                     "scatter_dist_km": [p[0] for p in scat],
                     "scatter_width": [p[1] for p in scat]}
        print(f"{name:14s} delta SUM={delta_sum:.0f} m  ({len(chan)} ch "
              f"{[round(c) for c in sorted(chan, reverse=True)]})  per-chan profile "
              f"median~{int(np.median(wmed)) if wmed else 0} m, {len(dk)} bins")

    json.dump(out, open(OUT, "w"), indent=1)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
