"""
Compare the multi-channel experiments against the lumped Colville baseline.

Runs compared (all Colville, year 2, open water only):
  baseline   runs/definitive/colville            lumped, B_lb = delta sum, B_ub = per-channel
  mc_lumped  .../experiment_multichannel/...     lumped, CGEM_MULTICHANNEL=on (B is total
                                                 conveyance; dispersion uses per-thread width)
  dist_main  .../experiment_multichannel/...     the 1310 m distributary, 96.7% of Q
  dist_minor .../experiment_multichannel/...     the  240 m distributary,  3.3% of Q

Reports the two things the experiments were designed to move: mouth salinity (the
open problem -- observed grabs are 8-32 PSU against a lumped model near 0) and the
carbon flux, both per-area and basin-integrated (where the surface-area correction
enters). Usage: python tools/compare_multichannel.py
"""
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

ROOT = Path(__file__).resolve().parent.parent
_EXP = ROOT / "runs/experiment_multichannel"
# Thread-count sweep. n=1.00 IS the baseline (single lumped channel); 2.49 is the
# SWORD-observed value (raw prismatic median 1052 m / B_ub 423 m) and 1.87/3.20 bracket
# it using that median's SWORD IQR (792-1354 m). n=3.50 is the superseded first guess,
# reconstructed as B_ub x median n_chan_mod -- kept as the top of the curve.
SWEEP = {1.00: ROOT / "runs/definitive/colville/output.nc",
         1.87: _EXP / "n187/output.nc",
         2.49: _EXP / "n249/output.nc",
         3.20: _EXP / "n320/output.nc",
         3.50: _EXP / "mc_lumped/output.nc"}
RUNS = {
    "baseline":   ROOT / "runs/definitive/colville/output.nc",
    "mc_lumped":  _EXP / "n249/output.nc",          # the observed-thread-count run
    "dist_main":  _EXP / "dist_main2/output.nc",
    "dist_minor": _EXP / "dist_minor2/output.nc",
}
# Q shares, from sites/_colville_delta.py -- used to area-weight the two distributaries
# back into one delta-average for comparison against the lumped runs.
SHARE = {"dist_main": 0.9665, "dist_minor": 0.0335}

YEAR2_S = 365 * 86400.0
MOUTH = 1          # first interior cell (index 0 is the boundary itself)
SEC_PER_YR = 365 * 86400.0
MMOL_C_TO_G = 12.011e-3


def load(path):
    if not path.exists():
        raise SystemExit(f"missing run: {path}\n(run the experiment first)")
    d = Dataset(path)
    t = d.variables["time"][:]
    y2 = t >= YEAR2_S
    ice = d.variables["ice_frac"][y2, :]
    doy = (t[y2] - YEAR2_S) / 86400.0
    return {
        "summer": (doy >= 180) & (doy <= 260),
        "S": np.asarray(d.variables["S"][y2, :]),
        "FCO2": np.asarray(d.variables["FCO2"][y2, :]),
        "U": np.asarray(d.variables["velocity"][y2, :]),
        "depth": np.asarray(d.variables["depth"][y2, :]),
        "width": np.asarray(d.variables["width"][y2, :]),
        "open": np.asarray(ice) < 0.5,
        "n": int(y2.sum()),
    }


def summarize(name, r):
    """Year-2 summary for one run.

    Two flux conventions are reported because they answer different questions and
    differ by ~3x on these rivers:
      FCO2_open   mean over OPEN-WATER cells/times only -- the intensity of outgassing
                  while the river is actually ice-free.
      FCO2_ann    mean over ALL of year 2 -- ice-covered cells contribute zero (biogeo
                  scales gas exchange by 1-ice_frac), so this is the annualized rate.
    FCO2 is stored per VOLUME [mmol C m^-3 s^-1]; multiplying by depth recovers the
    per-AREA flux, which is what RCO2 was before biogeo divided it by DEPTH.
    """
    op = r["open"]
    Sm = r["S"][:, MOUTH][op[:, MOUTH]]
    # Salt intrusion length: furthest cell upstream reaching 1 PSU, per timestep.
    # Reported as the open-water mean and the seasonal maximum.
    above = r["S"] > 1.0
    Lint = np.where(above.any(axis=1), above.shape[1] - 1 - above[:, ::-1].argmax(axis=1), 0) * 0.2
    Lo = Lint[op[:, MOUTH]]
    # Summer window (days 180-260 of year 2) -- the season the observed grabs that the
    # lumped model cannot reproduce were taken in, and when high flow flushes the mouth.
    su = r["summer"]
    Ss = r["S"][:, MOUTH][su]
    # FCO2 carries NaN in a small number of FULLY ice-covered cells (ice_frac == 1,
    # ~0.7% of year-2 cells here) where biogeo's openfac is 0. The flux there is
    # physically zero, so zero-fill rather than nan-skip -- nan-skipping would drop
    # them from the annual mean and inflate it.
    per_area = np.nan_to_num(r["FCO2"] * r["depth"], nan=0.0)   # mmol C m^-2 s^-1
    to_g = SEC_PER_YR * MMOL_C_TO_G
    f_open = float(np.mean(per_area[op])) * to_g if op.any() else np.nan
    f_ann = float(np.mean(per_area)) * to_g
    # basin-integrated: local water-surface area x the annualized per-area flux
    dx = 200.0
    cell_area = np.mean(r["width"], axis=0) * dx         # m^2 per cell
    f_tot = float(np.sum(np.mean(per_area, axis=0) * cell_area)) * to_g / 1e6  # tC/yr
    return {
        "S_mouth_med": float(np.median(Sm)) if Sm.size else np.nan,
        "S_mouth_p95": float(np.percentile(Sm, 95)) if Sm.size else np.nan,
        "S_mouth_max": float(Sm.max()) if Sm.size else np.nan,
        "S_frac_gt1": float((Sm > 1).mean() * 100) if Sm.size else np.nan,
        "U_mouth": float(np.mean(np.abs(r["U"][:, MOUTH][op[:, MOUTH]]))) if Sm.size else np.nan,
        "S_summer_med": float(np.median(Ss)) if Ss.size else np.nan,
        "S_summer_max": float(Ss.max()) if Ss.size else np.nan,
        "Lint_mean": float(Lo.mean()) if Lo.size else np.nan,
        "Lint_max": float(Lint.max()),
        "FCO2_open": f_open,
        "FCO2_ann": f_ann,
        "FCO2_total": f_tot,          # tonnes C yr^-1
        "surf_km2": float(np.sum(cell_area) / 1e6),
    }


def main():
    res = {k: summarize(k, load(p)) for k, p in RUNS.items()}

    # area-weighted delta average of the two distributaries
    wsum = {k: res[k]["surf_km2"] for k in SHARE}
    tot_a = sum(wsum.values())
    def wmean(key):
        return sum(res[k][key] * wsum[k] for k in SHARE) / tot_a
    res["delta_avg"] = {
        "S_mouth_med": wmean("S_mouth_med"), "S_mouth_p95": wmean("S_mouth_p95"),
        "S_mouth_max": max(res[k]["S_mouth_max"] for k in SHARE),
        "S_frac_gt1": wmean("S_frac_gt1"), "U_mouth": wmean("U_mouth"),
        "S_summer_med": wmean("S_summer_med"),
        "S_summer_max": max(res[k]["S_summer_max"] for k in SHARE),
        "Lint_mean": wmean("Lint_mean"),
        "Lint_max": max(res[k]["Lint_max"] for k in SHARE),
        "FCO2_open": wmean("FCO2_open"), "FCO2_ann": wmean("FCO2_ann"),
        "FCO2_total": sum(res[k]["FCO2_total"] for k in SHARE),
        "surf_km2": tot_a,
    }

    hdr = (f"{'run':12s}{'U_mouth':>9s}{'S_med':>7s}{'S_p95':>7s}{'S_max':>7s}{'S>1%':>7s}"
           f"{'gC/m2/yr':>10s}{'(annual)':>10s}{'tC/yr':>9s}{'surf km2':>10s}")
    print("\nColville — year 2.  Salinity at the mouth cell, open water only.\n" + "=" * len(hdr))
    print(hdr)
    print("-" * len(hdr))
    for k in ["baseline", "mc_lumped", "dist_main", "dist_minor", "delta_avg"]:
        r = res[k]
        print(f"{k:12s}{r['U_mouth']:9.3f}{r['S_mouth_med']:7.2f}{r['S_mouth_p95']:7.2f}"
              f"{r['S_mouth_max']:7.2f}{r['S_frac_gt1']:7.1f}{r['FCO2_open']:10.1f}"
              f"{r['FCO2_ann']:10.1f}{r['FCO2_total']:9.0f}{r['surf_km2']:10.2f}")
    hdr2 = (f"{'run':12s}{'S_sum_med':>11s}{'S_sum_max':>11s}{'Lint_mean km':>14s}{'Lint_max km':>13s}")
    print("\nSalt intrusion (summer = days 180-260; Lint = furthest cell above 1 PSU)\n"
          + "=" * len(hdr2))
    print(hdr2)
    print("-" * len(hdr2))
    for k in ["baseline", "mc_lumped", "dist_main", "dist_minor", "delta_avg"]:
        r = res[k]
        print(f"{k:12s}{r['S_summer_med']:11.2f}{r['S_summer_max']:11.2f}"
              f"{r['Lint_mean']:14.2f}{r['Lint_max']:13.2f}")
    print()
    b = res["baseline"]
    for k in ["mc_lumped", "delta_avg"]:
        r = res[k]
        print(f"{k:>10s} vs baseline:  per-area flux {r['FCO2_ann']/b['FCO2_ann']:.2f}x   "
              f"basin total {r['FCO2_total']/b['FCO2_total']:.2f}x   "
              f"surface area {r['surf_km2']/b['surf_km2']:.2f}x")
    # --- sensitivity of the salinity gain to the one uncertain number ------------
    have = {n: p for n, p in SWEEP.items() if p.exists()}
    if len(have) > 1:
        hdr3 = (f"{'n_chan_up':>10s}{'B_up (m)':>10s}{'S_sum_med':>11s}{'S_sum_max':>11s}"
                f"{'Lint_mean km':>14s}{'gC/m2/yr':>10s}{'tC/yr':>9s}")
        print("\nSensitivity to the prismatic thread count (the one uncertain number)\n"
              + "=" * len(hdr3))
        print(hdr3)
        print("-" * len(hdr3))
        for n in sorted(have):
            r = summarize(f"n{n}", load(have[n]))
            note = ""
            if abs(n - 2.49) < 0.01:
                note = "  <- SWORD observed"
            elif abs(n - 1.0) < 0.01:
                note = "  <- baseline (lumped)"
            elif abs(n - 3.5) < 0.01:
                note = "  <- superseded first guess"
            print(f"{n:10.2f}{423.0*n:10.0f}{r['S_summer_med']:11.2f}{r['S_summer_max']:11.2f}"
                  f"{r['Lint_mean']:14.2f}{r['FCO2_ann']:10.1f}{r['FCO2_total']:9.0f}{note}")
        print("\n1.87 and 3.20 bracket the observed 2.49 using the SWORD IQR of the raw\n"
              "prismatic width (792-1354 m). They are the observational spread, not a\n"
              "free parameter range.")

    print("\nObserved mouth-salinity grabs for reference: 8-32 PSU (all years; "
          "possibly seaward-lagoon sampled).")
    return res


if __name__ == "__main__":
    main()
