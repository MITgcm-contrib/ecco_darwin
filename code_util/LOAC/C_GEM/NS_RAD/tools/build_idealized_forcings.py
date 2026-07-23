#!/usr/bin/env python3
"""
Build the analytic, time-varying forcings and boundary conditions for the IDEALIZED
North Slope river (sites/idealized.py).

This is a VERIFICATION FIXTURE, not a real river. Every series here is a closed-form
function of day-of-year, chosen so the model is exercised end to end over a full
seasonal cycle with a KNOWN, reproducible input:

  * a spring FRESHET pulse in discharge large enough to trigger the hydraulic ice
    break-up (Q_peak ~ 8x the annual mean, well above config.BREAKUP_Q_FACTOR);
  * an air-temperature sinusoid that crosses 0 C twice, so the prognostic ice model
    freezes up in autumn and is available to break up on the freshet;
  * TIME-VARYING BOUNDARY CHEMISTRY tied to the hydrograph -- riverine DOC flushed
    high on the freshet, alkalinity/DIC diluted by snowmelt, nitrate drawn down in
    summer -- plus a marine salinity that freshens under summer sea-ice melt. These
    are what exercise the new file-driven BOUNDARY_FORCING mechanism (both the
    upstream `cub` and the downstream `clb` ends).

Everything is derived from the PARAMS block below, so the fixture is fully
reproducible: re-running this script reproduces the CSVs byte-for-byte.

    python tools/build_idealized_forcings.py            # write CSVs
    python tools/build_idealized_forcings.py --plot     # + a preview PNG

Output: forcing/idealized_*.csv  (single column, 365 daily values,
the format exfread expects). Filenames are referenced from sites/idealized.py.
"""
import argparse
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FORCINGS = os.path.normpath(os.path.join(HERE, "..", "forcing"))
N = 365
DAY = np.arange(N, dtype=float)          # day-of-year 0..364
TWO_PI = 2.0 * np.pi

# --- all knobs in one place; the fixture is fully defined by these -------------
PARAMS = dict(
    # discharge [m^3 s^-1]
    Q_BASE=8.0, Q_FRESHET=380.0, Q_PEAK_DAY=158.0, Q_SIG=12.0,
    Q_SUMMER=45.0, Q_OPEN=140.0, Q_CLOSE=300.0,
    # air temperature [degC]:  mean + amp*(-cos), coldest ~ AIR_PHASE
    AIR_MEAN=-11.0, AIR_AMP=21.0, AIR_PHASE=15.0,
    # sea (marine boundary) temperature [degC], floored at seawater freezing
    SEA_MEAN=1.0, SEA_AMP=4.0, SEA_PHASE=30.0, SEA_FREEZE=-1.8,
    # river temperature from air via the model's own regression (build_river_temp.py)
    RIVERT_SLOPE=1.432, RIVERT_T0=-4.0, RIVERT_WIN=10,
    # shortwave [W m^-2], polar: near zero in winter
    SOLAR_BASE=150.0, SOLAR_AMP=175.0, SOLAR_SOLSTICE=172.0,
    # wind [m s^-1]
    WIND_MEAN=6.0, WIND_AMP=1.5,
    # relative humidity [frac]
    RH_MEAN=0.84, RH_AMP=0.04,
    # atmospheric pCO2 [uatm]
    PCO2_MEAN=415.0, PCO2_AMP=8.0,
    # --- time-varying boundaries --------------------------------------------
    # riverine (cub) chemistry [mmol m^-3], tied to a broadened freshet shape
    CHEM_SIG=22.0,                         # freshet-response width (> Q_SIG)
    TOC_BASE=250.0, TOC_PULSE=500.0,       # labile DOC flushed HIGH on the freshet
    RDOC_BASE=250.0, RDOC_PULSE=550.0,     # refractory/chromophoric DOC, freshet-flushed
    ALK_BASE=900.0, ALK_DILUTE=550.0,      # alkalinity DILUTED by snowmelt
    DIC_BASE=880.0, DIC_DILUTE=540.0,      # DIC diluted with it (stays < ALK)
    NO3_BASE=4.0, NO3_AMP=3.0,             # nitrate winter-high, summer drawdown
    # marine (clb) salinity [PSU]: summer sea-ice-melt freshening
    S_BASE=30.0, S_MELT=6.0,
)


def _seasonal(mean, amp, phase, sign=-1.0):
    """mean + sign*amp*cos(2 pi (d - phase)/365).  sign=-1 => minimum at `phase`."""
    return mean + sign * amp * np.cos(TWO_PI * (DAY - phase) / N)


def _open_water(p):
    """0..1 half-sine hump over the open-water window [Q_OPEN, Q_CLOSE]."""
    x = (DAY - p["Q_OPEN"]) / (p["Q_CLOSE"] - p["Q_OPEN"])
    return np.clip(np.sin(np.pi * np.clip(x, 0.0, 1.0)), 0.0, 1.0)


def _freshet(p, sig):
    """Unit-height Gaussian freshet shape (0..1), width `sig`."""
    return np.exp(-0.5 * ((DAY - p["Q_PEAK_DAY"]) / sig) ** 2)


def _rolling_mean_wrap(x, win):
    """Trailing `win`-day mean with year wrap-around (climatology is periodic)."""
    pad = np.concatenate([x[-win:], x])
    c = np.cumsum(pad)
    out = (c[win:] - c[:-win]) / win
    return out


def build(p=PARAMS):
    """Return {filename: 365-array} for every idealized series."""
    # discharge: baseflow + Gaussian freshet + summer recession hump
    Q = (p["Q_BASE"]
         + p["Q_FRESHET"] * _freshet(p, p["Q_SIG"])
         + p["Q_SUMMER"] * _open_water(p))

    # air temperature, and river temperature from it (model's own relation)
    air = _seasonal(p["AIR_MEAN"], p["AIR_AMP"], p["AIR_PHASE"])
    air_sm = _rolling_mean_wrap(air, p["RIVERT_WIN"])
    # T_river = max(0, 1.432*(T_air_10day + 4.0)) -- the model's own regression
    # (tools/build_river_temp.py), with RIVERT_T0 = -4 C so (air_sm - T0) = air_sm + 4.
    river_t = np.maximum(0.0, p["RIVERT_SLOPE"] * (air_sm - p["RIVERT_T0"]))
    sea_t = np.maximum(p["SEA_FREEZE"], _seasonal(p["SEA_MEAN"], p["SEA_AMP"], p["SEA_PHASE"]))

    # meteorology
    solar = np.maximum(0.0, _seasonal(p["SOLAR_BASE"], p["SOLAR_AMP"],
                                      p["SOLAR_SOLSTICE"], sign=+1.0))
    wind = np.maximum(0.5, p["WIND_MEAN"] + p["WIND_AMP"] * np.sin(TWO_PI * DAY / N))
    rh = _seasonal(p["RH_MEAN"], p["RH_AMP"], p["AIR_PHASE"], sign=+1.0)  # winter-high
    pco2 = _seasonal(p["PCO2_MEAN"], p["PCO2_AMP"], p["AIR_PHASE"], sign=+1.0)

    # time-varying boundaries
    fchem = _freshet(p, p["CHEM_SIG"])          # 0..1 broadened freshet response
    ow = _open_water(p)
    toc = p["TOC_BASE"] + p["TOC_PULSE"] * fchem
    rdoc = p["RDOC_BASE"] + p["RDOC_PULSE"] * fchem   # refractory DOC, freshet-flushed
    alk = p["ALK_BASE"] - p["ALK_DILUTE"] * fchem
    dic = p["DIC_BASE"] - p["DIC_DILUTE"] * fchem
    no3 = np.maximum(0.1, _seasonal(p["NO3_BASE"], p["NO3_AMP"], p["AIR_PHASE"], sign=+1.0))
    s_marine = p["S_BASE"] - p["S_MELT"] * ow

    return {
        "idealized_discharge_m3sec.csv": Q,
        "idealized_airtemp_degC.csv": air,
        "idealized_river_watertemp_degC.csv": river_t,
        "idealized_seatemp_degC.csv": sea_t,
        "idealized_solar_Wm2.csv": solar,
        "idealized_wind_msec.csv": wind,
        "idealized_relhum_frac.csv": rh,
        "idealized_pCO2_uatm.csv": pco2,
        "idealized_S_clb_marine.csv": s_marine,
        "idealized_TOC_cub_river.csv": toc,
        "idealized_RDOC_cub_river.csv": rdoc,
        "idealized_NO3_cub_river.csv": no3,
        "idealized_ALK_cub_river.csv": alk,
        "idealized_DIC_cub_river.csv": dic,
    }


def write(series, outdir=FORCINGS):
    for fname, arr in series.items():
        assert arr.shape == (N,), f"{fname}: {arr.shape}, want ({N},)"
        assert np.all(np.isfinite(arr)), f"{fname}: non-finite value"
        # single column, 365 values, no header, no trailing newline -- matches the
        # existing forcing CSVs and exfread's genfromtxt reader.
        with open(os.path.join(outdir, fname), "w") as f:
            f.write("\n".join(f"{x:.6f}" for x in arr))
        print(f"  {fname:38s} min={arr.min():9.3f}  max={arr.max():9.3f}")


def summarize(series, p=PARAMS):
    Q = series["idealized_discharge_m3sec.csv"]
    air = series["idealized_airtemp_degC.csv"]
    qmean = Q.mean()
    print("\nfixture summary")
    print(f"  discharge: mean {qmean:.1f}, peak {Q.max():.0f} m3/s "
          f"(peak/mean = {Q.max()/qmean:.1f}x; break-up trigger is "
          f"{p_breakup()}x mean)")
    below0 = np.where(air < 0)[0]
    above0 = np.where(air >= 0)[0]
    if above0.size:
        print(f"  air > 0 C on days {above0.min()}-{above0.max()} "
              f"({above0.size} d open); frozen the rest")
    print(f"  freshet peak day {int(np.argmax(Q))}; air-0C crossings frame the "
          "ice season")


def p_breakup():
    # BREAKUP_Q_FACTOR lives in config, but importing config needs a site env var;
    # keep this script standalone and just report the shipped default.
    return 3.0


def plot(series, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    order = list(series)
    fig, axes = plt.subplots(5, 3, figsize=(13, 12))
    for ax, name in zip(axes.ravel(), order):
        ax.plot(DAY, series[name], lw=1.4)
        ax.set_title(name.replace("idealized_", "").replace(".csv", ""), fontsize=9)
        ax.set_xlim(0, 364)
        ax.grid(alpha=0.3)
    for ax in axes.ravel()[len(order):]:
        ax.axis("off")
    fig.suptitle("Idealized North Slope river -- analytic forcings & boundaries",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(path, dpi=120)
    print(f"\npreview -> {path}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--plot", action="store_true", help="also write a preview PNG")
    ap.add_argument("--outdir", default=FORCINGS)
    args = ap.parse_args()
    s = build()
    print(f"writing {len(s)} idealized forcing files to {args.outdir}")
    write(s, args.outdir)
    summarize(s)
    if args.plot:
        plot(s, os.path.join(args.outdir, "..", "docs", "idealized_forcings_preview.png"))
