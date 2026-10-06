# Where things live: configurations and projects

All under `~/Documents/research/` unless noted. Project `CLAUDE.md` files and the memory dirs
`~/.claude/projects/-Users-<you>-Documents-research-ECCO-*/memory/` hold the latest detail.

## Local MITgcm / darwin3 clones (keep separate; each has its own branch work)

| Path | Branch / purpose |
|---|---|
| `ECCO/BBL/MITgcm` | pkg/bbl development (detrainment etc.) on MITgcm master |
| `ECCO/BBL/MITgcm_c68g` | pkg/bbl backported to checkpoint68g (`bbl-c68g`) for ECCO v4r5–r7 |
| `ECCO/wetting_drying/MITgcm` | pkg/wad, pkg/overflood, pkg/mangrove, pkg/sediment development, incl. sea ice |
| `ECCO/MITgcm_wad_checkin` | pkg/wad clean check-in branch (`wad`, no sea ice) |
| `ECCO/sea_ice_BCs/MITgcm` | OBCS sea-ice sponge (`ALLOW_OBCS_SEAICE_SPONGE`) |
| `~/Documents/GitHub/darwin3` | darwin3, branch `darwin` (default; reference source) |
| `research/debug/darwin3` | darwin3 `backport_ckpt68y` |
| `research/emulation/darwin3` | another darwin3 tree |
| `research/darwin3` | not a clone: a Darwin analysis project (m_files, python, model_setup, figures) |

## Configurations

| Config | Key facts | Where |
|---|---|---|
| ECCO v4r6 LLC90 forward | checkpoint68g; 113 ranks × 30×30, Nr=50, dt=3600 s, z* (`select_rStar=2`, `nonlinFreeSurf=4`), JMD95Z; seaice, shelfice, salt_plume, ggl90, gmredi, exch2, autodiff/ctrl/ecco; 1992–2025 (`nTimeSteps=298031`), ~36–39 h | `ECCO/offline/code_v4r6_forward`, `run_v4r6_forward` |
| Offline Darwin on LLC90 (V4r6) | open mismatches vs Oliver Jahn's set-up are tracked in the ECCO-offline memory dir; darwin3 `backport_ckpt68y`; pkg/offline reads the daily v4r6 archive (converted from MDS); custom `ggl90_check.F`; 36 tracers (6 phyto + 4 zoo + Chl), 13-band radtrans, OASIM; published as `ecco_darwin/offline/V4r6_darwin_offline` | `ECCO/offline/{code_offline_ggl90,code_6+4+0_llc90_ggl90,run_6+4+0_llc90_offline,repo}`; `COMPILE_AND_RUN.md`, `EXTERNAL_DATA_NEEDED.md` |
| Steph's 1° offline Darwin | 360×160×23 lat-lon, 96 ranks, 360-day calendar, KPP, ECCO it3.73 physics | `ECCO/offline/from_steph`; `ECCO/CLAUDE.md` |
| ECCO-Darwin v05 / v06 | v05: LLC90 V4r5, 31 tracers, nTimeSteps 289282; v06: 36 tracers + SOLVESAPHE + rivers, 1° V4r5 (`debug/ecco_darwin/v06/1deg`, nTimeSteps 280498) and LLC270; 3° testbed 128×64×15 | Pleiades `/nobackup/<nas_user>/v05_1deg_V4r5`; ecco_darwin repo; `research/debug/` |
| pkg/bbl tests | `global_oce_latlon` (4°), `global_ocean.cs32x15` (+seaice), cone/polynya idealized (polynya in `BBL/movies/polynya`), isomip shelfice, LLC90 v05 10-yr ctl vs bbl, adjoint tests | `ECCO/BBL/{tests,plans,docs,figures}` |
| Wetting/drying | Colville delta (Alaska Albers 100 m, 640×560×10, OBCS Flather, KPP, Manning drag); LR17 Pertuis Charentais (80 m, IBI boundaries); idealized wad_* experiments | `ECCO/wetting_drying/{colville_model,LR17_model,runs,scripts}` |
| Sea-ice open boundaries | idealized twin, synmac on LLC270 grid, downscaling and lat-lon analysis | `ECCO/sea_ice_BCs/{idealized,synmac,downscale,latlon_analysis}` |
| Mackenzie delta regional | Mac270/llc4320 regional; 1-km lat-lon (882 ranks) | `research/emulation/ecco_darwin/regions/mac_delta`; pfe `/nobackup/<nas_user>/sea_ice_BCs_latlon` |

Packages seen across configs: exf, cal, obcs, ptracers, gchem, darwin, radtrans, offline,
diagnostics, mnc, layers, salt_plume, ggl90, kpp, gmredi, seaice, shelfice, rbcs, profiles,
down_slope (mutually exclusive with bbl), bbl, wad.

## Data sources

ECCO V4r6 ancillary data (PO.DAAC; copy on Pleiades under `/nobackup/owang/...`), GLODAPv2,
WOA/Levitus, OASIM (Oliver Jahn), ERA5 (cdsapi), NOAA GFS-Wave (AWS), IBI, SHOM tide gauges.
For the sea-ice boundaries, use 6-hourly means of ice velocity, never snapshots;
`SEAICEnonLinIterMax` matters as much as the boundary values.
