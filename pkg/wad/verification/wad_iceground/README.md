# wad_iceground: Grounded sea ice that melts and floats

A closed tidal flat covered by 0.8 m of sea ice that loads the surface
(real fresh-water flux). Where the water is shallower than the ice draft
(0.71 m) the column sits at the WAD film: the ice is grounded. A warm,
constant atmosphere melts the ice; as it thins, water flows back under it
and the grounding line moves up the flat.

## Set-up

- Closed basin, 60 × 10 columns of 50 m, 5 r* levels of 0.6 m; model rest
  level 1 m above MSL.
- Bed from −2 m MSL (west) to +0.3 m (east); ice everywhere (AREA 1).
- Atmosphere (`data.exf`): 15 °C, 300 W/m² short and long wave, 5 m/s wind.
- Δt = 60 s, 14400 steps = 10 days. Linear EOS, f = 0.
- Packages: `gfd diagnostics exf cal seaice wad`; `wadIceTopMelt=.TRUE.`.

## Variants

| Variant | Change |
|---|---|
| `dyn` | ice dynamics on (LSR, basal drag), ice only east of 1 km over a deeper shelf, open water to drift into, Coriolis on, ice load passed through `wadLoadDepth = 0.3 m`; Δt = 20 s, 43200 steps |

## Inputs

All binary inputs are in the repository; `input/gendata.py` regenerates
them (`bathy.bin`, `eta.bin`, `heff.bin`, `area.bin`, and the `_part` /
`bathy_deep` files used by `dyn`).

## What to check

- `%WAD_MON budgetErr/vol0` ~1e-15 (volume).
- Total salt (Σ S·h·A) constant to ~1e-5 over the run: the basin is closed
  and the melt water is fresh. Before the 2026-10-06 fix of the real
  fresh-water term on drying columns this case gained 51 % salt in 10 days.
