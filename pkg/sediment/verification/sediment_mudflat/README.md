# sediment_mudflat: river sediment on a drying macrotidal mudflat

`pkg/sediment` with `pkg/wetdry`. This is the `wetdry_mudflat` experiment
of `../../wetdry` (same grid, tide, river and code) with two sediment
classes added. Over one M2 tide, a river delivers mud and fine sand to the
head of the main creek. They settle, deposit on the flats, and erode
again in the creek on the ebb.

## Set-up

As `wetdry_mudflat` (see `../../wetdry/verification/wetdry_mudflat/README.md`):

- 100 × 60 columns of 100 m (4 × 2 tiles of 25 × 30), 5 levels of 2.5 m, r\*.
- Offshore basin 10 m deep, a flat that rises to +2 m, a creek network.
- 2 m M2 tide at the western boundary (Flather); a 20 m^3/s river
  (`addMass`) at the head of the main creek.
- Δt = 6 s, 7452 steps = 1 M2 cycle.
- Packages: `gfd obcs kpp ptracers rbcs diagnostics wetdry sediment`.

Sediment (`input/data.sediment`; illustrative values in typical ranges,
not calibrated):

- class 1, mud (pTracer 3): w_s = 0.5 mm/s, particle τ_ce = 0.05 Pa,
  0.3 kg/m^3 in the river (6 kg/s), cohesive
- class 2, fine sand (pTracer 4): w_s = 10 mm/s, τ_ce = 0.15 Pa,
  0.05 kg/m^3 in the river (1 kg/s)
- clean bed at t = 0; layered bed (8 layers, porosity 0.9, consolidation
  1 d, swelling 10 d), continuous deposition
- pTracers 1 and 2 mark river water and flat water, as in `wetdry_mudflat`

## Variants

| Variant | Change |
|---|---|
| `waves` | adds a parametric wind sea (`sedWaveMode=2`, Young and Verhagen 1996): a steady 6 m/s onshore wind over a 10 km fetch, depth-limited breaking at Hs = 0.5 h |

## Inputs

All binary inputs are in the repository. To make them again, run
`python3 gendata.py` in `input/` (the same script as `wetdry_mudflat`).

## Build and run

Install pkg/wetdry and pkg/sediment first (see `../../README.md`). Then,
from MITgcm's `verification/` directory:

```sh
./testreport -t sediment_mudflat
```

By hand: as for `sediment_seiche` (see its README).

## What to check

- `%WETDRY_MON`: volume budget error at round-off, as in `wetdry_mudflat`.
- `%SEDIMENT_MON`: after one tide, river input = suspended + bed for each
  class, and the step budget error stays at round-off (~3e-12 kg).
  Main run: of 2.68e5 kg of river mud, 2.22e5 kg is on the bed and
  0.47e5 kg in suspension; of 4.47e4 kg of sand, all but ~40 kg is on
  the bed.
- `waves`: the waves keep the flats stirred. Almost all river sediment
  (2.68e5 kg of mud, 4.47e4 kg of sand) is in suspension at the end of
  the tide.

Compare the `%MON` lines with `results/output*.txt` (testreport does
this).

Reference outputs were made with testreport (default `-ieee`) on arm64
macOS with conda-forge gfortran (optfile `darwin_arm64_gfortran_conda`),
using pkg/wetdry from `../../wetdry/src`. The development branch has a
newer pkg/wetdry (with sea-ice options); there the same runs agree with
these to 4–5 digits.
