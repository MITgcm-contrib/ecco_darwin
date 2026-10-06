# verification/wad_thacker_1d@wadcheckin

**vs its upstream base (merge-base):** new (not in its upstream base)

## README (first 40 lines)
```
# wad_thacker_1d: Thacker (1981) oscillating basin

A frictionless parabolic basin in 1-D (x-z section) whose planar surface
sloshes with an exact analytic solution (Thacker 1981): the shoreline moves
over the dry beach on both sides every period. It is the accuracy test of
the package: volume and salt must be exact, and the surface error against
the exact solution measures the wetting-drying scheme alone.

## Set-up

- 250 columns of 40 m (5 tiles of 50), 1 level of 20 m, r* (`select_rStar=2`).
- Basin half-width 3 km, centre depth 10 m below Thacker's mean level; the
  model rest level is 10 m above it (`DATUM`), so the whole beach is ocean.
- Δt = T/272 ≈ 4.95 s (period T = 1346 s), 1360 steps = 5 periods; surface
  written every T/8.
- Packages: `gfd wad diagnostics`.

## Variants

| Variant | Change |
|---|---|
| `carry` | `wadCarryVel=.TRUE.`: a face that opens takes the upstream depth-mean velocity |
| `zstar` | surface-level non-linear free surface (`select_rStar=0`, `hFacInf=0.001`) |

## Inputs

All binary inputs are in the repository; to regenerate them: `input/analytic.py gen` writes `bathy_thacker.bin` and `eta_thacker.bin`
(run it in `input/`).

## Build and run

With testreport (runs the main case and every `input.*` variant and compares
with `results/`), from MITgcm's `verification/` directory:

```sh
./testreport -t wad_thacker_1d            # add -of <optfile>, -mpi as usual
```

By hand:

... (truncated)
```

## Build variants (code*/)
- **code**: packages.conf = `gfd wad diagnostics`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor wad diagnostics
  - SIZE.h: grid 250x1x1; sNx=50, sNy=1, OLx=3, OLy=3, nSx=5, nSy=1, nPx=1, nPy=1, Nr=1
  - option/size headers: CPP_OPTIONS.h

## Input variants (input*/)
- **input**: data.pkg on: useWAD
  - data: deltaT=4.947464852, nTimeSteps=1360, nIter0=0, nonlinFreeSurf=4, select_rStar=2, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., saltAdvScheme=33, momStepping=.TRUE.
  - namelist files: data data.pkg data.wad eedata
- **input.carry**: data.pkg on: -
  - namelist files: data.wad
- **input.zstar**: data.pkg on: -
  - data: deltaT=4.947464852, nTimeSteps=1360, nIter0=0, nonlinFreeSurf=4, select_rStar=0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., saltAdvScheme=33, momStepping=.TRUE.
  - namelist files: data

## Reference results
`output.carry.txt` `output.txt` `output.zstar.txt`

Run: `cd verification; ./testreport -of <optfile> -t wad_thacker_1d` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
