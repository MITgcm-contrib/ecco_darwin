# verification/wad_mudflat@wadcheckin

**vs its upstream base (merge-base):** new (not in its upstream base)

## README (first 40 lines)
```
# wad_mudflat: Macrotidal mudflat with tidal creeks and a river

An idealized macrotidal mudflat (2 m M2 tide) cut by a network of tidal
creeks, with a 20 m³/s river at the head of the main creek and passive
tracers marking river and flat water. Most of the flat dries every tide.
The variants add the real fresh-water flux.

## Set-up

- 100 × 60 columns of 100 m (4 × 2 tiles of 25 × 30), 5 levels of 2.5 m, r*.
- Offshore basin 10 m deep (x < 2 km), slope to −3 m, flat rising to +2 m
  at 10 km (~1:1400); a main creek along y = 3 km with 6 side branches.
- Tide at the western boundary (Flather, `code/obcs_calc.F`); river as
  `addMass` (`selectAddFluid=1`, fresh); river tracer relaxed by pkg/rbcs.
- Δt = 6 s, 7452 steps = 1 M2 cycle; implicit bottom drag
  (`selectImplicitDrag=2`).
- Packages: `gfd obcs kpp ptracers rbcs diagnostics wad`.

## Variants

| Variant | Change |
|---|---|
| `rain` | `useRealFreshWaterFlux` with extreme surface fresh water (`empmr.bin`: 20 mm/h evaporation, a 50 mm/h storm in hours 5–7): tests the evaporation limiter on the drying flats |
| `rfwf0` | `useRealFreshWaterFlux` with zero fresh-water flux: must reproduce the main case |

## Inputs

All binary inputs are in the repository; to regenerate them: `input/gendata.py` writes the grid, initial fields, river (`addMass.bin`),
tracers and rbcs files in `input/`.

## Build and run

With testreport (runs the main case and every `input.*` variant and compares
with `results/`), from MITgcm's `verification/` directory:

```sh
./testreport -t wad_mudflat            # add -of <optfile>, -mpi as usual
```

By hand:
... (truncated)
```

## Build variants (code*/)
- **code**: packages.conf = `gfd obcs kpp ptracers rbcs diagnostics wad`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor obcs kpp ptracers rbcs diagnostics wad
  - SIZE.h: grid 100x60x5; sNx=25, sNy=30, OLx=3, OLy=3, nSx=4, nSy=2, nPx=1, nPy=1, Nr=5
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h, PTRACERS_SIZE.h
  - modified/extra source: balzano_get_etan.F, obcs_calc.F, obcs_wad_offshore.F

## Input variants (input*/)
- **input**: data.pkg on: useOBCS, useKPP, usePTRACERS, useRBCS, useWAD
  - data: nIter0=0, nonlinFreeSurf=4, eosType='LINEAR', usingCartesianGrid=.TRUE., momStepping=.TRUE.
  - namelist files: data data.kpp data.obcs data.pkg data.ptracers data.rbcs data.wad eedata
- **input.rain**: data.pkg on: -
  - data: nIter0=0, nonlinFreeSurf=4, eosType='LINEAR', usingCartesianGrid=.TRUE., momStepping=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data
- **input.rfwf0**: data.pkg on: -
  - data: nIter0=0, nonlinFreeSurf=4, eosType='LINEAR', usingCartesianGrid=.TRUE., momStepping=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data

## Reference results
`output.rain.txt` `output.rfwf0.txt` `output.txt`

Run: `cd verification; ./testreport -of <optfile> -t wad_mudflat` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
