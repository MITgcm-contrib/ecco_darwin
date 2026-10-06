# verification/wad_balzano@wadcheckin

**vs its upstream base (merge-base):** new (not in its upstream base)

## README (first 40 lines)
```
# wad_balzano: Balzano (1998) tidal-flat tests: slope, step, pool

A sloping beach driven by an M2-like tide at an open boundary, after
Balzano's (1998) benchmark tests: the shoreline must follow the tide over a
plain slope, flow over a 1 m bed step without negative depths, and leave a
pool trapped behind a sill at the sill crest (plus at most `wadCritDepth`)
at low water.

## Set-up

- 140 columns of 100 m (5 tiles of 28), 1 level of 8 m, r*.
- Bed from −5 m at the boundary to 0 m at 13.8 km; land at the end.
- A 2 m, 12 h tide imposed at the western boundary (`code/obcs_calc.F`,
  `code/balzano_get_etan.F`); model rest level 2.5 m above mean sea level.
- Δt = 10 s, 12960 steps = 3 tidal cycles of 12 h; starts at high water.
- Packages: `gfd obcs wad diagnostics`.

## Variants

| Variant | Change |
|---|---|
| `step` | a 1 m bed step half way up the beach (`input.step/bathy.bin`) |
| `pool` | a 0.5 m sill and a 1 m deep pool on the upper beach (`input.pool/bathy.bin`) |
| `zstar` | surface-level non-linear free surface (`select_rStar=0`, `hFacInf=0.005`) |

## Inputs

All binary inputs are in the repository; to regenerate them: `input/gendata.py` writes `input/bathy.bin`, `input.step/bathy.bin`,
`input.pool/bathy.bin` and `input/eta_init.bin`.

## Build and run

With testreport (runs the main case and every `input.*` variant and compares
with `results/`), from MITgcm's `verification/` directory:

```sh
./testreport -t wad_balzano            # add -of <optfile>, -mpi as usual
```

By hand:
... (truncated)
```

## Build variants (code*/)
- **code**: packages.conf = `gfd obcs wad diagnostics`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor obcs wad diagnostics
  - SIZE.h: grid 140x1x1; sNx=28, sNy=1, OLx=3, OLy=3, nSx=5, nSy=1, nPx=1, nPy=1, Nr=1
  - option/size headers: CPP_OPTIONS.h
  - modified/extra source: balzano_get_etan.F, obcs_calc.F

## Input variants (input*/)
- **input**: data.pkg on: useOBCS, useWAD
  - data: deltaT=10., nTimeSteps=12960, nIter0=0, nonlinFreeSurf=4, select_rStar=2, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., saltAdvScheme=33, momStepping=.TRUE.
  - namelist files: data data.obcs data.pkg data.wad eedata
- **input.pool**: data.pkg on: -
- **input.step**: data.pkg on: -
- **input.zstar**: data.pkg on: -
  - data: deltaT=10., nTimeSteps=12960, nIter0=0, nonlinFreeSurf=4, select_rStar=0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., saltAdvScheme=33, momStepping=.TRUE.
  - namelist files: data

## Reference results
`output.pool.txt` `output.step.txt` `output.txt` `output.zstar.txt`

Run: `cd verification; ./testreport -of <optfile> -t wad_balzano` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
