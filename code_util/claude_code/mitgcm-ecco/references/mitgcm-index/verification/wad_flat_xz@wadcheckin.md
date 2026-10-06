# verification/wad_flat_xz@wadcheckin

**vs its upstream base (merge-base):** new (not in its upstream base)

## README (first 40 lines)
```
# wad_flat_xz: Stratified beach at rest, and under a tide

An x-z section of a stratified ocean over a shelf and a beach whose top is
dry at rest. With no forcing every velocity is spurious: the test measures
the pressure-gradient error of thin r* columns next to dry cells (expect
max |u| below 1 mm/s, decaying). The `tide` variant drives the same section
with an M2 tide and KPP.

## Set-up

- 100 columns of 100 m (5 tiles of 20), 10 levels of 2.5 m, r*.
- 20 m shelf for x < 5 km, beach rising to +2 m at 10 km; model rest level
  3 m above sea level.
- Temperature stratified with depth only (20 °C at the surface to 10 °C at
  20 m).
- Main case: Δt = 30 s, 2880 steps = 1 day at rest.
- Packages: `gfd obcs kpp wad`.

## Variants

| Variant | Change |
|---|---|
| `tide` | 1.5 m M2 tide at a western Flather boundary, KPP on; Δt = 6 s, one tidal cycle |

## Inputs

All binary inputs are in the repository; to regenerate them: `input/gendata.py` writes `bathy.bin`, `eta.bin` and `T.bin`.

## Build and run

With testreport (runs the main case and every `input.*` variant and compares
with `results/`), from MITgcm's `verification/` directory:

```sh
./testreport -t wad_flat_xz            # add -of <optfile>, -mpi as usual
```

By hand:

```sh
... (truncated)
```

## Build variants (code*/)
- **code**: packages.conf = `gfd obcs kpp wad`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor obcs kpp wad
  - SIZE.h: grid 100x1x10; sNx=20, sNy=1, OLx=3, OLy=3, nSx=5, nSy=1, nPx=1, nPy=1, Nr=10
  - option/size headers: CPP_OPTIONS.h
  - modified/extra source: balzano_get_etan.F, obcs_calc.F

## Input variants (input*/)
- **input**: data.pkg on: useWAD
  - data: nIter0=0, nonlinFreeSurf=4, eosType='LINEAR', usingCartesianGrid=.TRUE., momStepping=.TRUE.
  - namelist files: data data.pkg data.wad eedata
- **input.tide**: data.pkg on: useOBCS, useKPP, useWAD
  - data: nIter0=0, nonlinFreeSurf=4, eosType='LINEAR', usingCartesianGrid=.TRUE., momStepping=.TRUE.
  - namelist files: data data.kpp data.obcs data.pkg

## Reference results
`output.tide.txt` `output.txt`

Run: `cd verification; ./testreport -of <optfile> -t wad_flat_xz` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
