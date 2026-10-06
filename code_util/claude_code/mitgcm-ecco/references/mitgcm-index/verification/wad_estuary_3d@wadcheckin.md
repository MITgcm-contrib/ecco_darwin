# verification/wad_estuary_3d@wadcheckin

**vs its upstream base (merge-base):** new (not in its upstream base)

## README (first 40 lines)
```
# wad_estuary_3d: 3-D tidal estuary with drying flats

A 3-D estuary: a 10 m deep channel shoaling to 3 m at its head, with tidal
flats on both sides that dry at low water, stratified salinity, Coriolis,
KPP and GM/Redi, forced by a 1.5 m M2 tide through a Flather boundary. The
variants add passive tracers, the carry option and a cross-estuary
seiche.

## Set-up

- 40 × 30 columns of 200 m (2 × 2 tiles of 20 × 15), 5 levels of 2.6 m, r*.
- Channel along x (|y − 3 km| < 600 m); flats rising to +1.3 m at the side
  walls; no flats in the first 1 km so the open boundary is in deep water.
- Salinity 30 at the surface to 34 at 10 m; f = 1e-4 s⁻¹; model rest level
  3 m above sea level; starts at high water, at rest.
- Δt = 15 s (the Flather boundary needs c Δt/Δx < 1), 2981 steps = 1 M2
  cycle.
- Packages: `gfd obcs kpp gmredi ptracers diagnostics wad`
  (`code/obcs_calc.F`: Flather tide; `code/obcs_wad_offshore.F`: offshore
  water entering on the flood).

## Variants

| Variant | Change |
|---|---|
| `carry` | `wadCarryVel=.TRUE.` with `momImplVertAdv=.TRUE.` |
| `ptr` | two passive tracers: a uniform one (must stay uniform) and a dye patch |
| `seiche` | the high-water surface tilted by 1e-4 across the estuary (±0.29 m at the walls) and released: a seiche over the drying flats; 480 steps (2 h), monitor every 10 min |

## Inputs

All binary inputs are in the repository; to regenerate them: `input/gendata.py` writes `input/` (`bathy.bin`, `eta.bin`, `S.bin`),
`input.ptr/` (`ptr1.bin`, `ptr2.bin`) and `input.seiche/eta_seiche.bin`. `WAD_NX`, `WAD_NY`,
`WAD_NR` refine the grid for runs outside testreport.

## Build and run

With testreport (runs the main case and every `input.*` variant and compares
with `results/`), from MITgcm's `verification/` directory:

... (truncated)
```

## Build variants (code*/)
- **code**: packages.conf = `gfd obcs kpp gmredi ptracers diagnostics wad`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor obcs kpp gmredi ptracers diagnostics wad
  - SIZE.h: grid 40x30x5; sNx=20, sNy=15, OLx=3, OLy=3, nSx=2, nSy=2, nPx=1, nPy=1, Nr=5
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h, PTRACERS_SIZE.h
  - modified/extra source: balzano_get_etan.F, obcs_calc.F, obcs_wad_offshore.F

## Input variants (input*/)
- **input**: data.pkg on: useOBCS, useKPP, useGMRedi, useWAD
  - data: nIter0=0, nonlinFreeSurf=4, eosType='LINEAR', usingCartesianGrid=.TRUE., momStepping=.TRUE.
  - namelist files: data data.gmredi data.kpp data.obcs data.pkg data.wad eedata
- **input.carry**: data.pkg on: -
  - data: nIter0=0, nonlinFreeSurf=4, eosType='LINEAR', usingCartesianGrid=.TRUE., momStepping=.TRUE.
  - namelist files: data data.wad
- **input.ptr**: data.pkg on: useOBCS, useKPP, useGMRedi, usePTRACERS, useDiagnostics, useWAD
  - namelist files: data.diagnostics data.pkg data.ptracers
- **input.seiche**: data.pkg on: -
  - data: nIter0=0, nonlinFreeSurf=4, eosType='LINEAR', usingCartesianGrid=.TRUE., momStepping=.TRUE.
  - namelist files: data

## Reference results
`output.carry.txt` `output.ptr.txt` `output.seiche.txt` `output.txt`

Run: `cd verification; ./testreport -of <optfile> -t wad_estuary_3d` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
