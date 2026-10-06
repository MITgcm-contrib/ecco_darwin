# verification/shelfice_2d_remesh

## README (first 40 lines)
```
Simplified experiment to test pkg/shelfice vertical remeshing code
(simplified version of experiment verification_other/shelfice_remeshing/)
===============================================================================

Specific options:
* Use a 2-D (y z) slice southern ocean (lat-long grid) domain with
  flat bottom and Open Boundary at Northern edge; uniform horiz and
  vertical resolution.
* shelfice pkg is used with SHELFICEuseGammaFrict=T but without
  SHELFICEboundaryLayer. To avoid horizontal noise in solution next to a
  large jump in top grid-cell thickness, vertical viscosity and diffusivity
  are increased wherever the top grid-cell is too thin (pCellMix_select=20)
  and are treated implicitly (implicitDiffusion=implicitViscosity=T).
* Use Non-Linear Free surface formulation with real fresh-water flux.
* Ice-Shelf has a simple initial shape (trapezoidal) and evolves
  (SHELFICEMassStepping=T) as a result of melting and prescribed external
  tendency (SHELFICEMassDynTendFile).

IMPORTANT:
  In order to experience several remeshing event during a very short test run,
  the ice-mass tendency forcing has been set to to un-realistically large value.

Input files:
* generated from matlab script: gendata.n

Sequence of runs:
* From resting initial conditions, model was integrated for 1 day (-> iter=2880)
  without SHELFICEMassDynTendFile.
  Then, with SHELFICEMassDynTendFile, run just 18.it to generated current pickup
  files (at iter=2898).
* During current short test run (20.iter long), 4 remeshing event occurs,
  ( grep -A4 'SHI_REMESH at' output.txt )
  2 top-cell merge (it= 2900 & 2916) and 2 top-cell split (it= 2904 & 2914).
```

## Build variants (code*/)
- **code**: packages.conf = `gfd obcs shelfice diagnostics`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor obcs shelfice diagnostics
  - SIZE.h: grid 1x200x90; sNx=1, sNy=50, OLx=3, OLy=3, nSx=1, nSy=4, nPx=1, nPy=1, Nr=90; has SIZE.h_mpi
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h, DIAG_OPTIONS.h, OBCS_OPTIONS.h, SHELFICE_OPTIONS.h

## Input variants (input*/)
- **input**: data.pkg on: useOBCS, useShelfIce, useDiagnostics
  - data: deltaT=300.0, nTimeSteps=20, nIter0=2898, nonlinFreeSurf=4, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE., usingCartesianGrid=.FALSE., tempAdvScheme=77, saltAdvScheme=77, useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.diagnostics data.obcs data.pkg data.shelfice eedata

## Reference results
`output.txt`

Run: `cd verification; ./testreport -of <optfile> -t shelfice_2d_remesh` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
