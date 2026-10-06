# verification/tutorial_advection_in_gyre

## Build variants (code*/)
- **code**: packages.conf = `gfd ptracers diagnostics mnc`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor ptracers diagnostics mnc
  - SIZE.h: grid 60x60x1; sNx=30, sNy=30, OLx=4, OLy=4, nSx=2, nSy=2, nPx=1, nPy=1, Nr=1; has SIZE.h_mpi
  - option/size headers: DIAGNOSTICS_SIZE.h, PTRACERS_OPTIONS.h

## Input variants (input*/)
- **input**: data.pkg on: usePTRACERS, useDiagnostics, useMNC
  - data: deltaTmom=1200.0, deltaTtracer=1200.0, nTimeSteps=4, nIter0=259200, eosType='LINEAR', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.FALSE., usingCartesianGrid=.TRUE.
  - namelist files: data data.diagnostics data.mnc data.pkg data.ptracers eedata

## Reference results
`output.txt`

Run: `cd verification; ./testreport -of <optfile> -t tutorial_advection_in_gyre` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
