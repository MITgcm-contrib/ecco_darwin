# verification/advect_cs

## Build variants (code*/)
- **code**: packages.conf = `exch2 gfd diagnostics`
  - expanded: exch2 mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor diagnostics
  - SIZE.h: grid 64x96x1; sNx=32, sNy=32, OLx=4, OLy=4, nSx=2, nSy=3, nPx=1, nPy=1, Nr=1; has SIZE.h_mpi
  - option/size headers: DIAGNOSTICS_SIZE.h, GAD_OPTIONS.h
  - modified/extra source: ini_vel.F

## Input variants (input*/)
- **input**: data.pkg on: useDiagnostics
  - data: deltaT=2700., nIter0=0, endTime=518400., usingCurvilinearGrid=.TRUE., tempAdvScheme=33, saltAdvScheme=80, momStepping=.FALSE.
  - namelist files: data data.diagnostics data.pkg eedata prepare_run

## Reference results
`output.txt`

Run: `cd verification; ./testreport -of <optfile> -t advect_cs` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
