# verification/advect_xz

## Build variants (code*/)
- **code**: packages.conf = `gfd -mom_common -mom_fluxform -mom_vecinv diagnostics`
  - expanded: generic_advdiff debug mdsio rw monitor diagnostics
  - SIZE.h: grid 20x1x20; sNx=10, sNy=1, OLx=4, OLy=4, nSx=2, nSy=1, nPx=1, nPy=1, Nr=20; has SIZE.h_mpi
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h, GAD_OPTIONS.h

## Input variants (input*/)
- **input**: data.pkg on: -
  - data: deltaT=1200., endTime=240000., startTime=0., implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., tempAdvScheme=42, saltAdvScheme=81, momStepping=.FALSE.
  - namelist files: data data.pkg eedata
- **input.nlfs**: data.pkg on: useDiagnostics
  - data: deltaT=1200., endTime=240000., startTime=0., nonlinFreeSurf=4, select_rStar=2, implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., tempAdvScheme=77, saltAdvScheme=3, momStepping=.FALSE.
  - namelist files: data data.diagnostics data.pkg eedata
- **input.pqm**: data.pkg on: -
  - data: deltaT=1200., endTime=240000., startTime=0., implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., tempAdvScheme=51, saltAdvScheme=52, momStepping=.FALSE.
  - namelist files: data eedata

## Reference results
`output.nlfs.txt` `output.pqm.txt` `output.txt`

Run: `cd verification; ./testreport -of <optfile> -t advect_xz` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
