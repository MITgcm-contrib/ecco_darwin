# verification/dome

## Build variants (code*/)
- **code**: packages.conf = `gfd obcs diagnostics`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor obcs diagnostics
  - SIZE.h: grid 200x45x25; sNx=25, sNy=15, OLx=3, OLy=3, nSx=8, nSy=3, nPx=1, nPy=1, Nr=25; has SIZE.h_mpi
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h
  - modified/extra source: obcs_calc.F

## Input variants (input*/)
- **input**: data.pkg on: useOBCS, useDiagnostics
  - data: deltaT=300., nTimeSteps=20, nIter0=0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., tempAdvScheme=77, saltAdvScheme=77
  - namelist files: data data.diagnostics data.obcs data.pkg eedata

## Reference results
`output.txt`

Run: `cd verification; ./testreport -of <optfile> -t dome` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
