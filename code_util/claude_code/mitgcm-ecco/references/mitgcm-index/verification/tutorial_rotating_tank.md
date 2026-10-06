# verification/tutorial_rotating_tank

## Build variants (code*/)
- **code**: packages.conf = `gfd diagnostics mnc`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor diagnostics mnc
  - SIZE.h: grid 120x23x29; sNx=30, sNy=23, OLx=3, OLy=3, nSx=4, nSy=1, nPx=1, nPy=1, Nr=29; has SIZE.h_mpi
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h
  - modified/extra source: apply_forcing.F

## Input variants (input*/)
- **input**: data.pkg on: useMNC
  - data: deltaT=0.1, nTimeSteps=20, nIter0=0, eosType='LINEAR', implicitFreeSurface=.FALSE., nonHydrostatic=.TRUE., usingCylindricalGrid=.TRUE.
  - namelist files: data data.mnc data.pkg eedata

## Reference results
`output.txt`

Run: `cd verification; ./testreport -of <optfile> -t tutorial_rotating_tank` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
