# verification/advect_xy

## Build variants (code*/)
- **code**: packages.conf = `(inherits/none)`
  - SIZE.h: grid 20x20x1; sNx=20, sNy=10, OLx=3, OLy=3, nSx=1, nSy=2, nPx=1, nPy=1, Nr=1; has SIZE.h_mpi
  - option/size headers: CPP_OPTIONS.h, GAD_OPTIONS.h
  - modified/extra source: ini_salt.F, ini_theta.F, ini_vel.F

## Input variants (input*/)
- **input**: data.pkg on: -
  - data: deltaT=2500.0, endTime=200000., startTime=0, implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., tempAdvScheme=80, saltAdvScheme=33, momStepping=.FALSE.
  - namelist files: data data.pkg eedata
- **input.ab3_c4**: data.pkg on: -
  - data: deltaT=2750., endTime=275000., startTime=0, implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., tempAdvScheme=4, saltAdvScheme=4, momStepping=.FALSE.
  - namelist files: data data.pkg eedata

## Reference results
`output.ab3_c4.txt` `output.txt`

Run: `cd verification; ./testreport -of <optfile> -t advect_xy` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
