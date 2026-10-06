# verification/exp2

## README (first 40 lines)
```
Example: "4x4 Steady Global Simulation"
=======================================

The notes for this experiment are now on-line.  Please see the third
chapter ("Getting started with MITgcm") which is available at:

  http://mitgcm.org/public/r2_manual/latest/online_documents/node1.html

```

## Build variants (code*/)
- **code**: packages.conf = `gfd cd_code`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor cd_code
  - SIZE.h: grid 90x40x20; sNx=45, sNy=20, OLx=2, OLy=2, nSx=2, nSy=2, nPx=1, nPy=1, Nr=20; has SIZE.h_mpi
  - option/size headers: CD_CODE_OPTIONS.h, CPP_OPTIONS.h

## Input variants (input*/)
- **input**: data.pkg on: -
  - data: deltaTmom=2400.0, deltaTtracer=108000.0, endTime=2808000., startTime=0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE.
  - namelist files: data data.pkg eedata
- **input.rigidLid**: data.pkg on: -
  - data: deltaTmom=2400.0, deltaTtracer=108000.0, nTimeSteps=12, startTime=0, eosType='LINEAR', usingSphericalPolarGrid=.TRUE.
  - namelist files: data

## Reference results
`output.rigidLid.txt` `output.txt`

Run: `cd verification; ./testreport -of <optfile> -t exp2` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
