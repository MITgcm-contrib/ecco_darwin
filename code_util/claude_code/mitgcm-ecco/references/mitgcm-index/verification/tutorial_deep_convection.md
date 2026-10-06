# verification/tutorial_deep_convection

## README (first 40 lines)
```
Tutorial Example: "Surface Driven (Deep) Convection"
====================================================
(formerly "exp5" ;
 also "nonhydrostatic_deep_convection" in release.1 branch)

Configure and compile the code:
  cd build
  ../../../tools/genmake2 -mods ../code [-of my_platform_optionFile]
  make depend
  make
  cd ..

To run:
  cd run
  ln -s ../input/* .
  ln -s ../build/mitgcmuv .
  ./mitgcmuv > output.txt
  cd ..

There is comparison output in the directory:
  results/output.txt

Comments:
  The input data is real*8 and generated using the MATLAB script
  input/gendata.m

```

## Build variants (code*/)
- **code**: packages.conf = `gfd diagnostics`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor diagnostics
  - SIZE.h: grid 100x100x50; sNx=50, sNy=50, OLx=2, OLy=2, nSx=2, nSy=2, nPx=1, nPy=1, Nr=50; has SIZE.h_mpi
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h, MOM_COMMON_OPTIONS.h

## Input variants (input*/)
- **input**: data.pkg on: useDiagnostics
  - data: deltaT=20., nTimeSteps=3, nIter0=0, eosType='LINEAR', implicitFreeSurface=.TRUE., nonHydrostatic=.TRUE., usingCartesianGrid=.TRUE.
  - namelist files: data data.diagnostics data.pkg eedata
- **input.smag3d**: data.pkg on: useDiagnostics
  - data: deltaT=20., nTimeSteps=3, nIter0=0, eosType='LINEAR', implicitFreeSurface=.TRUE., nonHydrostatic=.TRUE., usingCartesianGrid=.TRUE., tempAdvScheme=77
  - namelist files: data data.pkg eedata

## Reference results
`output.smag3d.txt` `output.txt`

Run: `cd verification; ./testreport -of <optfile> -t tutorial_deep_convection` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
