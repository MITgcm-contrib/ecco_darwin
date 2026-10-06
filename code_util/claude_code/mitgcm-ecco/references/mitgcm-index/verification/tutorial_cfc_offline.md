# verification/tutorial_cfc_offline

## README (first 40 lines)
```
Tutorial Example: "Offline CFC Experiments"
============================================================
(formerly "cfc_offline" verification )

Configure and compile the code:
```
  cd build
  ../../../tools/genmake2 -mods ../code [-of my_platform_optionFile]
  make depend
  make
  cd ..
```

To run:
```
  cd run
  ln -s ../input/* .
  ./prepare_run
  ../build/mitgcmuv > output.txt
```

There is comparison output in the directory:
  results/output.txt

Comments:
  The input data is real*4

And to run the simpler (no CFC) offline test:
```
  cd run ; rm -f *
  ln -s ../input_tutorial/* .
  ln -s ../input/* .
  ./prepare_run
  ../build/mitgcmuv > output.tut
```
```

## Build variants (code*/)
- **code**: packages.conf = `gfd -mom_common -mom_fluxform -mom_vecinv gmredi offline ptracers gchem cfc`
  - expanded: generic_advdiff debug mdsio rw monitor gmredi offline ptracers gchem cfc
  - SIZE.h: grid 128x64x15; sNx=64, sNy=32, OLx=4, OLy=4, nSx=2, nSy=2, nPx=1, nPy=1, Nr=15; has SIZE.h_mpi
  - option/size headers: GMREDI_OPTIONS.h, PTRACERS_SIZE.h

## Input variants (input*/)
- **input**: data.pkg on: useGMRedi, usePTRACERS, useGCHEM, useOffLine
  - data: deltaTmom=900.0, deltaTtracer=43200.0, nTimeSteps=4, nIter0=4269600, eosType='POLY3', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE.
  - namelist files: data data.cfc data.gchem data.gmredi data.off data.pkg data.ptracers eedata prepare_run
- **input_tutorial**: data.pkg on: useGMRedi, usePTRACERS, useGCHEM, useOffLine
  - data: deltaTtracer=43200.0, nTimeSteps=4, nIter0=4269600, usingSphericalPolarGrid=.TRUE.
  - namelist files: data data.gchem data.gmredi data.off data.pkg data.ptracers eedata

## Reference results
`output.txt`

Run: `cd verification; ./testreport -of <optfile> -t tutorial_cfc_offline` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
