# verification/tutorial_reentrant_channel

## README (first 40 lines)
```
Tutorial Example: "Reentrant channel"
(Southern Ocean Reentrant Channel Example)
==========================================

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
  ln -s ../build/mitgcmuv .
  ./mitgcmuv > output.txt
  cd ..
```

There is comparison output in the directory:
results/output.txt

Comments:
  The input data is real*4 and generated using the MATLAB script gendata_50km.m.
```

## Build variants (code*/)
- **code**: packages.conf = `gfd gmredi rbcs layers diagnostics`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor gmredi rbcs layers diagnostics
  - SIZE.h: grid 20x40x49; sNx=20, sNy=10, OLx=4, OLy=4, nSx=1, nSy=4, nPx=1, nPy=1, Nr=49; has SIZE.h_mpi
  - option/size headers: DIAGNOSTICS_SIZE.h, LAYERS_SIZE.h
  - modified/extra source: SIZE.h_eddy

## Input variants (input*/)
- **input**: data.pkg on: useGMRedi, useRBCS, useLayers, useDiagnostics
  - data: deltaT=1000.0, nTimeSteps=10, nIter0=0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., tempAdvScheme=7
  - namelist files: data data.diagnostics data.gmredi data.layers data.pkg data.rbcs eedata

## Reference results
`output.txt`

Run: `cd verification; ./testreport -of <optfile> -t tutorial_reentrant_channel` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
