# verification/matrix_example

## README (first 40 lines)
```
Example to test pkg/matrix

Instructions for building:

1) mkdir build
2) cd build
3) run genmake with: -mods=../code and -enable=matrix
e.g.,
../../../tools/genmake2  "-mods=../code" "-optfile=../../../tools/build_options/darwin_ppc_xlf_spk" "-enable=matrix"

4) make depend
5) make

Instructions for running:

1) mkdir run
2) cd run
3) cp -p ../input/* .
4) ./mitgcmuv > output.txt

Results:

compare with results/output.txt

```

## Build variants (code*/)
- **code**: packages.conf = `gfd matrix`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor matrix
  - SIZE.h: grid 32x32x1; sNx=16, sNy=8, OLx=3, OLy=3, nSx=2, nSy=4, nPx=1, nPy=1, Nr=1; has SIZE.h_mpi
  - option/size headers: PTRACERS_SIZE.h

## Input variants (input*/)
- **input**: data.pkg on: usePTRACERS, useMATRIX
  - data: deltaTmom=20000.0, deltaTtracer=20000.0, nTimeSteps=10, nIter0=200000, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE.
  - namelist files: data data.matrix data.pkg data.ptracers eedata

## Reference results
`output.txt`

Run: `cd verification; ./testreport -of <optfile> -t matrix_example` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
