# verification/inverted_barometer

## README (first 40 lines)
```
Example: "4 Layer Double Gyre with Pressure Loading"
====================================================

Comments:
The input data is real*8 and can be generated with the supplied
matlab script gendata.m
To change the input precision to real*4 change readBinaryPrec=32 in
data as well as in gendata.m

The experiment follow roughly the analytical analysis
Wunsch and Stammer, Atmospheric loading and the oceanic "inverted
barometer" effect. Rev. Geophys., 35, pp. 79-107, 1997.

```

## Build variants (code*/)
- **code**: packages.conf = `gfd diagnostics mnc`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor diagnostics mnc
  - SIZE.h: grid 60x60x4; sNx=30, sNy=30, OLx=2, OLy=2, nSx=2, nSy=2, nPx=1, nPy=1, Nr=4; has SIZE.h_mpi
  - option/size headers: DIAGNOSTICS_SIZE.h

## Input variants (input*/)
- **input**: data.pkg on: useDiagnostics, useMNC
  - data: deltaTmom=1200.0, deltaTtracer=1200.0, endTime=48000., startTime=0., eosType='LINEAR', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.FALSE., usingCartesianGrid=.TRUE.
  - namelist files: data data.diagnostics data.mnc data.pkg eedata

## Reference results
`output.txt`

Run: `cd verification; ./testreport -of <optfile> -t inverted_barometer` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
