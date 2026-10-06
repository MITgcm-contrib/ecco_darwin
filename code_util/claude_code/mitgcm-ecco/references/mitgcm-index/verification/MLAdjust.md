# verification/MLAdjust

## README (first 40 lines)
```
Simple set-up to test flow-dependent horizontal viscosity implementation.

 Domain size is 50 x 26 x 40 grid-cells,
 with uniform resolution dx=dy= 1.km , dz = 5.m

 Zonally re-entrant, flat bottom channel (closed by Northern Wall @ j=26)

 input files are real*8 (see matlab script input/gendata.m )

 start from initial density field (given by initial Temp), no forcing

 Test Exp. |  Momentum   |  Viscosity  | useFullLeith | Biharmonic | side-drag   |
           | formulation | formulation |              | vs harmonic| (no_slip BC)|
----------------------------------------------------------------------------------
 standard  | Vector-Inv. |  Vort-Div   |  FullLeith   |   viscC4   |    No       |
(dir=input)|             |             |              |            |             |
----------------------------------------------------------------------------------
 .A4FlxF   |  Flux-Form  |  FLux-Form  |  FullLeith   |   viscC4   |    Yes      |
----------------------------------------------------------------------------------
 .AhFlxF   |  Flux-Form  |  Flux-Form  |    No        |   viscC2   |    No       |
----------------------------------------------------------------------------------
 .AhVrDv   | Vector-Inv. |  Vort-Div   |  FullLeith   |   viscC2   |    Yes      |
----------------------------------------------------------------------------------
 .AhStTn   | Vector-Inv. | Strain-Tens |  FullLeith   |   viscC2   |    Yes      |
----------------------------------------------------------------------------------
 .QGLeith  | Vector-Inv. |  Vort-Div   |  FullLeith   |   viscC2   |    Yes      |
----------------------------------------------------------------------------------
 .QGLthGM  | Vector-Inv. |  Vort-Div   |  FullLeith   |   viscC2   |    Yes      |
----------------------------------------------------------------------------------

Notes:
1) Stain-Tension viscosity formulation is used when setting
     useStrainTensionVisc=.TRUE.,
   and currently only available with Vector-Invariant momentum.
   Default is .False., to use Vorticity & Divergence formulation
   (pkg/mom_vecinv) or Flux-Form formulation (pkg/mom_fluxform)
2) test .A4FlxF (input.A4FlxF) starts from a pickup-file,
   all other test exp. start from iter0=0
3) test with Biharmonic visc. are using:
    viscC4leith = viscC4leithD = 1.85,
... (truncated)
```

## Build variants (code*/)
- **code**: packages.conf = `exch2 gfd gmredi flt diagnostics mnc`
  - expanded: exch2 mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor gmredi flt diagnostics mnc
  - SIZE.h: grid 50x26x40; sNx=25, sNy=13, OLx=3, OLy=3, nSx=2, nSy=2, nPx=1, nPy=1, Nr=40; has SIZE.h_mpi
  - option/size headers: DIAGNOSTICS_SIZE.h, GMREDI_OPTIONS.h, MOM_COMMON_OPTIONS.h

## Input variants (input*/)
- **input**: data.pkg on: useMNC, useDiagnostics
  - data: deltaT=1200., nTimeSteps=12, nIter0=0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., tempAdvScheme=33, saltAdvScheme=33
  - namelist files: data data.diagnostics data.mnc data.pkg eedata
- **input.A4FlxF**: data.pkg on: -
  - data: deltaT=1200., nTimeSteps=12, nIter0=36, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., tempAdvScheme=33, saltAdvScheme=33
  - namelist files: data eedata
- **input.AhFlxF**: data.pkg on: -
  - data: deltaT=1200., nTimeSteps=12, nIter0=0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., tempAdvScheme=33, saltAdvScheme=33
  - namelist files: data eedata
- **input.AhStTn**: data.pkg on: -
  - data: deltaT=1200., nTimeSteps=12, nIter0=0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., tempAdvScheme=33, saltAdvScheme=33
  - namelist files: data eedata
- **input.AhVrDv**: data.pkg on: -
  - data: deltaT=1200., nTimeSteps=12, nIter0=0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., tempAdvScheme=33, saltAdvScheme=33
  - namelist files: data eedata
- **input.QGLeith**: data.pkg on: useMNC, useDiagnostics
  - data: deltaT=1200., nTimeSteps=12, nIter0=0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., tempAdvScheme=33, saltAdvScheme=33
  - namelist files: data data.pkg eedata
- **input.QGLthGM**: data.pkg on: useGMRedi, useMNC, useDiagnostics
  - data: deltaT=1200., nTimeSteps=12, nIter0=0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., tempAdvScheme=33, saltAdvScheme=33
  - namelist files: data data.diagnostics data.gmredi data.pkg eedata

## Reference results
`output.A4FlxF.txt` `output.AhFlxF.txt` `output.AhStTn.txt` `output.AhVrDv.txt` `output.QGLeith.txt` `output.QGLthGM.txt` `output.txt`

Run: `cd verification; ./testreport -of <optfile> -t MLAdjust` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
