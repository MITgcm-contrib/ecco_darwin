# verification/adjustment.cs-32x32x1

## README (first 40 lines)
```
Simple 1 layer, Barotropic adjustment on the Sphere, using the
                cubed-sphere grid 32x32x32
Contains also a "minimal" test case (just compile eesupp/src + pkgs)
                that does not do much.
=================================================================

General Description:
* using the same executable, 2 set-up can be tested, corresponding
  to input dir. "input" and "input.nlfs".
* Default set-up (input & parameter files in dir "input"):
  Oceanic set-up initially at rest, with flat bottom and a large
  quasi-rectangular continent. An initial free-surface large-scale
  anomaly centered at the equator triggers a barotropic adjustment
  and generates External Inertial-Gravity waves (Poincare waves).
  Use linear Free-Surface and linear dynamics (no momentum advection)
* Additional set-up (input & parameter files in dir "input.nlfs"):
  Atmospheric set-up, without orography, initially at rest.
  An initial large-scale surface pressure anomaly generated pure
  external gravity waves.
  Use non-linear Free-Surface, linear dynamics (no momentum advection)
  without rotation.

IMPORTANT: For the purpose of testing multiple tiles and "blank-tiles":
* Use multiple tiles (8) per cube-face (tile size: 16x8),
  which results in a total of 48 tiles:
    code/SIZE.h
* The oceanic-set-up contains 4 empty tiles (tiles: 11,12,13,14)
  associated with the large continent.
  This gives the opportunity to test the "blank-tiles" option of
  the EXCH2 pkg. An MPI version of this set-up is available to
  test this "blank-tiles" option:
    code/SIZE.h_mpi          : to replace code/SIZE.h
    input/data.exch2.mpi     : to be rename to data.exch2
    code/CPP_EEOPTIONS.h_mpi : to replace eesupp/inc/CPP_EEOPTIONS.h
  However, this particular (MPI) executable cannot be used for the
  atmospheric set-up (input.nlfs) and, in this case, an error
  from S/R EXCH2_CHECK_DEPTHS will stop the execution.

Forcing: none
Input Files (initial conditions):
... (truncated)
```

## Build variants (code*/)
- **code**: packages.conf = `exch2 gfd -generic_advdiff diagnostics`
  - expanded: exch2 mom_common mom_fluxform mom_vecinv debug mdsio rw monitor diagnostics
  - SIZE.h: grid 768x8x1; sNx=16, sNy=8, OLx=2, OLy=2, nSx=48, nSy=1, nPx=1, nPy=1, Nr=1; has SIZE.h_mpi
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h, W2_OPTIONS.h
  - modified/extra source: CPP_EEOPTIONS.h_mpi
- **code_min**: packages.conf = `exch2 debug`
  - expanded: exch2 debug
  - SIZE.h: grid 768x8x1; sNx=16, sNy=8, OLx=2, OLy=2, nSx=48, nSy=1, nPx=1, nPy=1, Nr=1; has SIZE.h_mpi
  - option/size headers: CPP_OPTIONS.h
  - modified/extra source: main.F

## Input variants (input*/)
- **input**: data.pkg on: useDiagnostics
  - data: deltaT=900., nTimeSteps=24, nIter0=0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCurvilinearGrid=.TRUE., buoyancyRelation='OCEANIC'
  - namelist files: data data.diagnostics data.exch2.mpi data.pkg eedata prepare_run
- **input.nlfs**: data.pkg on: -
  - data: deltaT=180.0, nTimeSteps=20, nIter0=0, nonlinFreeSurf=3, eosType='IDEALG', implicitFreeSurface=.TRUE., usingCurvilinearGrid=.TRUE., buoyancyRelation='ATMOSPHERIC'
  - namelist files: data data.pkg
- **input_min**: data.pkg on: -
  - namelist files: data.exch2.mpi eedata

## Reference results
`output.nlfs.txt` `output.txt`

Run: `cd verification; ./testreport -of <optfile> -t adjustment.cs-32x32x1` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
