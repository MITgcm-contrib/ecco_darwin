# verification/global_ocean.cs32x15

## README (first 40 lines)
```
global ocean using the cubed-sphere grid 32x32x32 with 15 levels
=================================================================

Specific option:
* Use Non-Linear Free surface formulation with z* coordinate
   with real fresh-water flux.
* Oceanic set-up on the cubed-sphere grid using the vector-invariant
   formulation.

Forcing :
 use Monthly mean climatological forcing (except P-E-R, annual mean).
 same data set as global-ocean lat-long experiments but interpolated
  on CS-32 grid.

Comments:
* bathymetry :
 designed to be coupled to Atmospheric model, therefore includes
 most of the semi-enclosed sea (Mediterranean, Black-Sea, Red-Sea,
   Hudson Bay ...)
 bathy_cs32.bin: initial bathymetry
   h < 0 is meant to stay wet-point whatever delZ(1) is ; Consequently
   the global ocean area is not affected by the vertical resolution.
 bathy_Hmin50.bin: bathymetry file used in the current set-up
    generated from bathy_cs32.bin using matlab script mk_bathy4gcm.m
 mk_bathy4gcm.m matlab script that deepen all shallow point up to 50m.
* global integral of E-P-R and annual mean net Q flux are zero.
* package thSIce and bulk_forc are included but not used in the standard
  set-up.

* additional forcing fields and parameter files are provided (in input.thsice)
  in order to illustrate the use of thSIce pkg.
  the output of a short run (20.iter) is given in results/output_thsice.txt

October 1rst, 2005:
* input.viscA4/data has been added to test biharmonic viscosity on CS-grid
  with side-drag. However, this set of parameters has only be used for
  short tests and is not recommended to begin with.
```

## Build variants (code*/)
- **code**: packages.conf = `exch2 gfd gmredi ggl90 bulk_force exf -cal seaice thsice diagnostics mnc`
  - expanded: exch2 mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor gmredi ggl90 bulk_force exf seaice thsice diagnostics mnc
  - SIZE.h: grid 384x16x15; sNx=32, sNy=16, OLx=4, OLy=4, nSx=12, nSy=1, nPx=1, nPy=1, Nr=15; has SIZE.h_mpi
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h, DIAG_OPTIONS.h, EXF_OPTIONS.h, GGL90_OPTIONS.h, SEAICE_OPTIONS.h
- **code_ad**: packages.conf = `exch2 gfd -mom_fluxform gmredi exf seaice thsice diagnostics adjoint`
  - expanded: exch2 gfd -mom_fluxform gmredi exf seaice thsice diagnostics autodiff cost ctrl grdchk
  - SIZE.h: grid 384x16x15; sNx=32, sNy=16, OLx=4, OLy=4, nSx=12, nSy=1, nPx=1, nPy=1, Nr=15; has SIZE.h_mpi
  - option/size headers: AUTODIFF_OPTIONS.h, COST_OPTIONS.h, CPP_OPTIONS.h, CTRL_OPTIONS.h, CTRL_SIZE.h, DIAGNOSTICS_SIZE.h, DIAG_OPTIONS.h, EXF_OPTIONS.h, GMREDI_OPTIONS.h, SEAICE_OPTIONS.h, THSICE_OPTIONS.h
  - modified/extra source: cost_test.F, tamc.h
- **code_alt**: packages.conf = `(inherits/none)`
  - modified/extra source: code.192t_8x4
- **code_tap**: packages.conf = `exch2 gfd -mom_fluxform gmredi exf seaice thsice diagnostics tapenade adjoint`
  - expanded: exch2 gfd -mom_fluxform gmredi exf seaice thsice diagnostics tapenade autodiff cost ctrl grdchk
  - SIZE.h: grid 384x16x15; sNx=32, sNy=16, OLx=4, OLy=4, nSx=12, nSy=1, nPx=1, nPy=1, Nr=15; has SIZE.h_mpi
  - option/size headers: AUTODIFF_OPTIONS.h, COST_OPTIONS.h, CPP_OPTIONS.h, CTRL_OPTIONS.h, CTRL_SIZE.h, DIAGNOSTICS_SIZE.h, DIAG_OPTIONS.h, EXF_OPTIONS.h, GMREDI_OPTIONS.h, SEAICE_OPTIONS.h, THSICE_OPTIONS.h
  - modified/extra source: cost_test.F

## Input variants (input*/)
- **input**: data.pkg on: useGMRedi, useDiagnostics, useMNC
  - data: deltaTmom=1200., deltaTtracer=86400., nTimeSteps=20, nIter0=72000, nonlinFreeSurf=4, select_rStar=2, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingCurvilinearGrid=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.diagnostics data.gmredi data.mnc data.pkg eedata prepare_run
- **input.icedyn**: data.pkg on: useGMRedi, useEXF, useSEAICE, useThSIce, useDiagnostics
  - data: deltaTmom=1200., deltaTtracer=86400., nTimeSteps=10, nIter0=36000, nonlinFreeSurf=4, select_rStar=2, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingCurvilinearGrid=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.diagnostics data.exf data.gmredi data.ice data.pkg data.seaice
- **input.in_p**: data.pkg on: useEXF, useSEAICE, useGGL90, useDiagnostics
  - data: deltaTmom=1200., deltaTtracer=86400., nTimeSteps=10, nIter0=0, nonlinFreeSurf=4, select_rStar=2, eosType='TEOS10', implicitFreeSurface=.TRUE., usingCurvilinearGrid=.TRUE., buoyancyRelation='OCEANICP', useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.diagnostics data.ggl90 data.pkg prepare_run
- **input.seaice**: data.pkg on: useGMRedi, useEXF, useSEAICE, useDiagnostics
  - data: deltaTmom=1200., deltaTtracer=86400., nTimeSteps=10, nIter0=36000, nonlinFreeSurf=4, select_rStar=2, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingCurvilinearGrid=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.diagnostics data.exf data.gmredi data.pkg data.seaice prepare_run
- **input.thsice**: data.pkg on: useGMRedi, useBulkforce, useThSIce, useDiagnostics
  - data: deltaTmom=1200., deltaTtracer=86400., nTimeSteps=20, nIter0=36000, nonlinFreeSurf=4, select_rStar=2, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingCurvilinearGrid=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.blk data.diagnostics data.gmredi data.ice data.pkg
- **input.viscA4**: data.pkg on: useGMRedi, useDiagnostics
  - data: nTimeSteps=10, nIter0=86400, eosType='JMD95Z', implicitFreeSurface=.TRUE., nonHydrostatic=.TRUE., usingCurvilinearGrid=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.diagnostics data.gmredi data.pkg prepare_run
- **input_ad**: data.pkg on: useGMRedi, useDiagnostics, useGrdchk
  - data: deltaTmom=1200., deltaTtracer=86400., nTimeSteps=5, nIter0=72000, nonlinFreeSurf=4, select_rStar=2, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingCurvilinearGrid=.TRUE., tempAdvScheme=30, saltAdvScheme=30, useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.autodiff data.cost data.ctrl data.diagnostics data.gmredi data.grdchk data.optim data.pkg eedata prepare_run
- **input_ad.seaice**: data.pkg on: useGMRedi, useEXF, useSEAICE, useDiagnostics, useGrdchk
  - data: deltaTmom=1200., deltaTtracer=86400., nTimeSteps=2, nIter0=36000, nonlinFreeSurf=4, select_rStar=2, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingCurvilinearGrid=.TRUE., tempAdvScheme=33, saltAdvScheme=33, useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.autodiff data.cal data.diagnostics data.exf data.pkg data.seaice prepare_run
- **input_ad.seaice_dynmix**: data.pkg on: useGMRedi, useEXF, useSEAICE, useGrdchk
  - data: deltaTmom=1200., deltaTtracer=86400., nTimeSteps=5, nIter0=36000, nonlinFreeSurf=2, select_rStar=1, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingCurvilinearGrid=.TRUE., tempAdvScheme=30, saltAdvScheme=30, useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.autodiff data.exf data.pkg data.seaice prepare_run
- **input_ad.thsice**: data.pkg on: useGMRedi, useEXF, useThSIce, useGrdchk
  - data: deltaTmom=1200., deltaTtracer=86400., nTimeSteps=5, nIter0=36000, nonlinFreeSurf=2, select_rStar=1, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingCurvilinearGrid=.TRUE., tempAdvScheme=30, saltAdvScheme=30, useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.ctrl data.exf data.gmredi data.grdchk data.ice data.pkg prepare_run
- **input_tap**: data.pkg on: -
  - namelist files: prepare_run

## Reference results
`output.icedyn.txt` `output.in_p.txt` `output.seaice.txt` `output.thsice.txt` `output.txt` `output.viscA4.txt` `output_adm.seaice.txt` `output_adm.seaice_dynmix.txt` `output_adm.thsice.txt` `output_adm.txt` `output_tap_adj.txt` `output_tap_tlm.txt` `output_tlm.seaice.txt.gz` `output_tlm.seaice_dynmix.txt.gz` `output_tlm.thsice.txt.gz` `output_tlm.txt.gz`

Run: `cd verification; ./testreport -of <optfile> -t global_ocean.cs32x15` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
