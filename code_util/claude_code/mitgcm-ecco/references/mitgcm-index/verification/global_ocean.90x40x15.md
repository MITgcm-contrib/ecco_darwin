# verification/global_ocean.90x40x15

## README (first 40 lines)
```
============================================
Example: "4x4 Global Simulation with Seasonal Forcing"
============================================
(see also similar set-up in: verification/tutorial_global_oce_latlon/)

From verification/global_ocean.90x40x15 dir:

Configure and compile the code:
  cd build
  ../../../tools/genmake2 -mods ../code [-of my_platform_optionFile]
 [make Clean]
  make depend
  make
  cd ..

To run:
  cd run
  ln -s ../input/* .
  ./prepare_run
  ln -s ../build/mitgcmuv .
  ./mitgcmuv > output.txt
  cd ..

There is comparison output in the directory:
  results/output.txt

There is comparison output in directory:
  (verification/global_ocean.90x40x15/) results

Comments:
o The input data is real*4.
o The surface fluxes are derived from monthly means of the NCEP climatology;
  - a matlab script is provided that created the surface flux data files from
    the original NCEP data: ncep2global_ocean.m in the diags_matlab directory,
    needs editing to adjust search paths.
o matlab scripts that make a simple diagnostic (barotropic stream function,
  overturning stream functions, averaged hydrography etc.) is provided in
  verification/tutorial_global_oce_latlon/diags_matlab:
  - mit_loadglobal is the toplevel script that run all other scripts
  - mit_globalmovie animates theta, salinity, and 3D-velocity field for
... (truncated)
```

## Build variants (code*/)
- **code**: packages.conf = `exch2 oceanic -kpp cd_code down_slope ggl90 ptracers sbo diagnostics`
  - expanded: exch2 mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor gmredi cd_code down_slope ggl90 ptracers sbo diagnostics
  - SIZE.h: grid 90x40x15; sNx=10, sNy=10, OLx=3, OLy=3, nSx=9, nSy=4, nPx=1, nPy=1, Nr=15; has SIZE.h_mpi
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h, GAD_OPTIONS.h, GGL90_OPTIONS.h
- **code_ad**: packages.conf = `gfd cd_code gmredi sbo mnc autodiff cost ctrl grdchk`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor cd_code gmredi sbo mnc autodiff cost ctrl grdchk
  - SIZE.h: grid 90x40x15; sNx=45, sNy=20, OLx=3, OLy=3, nSx=2, nSy=2, nPx=1, nPy=1, Nr=15; has SIZE.h_mpi
  - option/size headers: COST_OPTIONS.h, CPP_OPTIONS.h, CTRL_OPTIONS.h, GMREDI_OPTIONS.h
  - modified/extra source: CPP_EEOPTIONS.h_mpi, tamc.h

## Input variants (input*/)
- **input**: data.pkg on: useGMRedi, useSBO, useDiagnostics
  - data: deltaTmom=1800., deltaTtracer=86400., nTimeSteps=10, nIter0=36000, nonlinFreeSurf=4, select_rStar=2, eosType='JMD95P', usingSphericalPolarGrid=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.diagnostics data.exch2.mpi data.gmredi data.pkg data.sbo eedata prepare_run
- **input.dwnslp**: data.pkg on: useGMRedi, useDOWN_SLOPE, usePTRACERS, useDiagnostics
  - data: deltaTmom=1800., deltaTtracer=86400., nTimeSteps=10, nIter0=36000, eosType='JMD95P', usingSphericalPolarGrid=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.diagnostics data.down_slope data.exch2.mpi data.gmredi data.pkg data.ptracers eedata
- **input.idemix**: data.pkg on: useGMRedi, useDiagnostics, useGGL90
  - data: deltaTmom=1800., deltaTtracer=86400., nTimeSteps=10, nIter0=0, eosType='JMD95P', usingSphericalPolarGrid=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.diagnostics data.exch2.mpi data.ggl90 data.gmredi data.pkg data.ptracers eedata
- **input_ad**: data.pkg on: useGMRedi, useGrdchk
  - data: deltaTmom=1200.0, deltaTtracer=43200.0, nTimeSteps=10, nIter0=0, nonlinFreeSurf=2, eosType='JMD95Z', usingSphericalPolarGrid=.TRUE., tempAdvScheme=30, saltAdvScheme=30, useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.autodiff data.cost data.ctrl data.gmredi data.grdchk data.optim data.pkg eedata prepare_run
- **input_ad.bottomdrag**: data.pkg on: useGMRedi, useGrdchk, useMNC
  - data: deltaTmom=1200.0, deltaTtracer=43200.0, nTimeSteps=10, nIter0=0, nonlinFreeSurf=4, eosType='JMD95Z', usingSphericalPolarGrid=.TRUE., tempAdvScheme=30, saltAdvScheme=30, useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.autodiff data.ctrl data.grdchk data.mnc data.pkg
- **input_ad.kapgm**: data.pkg on: -
  - namelist files: data.ctrl data.grdchk
- **input_ad.kapredi**: data.pkg on: -
  - namelist files: data.ctrl data.grdchk

## Reference results
`output.dwnslp.txt` `output.idemix.txt` `output.txt` `output_adm.bottomdrag.txt` `output_adm.kapgm.txt` `output_adm.kapredi.txt` `output_adm.txt` `output_tlm.bottomdrag.txt.gz` `output_tlm.kapgm.txt.gz` `output_tlm.kapredi.txt.gz` `output_tlm.txt.gz`

Run: `cd verification; ./testreport -of <optfile> -t global_ocean.90x40x15` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
