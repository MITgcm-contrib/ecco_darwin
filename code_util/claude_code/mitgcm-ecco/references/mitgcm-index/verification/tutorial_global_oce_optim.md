# verification/tutorial_global_oce_optim

## README (first 40 lines)
```
Tutorial Example: "Global Ocean State Estimation"
(Global Ocean State Estimation at 4.o Resolution)
=================================================

Configure and compile the code:
  cd build
  ../../../tools/genmake2 -mods ../code_ad [-of my_platform_optionFile]
  make depend
  make adall
  cd ..

To run:
  cd run
  ln -s ../input_ad/* .
  ./prepare_run
  ln -s ../build/mitgcmuv_ad .
  ./mitgcmuv_ad > output_adm.txt
  cd ..

There is comparison output in the directory:
  results/output_adm.txt

grep for grdchk output:
  grep 'grdchk output' output_adm.txt

Comments:
```

## Build variants (code*/)
- **code_ad**: packages.conf = `gfd cd_code gmredi autodiff cost ctrl grdchk`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor cd_code gmredi autodiff cost ctrl grdchk
  - SIZE.h: grid 90x40x15; sNx=45, sNy=20, OLx=2, OLy=2, nSx=2, nSy=2, nPx=1, nPy=1, Nr=15; has SIZE.h_mpi
  - option/size headers: AUTODIFF_OPTIONS.h, COST_OPTIONS.h, CPP_OPTIONS.h, CTRL_OPTIONS.h, GMREDI_OPTIONS.h, MOM_COMMON_OPTIONS.h
  - modified/extra source: code_ad_diff.list, cost_hflux.F, cost_local.h, cost_temp.F, cost_weights.F, tamc.h

## Input variants (input*/)
- **input_ad**: data.pkg on: useGMRedi, useGrdchk
  - data: deltaTmom=1800., deltaTtracer=86400., nTimeSteps=10, nIter0=0, eosType='JMD95Z', usingSphericalPolarGrid=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.autodiff data.cost data.ctrl data.gmredi data.grdchk data.optim data.pkg eedata prepare_run

## Reference results
`output_adm.txt`

Run: `cd verification; ./testreport -of <optfile> -t tutorial_global_oce_optim` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
