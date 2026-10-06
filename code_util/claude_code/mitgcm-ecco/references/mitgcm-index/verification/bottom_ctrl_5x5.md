# verification/bottom_ctrl_5x5

## Build variants (code*/)
- **code_ad**: packages.conf = `gfd cd_code adjoint`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor cd_code autodiff cost ctrl grdchk
  - SIZE.h: grid 5x5x4; sNx=5, sNy=5, OLx=2, OLy=2, nSx=1, nSy=1, nPx=1, nPy=1, Nr=4
  - option/size headers: AUTODIFF_OPTIONS.h, COST_OPTIONS.h, CPP_OPTIONS.h, CTRL_OPTIONS.h, GAD_OPTIONS.h
  - modified/extra source: cost_test.F, dummy_in_hfac.F, tamc.h

## Input variants (input*/)
- **input_ad**: data.pkg on: useGrdchk
  - data: deltaTmom=3600.0, deltaTtracer=3600.0, nTimeSteps=100, nIter0=0, nonlinFreeSurf=0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.FALSE., usingCartesianGrid=.TRUE.
  - namelist files: data data.autodiff data.cost data.ctrl data.grdchk data.optim data.pkg eedata prepare_run
- **input_ad.facg2d**: data.pkg on: -
  - data: deltaTmom=3600.0, deltaTtracer=3600.0, nTimeSteps=100, nIter0=0, nonlinFreeSurf=0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.FALSE., usingCartesianGrid=.TRUE.
  - namelist files: data data.autodiff

## Reference results
`output_adm.facg2d.txt` `output_adm.txt` `output_tlm.facg2d.txt.gz` `output_tlm.txt.gz`

Run: `cd verification; ./testreport -of <optfile> -t bottom_ctrl_5x5` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
