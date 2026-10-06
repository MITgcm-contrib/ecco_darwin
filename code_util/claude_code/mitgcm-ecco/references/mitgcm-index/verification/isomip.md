# verification/isomip

## Build variants (code*/)
- **code**: packages.conf = `gfd ggl90 cd_code obcs shelfice steep_icecavity icefront diagnostics mnc`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor ggl90 cd_code obcs shelfice steep_icecavity icefront diagnostics mnc
  - SIZE.h: grid 50x100x30; sNx=25, sNy=25, OLx=3, OLy=3, nSx=2, nSy=4, nPx=1, nPy=1, Nr=30; has SIZE.h_mpi
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h, DIAG_OPTIONS.h
- **code_ad**: packages.conf = `gfd ggl90 cd_code shelfice steep_icecavity diagnostics mnc autodiff cost ctrl grdchk`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor ggl90 cd_code shelfice steep_icecavity diagnostics mnc autodiff cost ctrl grdchk
  - SIZE.h: grid 50x100x30; sNx=25, sNy=25, OLx=4, OLy=4, nSx=2, nSy=4, nPx=1, nPy=1, Nr=30; has SIZE.h_mpi
  - option/size headers: COST_OPTIONS.h, CPP_OPTIONS.h, CTRL_OPTIONS.h, CTRL_SIZE.h, DIAGNOSTICS_SIZE.h, SHELFICE_OPTIONS.h
  - modified/extra source: cost_test.F, tamc.h
- **code_tap**: packages.conf = `gfd monitor cd_code gmredi shelfice tapenade adjoint`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor cd_code gmredi shelfice tapenade autodiff cost ctrl grdchk
  - SIZE.h: grid 50x100x30; sNx=25, sNy=25, OLx=3, OLy=3, nSx=2, nSy=4, nPx=1, nPy=1, Nr=30; has SIZE.h_mpi
  - option/size headers: AUTODIFF_OPTIONS.h, COST_OPTIONS.h, CPP_OPTIONS.h, CTRL_OPTIONS.h, CTRL_SIZE.h, DIAGNOSTICS_SIZE.h, GMREDI_OPTIONS.h, SHELFICE_OPTIONS.h
  - modified/extra source: README_TAP_HACKS.txt, cost_test.F

## Input variants (input*/)
- **input**: data.pkg on: useShelfIce
  - data: deltaT=1800.0, nTimeSteps=20, nIter0=0, eosType='JMD95Z', implicitFreeSurface=.TRUE., nonHydrostatic=.FALSE., usingSphericalPolarGrid=.TRUE.
  - namelist files: data data.pkg data.shelfice eedata
- **input.htd**: data.pkg on: useShelfIce, useDiagnostics, useMNC
  - data: deltaT=1800., nTimeSteps=20, nIter0=8640, nonlinFreeSurf=4, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE.
  - namelist files: data data.diagnostics data.mnc data.pkg data.shelfice eedata
- **input.icefront**: data.pkg on: useDiagnostics, useShelfIce, useIcefront
  - namelist files: data.diagnostics data.icefront data.pkg data.shelfice eedata
- **input.obcs**: data.pkg on: useShelfIce, useGGL90, useOBCS, useDiagnostics
  - data: deltaT=1800.0, nTimeSteps=12, nIter0=0, nonlinFreeSurf=4, select_rStar=2, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.diagnostics data.ggl90 data.obcs data.pkg data.shelfice eedata
- **input.stic**: data.pkg on: useShelfIce, useSTIC, useDiagnostics
  - namelist files: data.diagnostics data.pkg data.shelfice data.stic
- **input_ad**: data.pkg on: useMNC, useShelfIce, useGrdchk
  - data: deltaT=1800.0, nTimeSteps=5, nIter0=8640, nonlinFreeSurf=2, eosType='JMD95Z', implicitFreeSurface=.TRUE., nonHydrostatic=.FALSE., usingSphericalPolarGrid=.TRUE.
  - namelist files: data data.autodiff data.cost data.ctrl data.grdchk data.mnc data.optim data.pkg data.shelfice eedata prepare_run
- **input_ad.htd**: data.pkg on: useShelfIce, useGrdchk
  - data: deltaT=1800.0, nTimeSteps=5, nIter0=8640, nonlinFreeSurf=2, eosType='JMD95Z', implicitFreeSurface=.TRUE., nonHydrostatic=.FALSE., usingSphericalPolarGrid=.TRUE., tempAdvScheme=7, saltAdvScheme=7
  - namelist files: data data.grdchk data.pkg data.shelfice
- **input_ad.stic**: data.pkg on: useMNC, useShelfIce, useSTIC, useGrdchk
  - namelist files: data.cal data.cost data.ctrl data.pkg data.shelfice data.stic prepare_run
- **input_tap**: data.pkg on: useShelfIce, useGrdchk
  - data: deltaT=1800.0, nTimeSteps=5, nIter0=8640, nonlinFreeSurf=0, eosType='JMD95Z', implicitFreeSurface=.TRUE., nonHydrostatic=.FALSE., usingSphericalPolarGrid=.TRUE.
  - namelist files: data data.pkg prepare_run

## Reference results
`output.htd.txt` `output.icefront.txt` `output.obcs.txt` `output.stic.txt` `output.txt` `output_adm.htd.txt` `output_adm.stic.txt` `output_adm.txt` `output_tap_adj.txt` `output_tap_tlm.txt`

Run: `cd verification; ./testreport -of <optfile> -t isomip` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
