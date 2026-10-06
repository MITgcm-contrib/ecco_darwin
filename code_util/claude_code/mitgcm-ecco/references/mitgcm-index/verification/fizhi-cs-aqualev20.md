# verification/fizhi-cs-aqualev20

## Build variants (code*/)
- **code**: packages.conf = `exch2 gfd -mom_fluxform shap_filt fizhi gridalt diagnostics mnc`
  - expanded: exch2 mom_common mom_vecinv generic_advdiff debug mdsio rw monitor shap_filt fizhi gridalt diagnostics mnc
  - SIZE.h: grid 192x32x20; sNx=32, sNy=32, OLx=2, OLy=2, nSx=6, nSy=1, nPx=1, nPy=1, Nr=20; has SIZE.h_mpi
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h, FIZHI_OPTIONS.h, fizhi_SIZE.h
  - modified/extra source: fizhi_init_fixed.F, fizhi_init_veg.F, fizhi_init_vegsurftiles.F, ini_salt.F, update_ocean_exports.F

## Input variants (input*/)
- **input**: data.pkg on: useFizhi, useGridAlt, useSHAP_FILT, useDiagnostics, useMNC
  - data: deltaT=120.0, nTimeSteps=7, nIter0=0, nonlinFreeSurf=4, select_rStar=2, eosType='IDEALG', implicitFreeSurface=.TRUE., usingCurvilinearGrid=.TRUE., buoyancyRelation='ATMOSPHERIC', saltAdvScheme=3
  - namelist files: data data.diagnostics data.fizhi data.mnc data.pkg data.sage data.shap eedata prepare_run

## Reference results
`output.txt`

Run: `cd verification; ./testreport -of <optfile> -t fizhi-cs-aqualev20` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
