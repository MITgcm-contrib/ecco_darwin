# verification/fizhi-cs-32x32x40

## Build variants (code*/)
- **code**: packages.conf = `exch2 gfd -mom_fluxform shap_filt fizhi gridalt diagnostics`
  - expanded: exch2 mom_common mom_vecinv generic_advdiff debug mdsio rw monitor shap_filt fizhi gridalt diagnostics
  - SIZE.h: grid 192x32x40; sNx=32, sNy=32, OLx=2, OLy=2, nSx=6, nSy=1, nPx=1, nPy=1, Nr=40; has SIZE.h_mpi
  - option/size headers: CPP_OPTIONS.h, fizhi_SIZE.h

## Input variants (input*/)
- **input**: data.pkg on: useFizhi, useGridAlt, useSHAP_FILT, useDiagnostics
  - data: deltaT=120.0, nTimeSteps=6, nIter0=0, nonlinFreeSurf=4, select_rStar=2, eosType='IDEALG', implicitFreeSurface=.TRUE., usingCurvilinearGrid=.TRUE., buoyancyRelation='ATMOSPHERIC', saltAdvScheme=2
  - namelist files: data data.diagnostics data.fizhi data.gcmo3 data.pkg data.sage data.shap eedata prepare_run

## Reference results
`output.txt`

Run: `cd verification; ./testreport -of <optfile> -t fizhi-cs-32x32x40` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
