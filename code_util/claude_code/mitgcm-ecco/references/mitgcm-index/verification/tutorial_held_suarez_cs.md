# verification/tutorial_held_suarez_cs

## Build variants (code*/)
- **code**: packages.conf = `exch2 gfd shap_filt diagnostics mnc`
  - expanded: exch2 mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor shap_filt diagnostics mnc
  - SIZE.h: grid 192x32x20; sNx=32, sNy=32, OLx=2, OLy=2, nSx=6, nSy=1, nPx=1, nPy=1, Nr=20; has SIZE.h_mpi
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h
  - modified/extra source: apply_forcing.F

## Input variants (input*/)
- **input**: data.pkg on: useSHAP_FILT, useDiagnostics
  - data: deltaT=450., nTimeSteps=16, startTime=124416000., nonlinFreeSurf=4, select_rStar=2, eosType='IDEALG', implicitFreeSurface=.TRUE., usingCurvilinearGrid=.TRUE., buoyancyRelation='ATMOSPHERIC'
  - namelist files: data data.diagnostics data.pkg data.shap eedata

## Reference results
`output.txt`

Run: `cd verification; ./testreport -of <optfile> -t tutorial_held_suarez_cs` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
