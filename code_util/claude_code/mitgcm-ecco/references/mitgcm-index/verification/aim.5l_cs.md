# verification/aim.5l_cs

## Build variants (code*/)
- **code**: packages.conf = `exch2 atmospheric -mom_fluxform aim_v23 land thsice diagnostics mnc`
  - expanded: exch2 mom_common mom_vecinv generic_advdiff debug mdsio rw monitor shap_filt aim_v23 land thsice diagnostics mnc
  - SIZE.h: grid 192x32x5; sNx=32, sNy=32, OLx=2, OLy=2, nSx=6, nSy=1, nPx=1, nPy=1, Nr=5; has SIZE.h_mpi
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h, DIAG_OPTIONS.h, SHAP_FILT_OPTIONS.h
  - modified/extra source: mom_vi_hfacz_diss.F, mom_vi_mask_vort3.F

## Input variants (input*/)
- **input**: data.pkg on: useAIM, useLand, useSHAP_FILT
  - data: deltaT=450.0, nTimeSteps=10, nIter0=69120, nonlinFreeSurf=4, select_rStar=2, eosType='IDEALG', implicitFreeSurface=.TRUE., usingCurvilinearGrid=.TRUE., buoyancyRelation='ATMOSPHERIC', saltAdvScheme=3
  - namelist files: data data.aimphys data.diagnostics data.land data.pkg data.shap eedata
- **input.thSI**: data.pkg on: useAIM, useLand, useThSIce, useSHAP_FILT, useDiagnostics
  - data: deltaT=450.0, nTimeSteps=10, nIter0=0, nonlinFreeSurf=4, select_rStar=2, eosType='IDEALG', implicitFreeSurface=.TRUE., usingCurvilinearGrid=.TRUE., buoyancyRelation='ATMOSPHERIC', saltAdvScheme=3
  - namelist files: data data.aimphys data.diagnostics data.ice data.land data.pkg data.shap

## Reference results
`output.thSI.txt` `output.txt`

Run: `cd verification; ./testreport -of <optfile> -t aim.5l_cs` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
