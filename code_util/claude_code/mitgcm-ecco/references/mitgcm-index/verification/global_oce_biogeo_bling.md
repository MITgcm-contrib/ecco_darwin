# verification/global_oce_biogeo_bling

## Build variants (code*/)
- **code**: packages.conf = `obsfit cal gfd cd_code gmredi ptracers gchem bling diagnostics mnc`
  - expanded: obsfit cal mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor cd_code gmredi ptracers gchem bling diagnostics mnc
  - SIZE.h: grid 128x64x15; sNx=64, sNy=32, OLx=4, OLy=4, nSx=2, nSy=2, nPx=1, nPy=1, Nr=15; has SIZE.h_mpi
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h, PTRACERS_SIZE.h
- **code_ad**: packages.conf = `exch2 gfd cd_code gmredi cal ptracers gchem bling profiles obsfit diagnostics ecco autodiff cost ctrl grdchk`
  - expanded: exch2 mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor cd_code gmredi cal ptracers gchem bling profiles obsfit diagnostics ecco autodiff cost ctrl grdchk
  - SIZE.h: grid 128x64x15; sNx=64, sNy=32, OLx=4, OLy=4, nSx=2, nSy=2, nPx=1, nPy=1, Nr=15; has SIZE.h_mpi
  - option/size headers: BLING_OPTIONS.h, CPP_OPTIONS.h, CTRL_OPTIONS.h, CTRL_SIZE.h, DIAGNOSTICS_SIZE.h, ECCO_OPTIONS.h, GMREDI_OPTIONS.h, PROFILES_SIZE.h, PTRACERS_SIZE.h
  - modified/extra source: MDSIO_BUFF_WH.h, tamc.h
- **code_tap**: packages.conf = `exch2 gfd cd_code gmredi cal ptracers gchem bling profiles diagnostics ecco autodiff cost ctrl grdchk tapenade`
  - expanded: exch2 mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor cd_code gmredi cal ptracers gchem bling profiles diagnostics ecco autodiff cost ctrl grdchk tapenade
  - SIZE.h: grid 128x64x15; sNx=64, sNy=32, OLx=4, OLy=4, nSx=2, nSy=2, nPx=1, nPy=1, Nr=15; has SIZE.h_mpi
  - option/size headers: AUTODIFF_OPTIONS.h, BLING_OPTIONS.h, CPP_OPTIONS.h, CTRL_OPTIONS.h, CTRL_SIZE.h, DIAGNOSTICS_SIZE.h, ECCO_OPTIONS.h, GMREDI_OPTIONS.h, PROFILES_SIZE.h, PTRACERS_SIZE.h

## Input variants (input*/)
- **input**: data.pkg on: useGMRedi, usePTRACERS, useGCHEM, useDiagnostics
  - data: deltaTmom=900., deltaTtracer=43200., nTimeSteps=4, nIter0=0, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE., tempAdvScheme=2, saltAdvScheme=2
  - namelist files: data data.bling data.diagnostics data.gchem data.gmredi data.pkg data.ptracers eedata prepare_run
- **input_ad**: data.pkg on: useGMRedi, usePTRACERS, useGCHEM, useECCO, useCAL, useGrdchk, usePROFILES, useDIAGNOSTICS
  - data: deltaT=900., nTimeSteps=4, nIter0=0, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE., tempAdvScheme=2, saltAdvScheme=2
  - namelist files: data data.autodiff data.bling data.cal data.cost data.ctrl data.diagnostics data.ecco data.err data.gchem data.gmredi data.grdchk data.optim data.pkg data.profiles data.ptracers eedata prepare_run
- **input_ad.obsfit**: data.pkg on: useGMRedi, usePTRACERS, useGCHEM, useECCO, useGrdchk, useOBSFIT, useDIAGNOSTICS
  - namelist files: data.ecco data.grdchk data.obsfit data.pkg
- **input_tap**: data.pkg on: -
  - namelist files: prepare_run

## Reference results
`output.txt` `output_adm.obsfit.txt` `output_adm.txt` `output_tap_adj.txt` `output_tlm.obsfit.txt.gz` `output_tlm.txt.gz`

Run: `cd verification; ./testreport -of <optfile> -t global_oce_biogeo_bling` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
