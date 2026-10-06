# verification/cpl_aim+ocn

## README (first 40 lines)
```
Atmosphere-Ocean coupled set-up example "cpl_aim+ocn"
================================================================================
using simplified atmospheric physics (AIM), in realistic configuration (orography
& continent) with land and seaice component, on cubed-sphere (cs-32) grid.

### Overview:
Uses "in-house" MITgcm coupler<br>
(`pkg/atm_ocn_coupler`, `pkg/compon_communic`, `pkg/atm_compon_interf`,
`pkg/ocn_compon_interf`)<br>
with each component config and customized src code in: `code_cpl`, `code_atm`,
`code_ocn` ;<br>
and input parameter files in: `input_cpl`, `input_atm`, `input_ocn`.

- Atmos set-up and parameter is similar to `aim_5l_cs/` experiment
- Ocean set-up and parameter is similar to `global_ocean.cs32x15/` experiment

Requires the use of MPI; as default, use 1 proc for each component.

### Instructions:
To help getting started with this coupled set-up, the bash script
[../../tools/run_cpl_test](https://github.com/MITgcm/MITgcm/blob/master/tools/run_cpl_test)
(a short option summary is displayed when run without argument)
is provided and detailed instructions follow.

To clean everything:

    ../../tools/run_cpl_test 0

Configure and compile, e.g., using gfortran optfile:

    ../../tools/run_cpl_test 1 -of ../../tools/build_options/linux_amd64_gfortran

To run primary setup, thermodynamic seaice only (no seaice dynamics):

    ../../tools/run_cpl_test 2
    ../../tools/run_cpl_test 3

Step 2 above copies input files and directories, step 3 runs the coupled model.

To run secondary test (with seaice dynamics as part of ocean component), using
... (truncated)
```

## Build variants (code*/)
- **code_atm**: packages.conf = `exch2 atm_compon_interf compon_communic gfd -mom_fluxform shap_filt aim_v23 land thsice diagnostics`
  - expanded: exch2 atm_compon_interf compon_communic mom_common mom_vecinv generic_advdiff debug mdsio rw monitor shap_filt aim_v23 land thsice diagnostics
  - SIZE.h: grid 192x32x5; sNx=32, sNy=32, OLx=2, OLy=2, nSx=6, nSy=1, nPx=1, nPy=1, Nr=5
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h
- **code_cpl**: packages.conf = `atm_ocn_coupler compon_communic`
  - expanded: atm_ocn_coupler compon_communic
  - modified/extra source: ATMSIZE.h, OCNSIZE.h
- **code_ocn**: packages.conf = `exch2 ocn_compon_interf compon_communic gfd gmredi thsice seaice salt_plume diagnostics`
  - expanded: exch2 ocn_compon_interf compon_communic mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor gmredi thsice seaice salt_plume diagnostics
  - SIZE.h: grid 192x32x15; sNx=32, sNy=32, OLx=4, OLy=4, nSx=6, nSy=1, nPx=1, nPy=1, Nr=15
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h, SEAICE_OPTIONS.h

## Input variants (input*/)
- **input_atm**: data.pkg on: useAIM, useLand, useThSIce, useSHAP_FILT
  - data: deltaT=450.0, nTimeSteps=40, nIter0=0, nonlinFreeSurf=4, select_rStar=2, eosType='IDEALG', implicitFreeSurface=.TRUE., usingCurvilinearGrid=.TRUE., buoyancyRelation='ATMOSPHERIC', saltAdvScheme=3
  - namelist files: data data.aimphys data.cpl data.ice data.land data.pkg data.shap eedata prepare_run
- **input_atm.icedyn**: data.pkg on: useAIM, useLand, useThSIce, useSHAP_FILT, useDiagnostics
  - namelist files: data.diagnostics data.ice data.pkg
- **input_cpl**: data.pkg on: -
  - namelist files: data.cpl
- **input_cpl.icedyn**: data.pkg on: -
  - namelist files: data.cpl
- **input_ocn**: data.pkg on: useGMRedi, useDiagnostics
  - data: deltaTmom=3600., deltaTtracer=3600., nTimeSteps=5, nIter0=0, nonlinFreeSurf=4, select_rStar=2, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingCurvilinearGrid=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.cpl data.diagnostics data.gmredi data.pkg eedata prepare_run
- **input_ocn.icedyn**: data.pkg on: useGMRedi, useThSIce, useSEAICE, useDiagnostics
  - namelist files: data.diagnostics data.ice data.pkg data.salt_plume data.seaice

## Reference results
`atmSTDOUT.0000` `atmSTDOUT.icedyn` `ocnSTDOUT.0000` `ocnSTDOUT.icedyn`

Run: `cd verification; ./testreport -of <optfile> -t cpl_aim+ocn` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
