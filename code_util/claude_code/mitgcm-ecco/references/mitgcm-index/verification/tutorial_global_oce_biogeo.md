# verification/tutorial_global_oce_biogeo

## README (first 40 lines)
```
Tutorial Example: "Biochemistry Tutorial"
=========================================
(formerly "dic_example" verification)

--- Instructions for Forward tests:

Configure and compile the code:
  cd build
  ../../../tools/genmake2 -mods ../code [-of my_platform_optionFile]
  make depend
  make
  cd ..

To run:
  cd run
  ln -s ../input/* .
  ./prepare_run
  ../build/mitgcmuv > output.txt
  cd ..

There is comparison output in the directory:
  results/output.txt

-- Instructions for Adjoint tests:
Note: This requires access to a TAF license.

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
  ../build/mitgcmuv_ad  > output_adm.txt
  cd ..

... (truncated)
```

## Build variants (code*/)
- **code**: packages.conf = `gfd cd_code gmredi ptracers gchem dic diagnostics mnc`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor cd_code gmredi ptracers gchem dic diagnostics mnc
  - SIZE.h: grid 128x64x15; sNx=64, sNy=32, OLx=4, OLy=4, nSx=2, nSy=2, nPx=1, nPy=1, Nr=15; has SIZE.h_mpi
  - option/size headers: DIAGNOSTICS_SIZE.h, DIAG_OPTIONS.h, PTRACERS_SIZE.h
- **code_ad**: packages.conf = `gfd cd_code gmredi ptracers gchem dic autodiff cost ctrl grdchk`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor cd_code gmredi ptracers gchem dic autodiff cost ctrl grdchk
  - SIZE.h: grid 128x64x15; sNx=64, sNy=32, OLx=4, OLy=4, nSx=2, nSy=2, nPx=1, nPy=1, Nr=15; has SIZE.h_mpi
  - option/size headers: COST_OPTIONS.h, CPP_OPTIONS.h, CTRL_OPTIONS.h, CTRL_SIZE.h, GMREDI_OPTIONS.h, PTRACERS_SIZE.h
  - modified/extra source: cost_tracer.F, tamc.h
- **code_tap**: packages.conf = `gfd monitor cd_code gmredi ptracers gchem dic tapenade adjoint`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor cd_code gmredi ptracers gchem dic tapenade autodiff cost ctrl grdchk
  - SIZE.h: grid 128x64x15; sNx=64, sNy=32, OLx=4, OLy=4, nSx=2, nSy=2, nPx=1, nPy=1, Nr=15; has SIZE.h_mpi
  - option/size headers: AUTODIFF_OPTIONS.h, COST_OPTIONS.h, CPP_OPTIONS.h, CTRL_OPTIONS.h, CTRL_SIZE.h, GMREDI_OPTIONS.h, PTRACERS_SIZE.h
  - modified/extra source: cost_tracer.F

## Input variants (input*/)
- **input**: data.pkg on: useGMRedi, usePTRACERS, useGCHEM, useDiagnostics
  - data: deltaTmom=900., deltaTtracer=43200., nTimeSteps=4, nIter0=5184000, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE., tempAdvScheme=2, saltAdvScheme=2
  - namelist files: data data.diagnostics data.dic data.gchem data.gmredi data.pkg data.ptracers eedata
- **input_ad**: data.pkg on: useGMRedi, usePTRACERS, useGCHEM, useGrdchk
  - data: deltaTmom=900., deltaTtracer=43200., nTimeSteps=4, nIter0=5184000, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE., tempAdvScheme=2, saltAdvScheme=2
  - namelist files: data data.autodiff data.cost data.ctrl data.dic data.gchem data.gmredi data.grdchk data.optim data.pkg data.ptracers eedata prepare_run
- **input_tap**: data.pkg on: -
  - namelist files: prepare_run

## Reference results
`output.txt` `output_adm.txt` `output_tap_adj.txt` `output_tap_tlm.txt` `output_tlm.txt.gz`

Run: `cd verification; ./testreport -of <optfile> -t tutorial_global_oce_biogeo` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
