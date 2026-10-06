# verification/cfc_example

## Build variants (code*/)
- **code**: packages.conf = `gfd cd_code gmredi ptracers gchem cfc layers diagnostics`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor cd_code gmredi ptracers gchem cfc layers diagnostics
  - SIZE.h: grid 128x64x15; sNx=64, sNy=32, OLx=4, OLy=4, nSx=2, nSy=2, nPx=1, nPy=1, Nr=15; has SIZE.h_mpi
  - option/size headers: DIAGNOSTICS_SIZE.h, GCHEM_OPTIONS.h, GMREDI_OPTIONS.h, LAYERS_OPTIONS.h, LAYERS_SIZE.h, PTRACERS_SIZE.h
  - modified/extra source: MDSIO_BUFF_3D.h

## Input variants (input*/)
- **input**: data.pkg on: useGMRedi, usePTRACERS, useGCHEM, useDiagnostics, useLayers
  - data: deltaTmom=900., deltaTtracer=43200., nTimeSteps=4, nIter0=4269600, eosType='POLY3', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE.
  - namelist files: data data.cfc data.diagnostics data.gchem data.gmredi data.layers data.pkg data.ptracers eedata prepare_run

## Reference results
`output.txt`

Run: `cd verification; ./testreport -of <optfile> -t cfc_example` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
