# verification/cheapAML_box

## Build variants (code*/)
- **code**: packages.conf = `gfd cheapaml diagnostics mnc`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor cheapaml diagnostics mnc
  - SIZE.h: grid 100x100x1; sNx=50, sNy=25, OLx=3, OLy=3, nSx=2, nSy=4, nPx=1, nPy=1, Nr=1; has SIZE.h_mpi
  - option/size headers: CHEAPAML_OPTIONS.h, CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h

## Input variants (input*/)
- **input**: data.pkg on: useCheapAML, useDiagnostics, useMNC
  - data: deltaT=1200., nIter0=0, endTime=28800., eosType='LINEAR', usingSphericalPolarGrid=.TRUE., tempAdvScheme=33
  - namelist files: data data.cheapaml data.diagnostics data.mnc data.pkg eedata
- **input.lanl**: data.pkg on: useCheapAML, useDiagnostics
  - data: deltaT=1200., nIter0=0, endTime=28800., eosType='LINEAR', usingSphericalPolarGrid=.TRUE., tempAdvScheme=33
  - namelist files: data data.cheapaml data.pkg

## Reference results
`output.lanl.txt` `output.txt`

Run: `cd verification; ./testreport -of <optfile> -t cheapAML_box` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
