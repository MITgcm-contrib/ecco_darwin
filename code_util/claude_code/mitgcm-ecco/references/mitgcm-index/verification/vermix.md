# verification/vermix

## Build variants (code*/)
- **code**: packages.conf = `gfd kpp pp81 my82 ggl90 opps diagnostics mnc`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor kpp pp81 my82 ggl90 opps diagnostics mnc
  - SIZE.h: grid 1x1x26; sNx=1, sNy=1, OLx=2, OLy=2, nSx=1, nSy=1, nPx=1, nPy=1, Nr=26
  - option/size headers: DIAGNOSTICS_SIZE.h, GGL90_OPTIONS.h, KPP_OPTIONS.h

## Input variants (input*/)
- **input**: data.pkg on: useKPP, useDiagnostics, useMNC
  - data: deltaT=1200., nTimeSteps=20, nIter0=0, eosType='MDJWF', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., tempAdvScheme=33
  - namelist files: data data.diagnostics data.kpp data.mnc data.pkg eedata
- **input.dd**: data.pkg on: -
  - data: deltaT=1200., nTimeSteps=20, nIter0=0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., tempAdvScheme=33
  - namelist files: data data.kpp
- **input.ggl90**: data.pkg on: useGGL90, useDiagnostics, useMNC
  - namelist files: data.diagnostics data.ggl90 data.pkg
- **input.gglLC**: data.pkg on: useGGL90, useDiagnostics
  - namelist files: data.diagnostics data.ggl90 data.pkg
- **input.my82**: data.pkg on: useMY82, useDiagnostics, useMNC
  - namelist files: data.diagnostics data.my82 data.pkg
- **input.opps**: data.pkg on: useOPPS, useDiagnostics, useMNC
  - namelist files: data.diagnostics data.opps data.pkg
- **input.pp81**: data.pkg on: usePP81, useDiagnostics, useMNC
  - namelist files: data.diagnostics data.pkg data.pp81

## Reference results
`output.dd.txt` `output.ggl90.txt` `output.gglLC.txt` `output.my82.txt` `output.opps.txt` `output.pp81.txt` `output.txt`

Run: `cd verification; ./testreport -of <optfile> -t vermix` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
