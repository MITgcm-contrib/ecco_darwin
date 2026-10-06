# verification/obcs_ctrl

## Build variants (code*/)
- **code_ad**: packages.conf = `gfd -mom_fluxform obcs exf cal diagnostics ecco autodiff cost ctrl grdchk`
  - expanded: mom_common mom_vecinv generic_advdiff debug mdsio rw monitor obcs exf cal diagnostics ecco autodiff cost ctrl grdchk
  - SIZE.h: grid 64x64x8; sNx=32, sNy=32, OLx=2, OLy=2, nSx=2, nSy=2, nPx=1, nPy=1, Nr=8; has SIZE.h_mpi
  - option/size headers: AUTODIFF_OPTIONS.h, CPP_OPTIONS.h, CTRL_OPTIONS.h, EXF_OPTIONS.h, OBCS_OPTIONS.h
  - modified/extra source: tamc.h

## Input variants (input*/)
- **input_ad**: data.pkg on: useECCO, useOBCS, useEXF, useDiagnostics, useGrdchk
  - data: deltaTmom=1200.0, deltaTtracer=1200.0, endTime=4800., startTime=0., eosType='LINEAR', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE., tempAdvScheme=30, saltAdvScheme=30
  - namelist files: data data.autodiff data.cal data.cost data.ctrl data.diagnostics data.ecco data.err data.exf data.grdchk data.obcs data.optim data.pkg eedata

## Reference results
`output_adm.txt`

Run: `cd verification; ./testreport -of <optfile> -t obcs_ctrl` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
