# verification/halfpipe_streamice

## Build variants (code*/)
- **code**: packages.conf = `gfd -mom_common -mom_fluxform -mom_vecinv -generic_advdiff streamice diagnostics`
  - expanded: debug mdsio rw monitor streamice diagnostics
  - SIZE.h: grid 40x20x1; sNx=20, sNy=20, OLx=3, OLy=3, nSx=2, nSy=1, nPx=1, nPy=1, Nr=1; has SIZE.h_mpi
  - option/size headers: DIAGNOSTICS_SIZE.h, STREAMICE_OPTIONS.h
- **code_ad**: packages.conf = `gfd -mom_common -mom_fluxform -mom_vecinv -generic_advdiff streamice diagnostics adjoint`
  - expanded: debug mdsio rw monitor streamice diagnostics autodiff cost ctrl grdchk
  - SIZE.h: grid 40x20x1; sNx=20, sNy=20, OLx=3, OLy=3, nSx=2, nSy=1, nPx=1, nPy=1, Nr=1; has SIZE.h_mpi
  - option/size headers: AUTODIFF_OPTIONS.h, COST_OPTIONS.h, CTRL_OPTIONS.h, CTRL_SIZE.h, DIAGNOSTICS_SIZE.h, STREAMICE_OPTIONS.h
  - modified/extra source: cost_test.F, tamc.h
- **code_tap**: packages.conf = `gfd -mom_common -mom_fluxform -mom_vecinv -generic_advdiff streamice diagnostics tapenade adjoint`
  - expanded: debug mdsio rw monitor streamice diagnostics tapenade autodiff cost ctrl grdchk
  - SIZE.h: grid 40x20x1; sNx=20, sNy=20, OLx=3, OLy=3, nSx=2, nSy=1, nPx=1, nPy=1, Nr=1; has SIZE.h_mpi
  - option/size headers: AUTODIFF_OPTIONS.h, COST_OPTIONS.h, CTRL_OPTIONS.h, CTRL_SIZE.h, DIAGNOSTICS_SIZE.h, STREAMICE_OPTIONS.h

## Input variants (input*/)
- **input**: data.pkg on: useStreamIce, useDiagnostics
  - data: deltaT=6307200., nTimeSteps=10, startTime=0., usingCartesianGrid=.TRUE., momStepping=.FALSE.
  - namelist files: data data.diagnostics data.pkg data.streamice data.streamice_geomSetup eedata
- **input_ad**: data.pkg on: useStreamIce, useDiagnostics, useGrdchk
  - data: deltaT=6307200., nTimeSteps=3, startTime=0., usingCartesianGrid=.TRUE., momStepping=.FALSE.
  - namelist files: data data.autodiff data.cost data.ctrl data.diagnostics data.grdchk data.optim data.pkg data.streamice eedata
- **input_tap**: data.pkg on: useStreamIce, useGrdchk
  - namelist files: data.pkg data.streamice prepare_run

## Reference results
`output.txt` `output_adm.txt` `output_tap_adj.txt.gz` `output_tap_tlm.txt.gz`

Run: `cd verification; ./testreport -of <optfile> -t halfpipe_streamice` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
