# verification/tutorial_dic_adjoffline

## Build variants (code*/)
- **code_ad**: packages.conf = `gfd -mom_common -mom_fluxform -mom_vecinv gmredi offline ptracers gchem dic mnc autodiff cost ctrl grdchk`
  - expanded: generic_advdiff debug mdsio rw monitor gmredi offline ptracers gchem dic mnc autodiff cost ctrl grdchk
  - SIZE.h: grid 128x64x15; sNx=32, sNy=32, OLx=4, OLy=4, nSx=4, nSy=2, nPx=1, nPy=1, Nr=15; has SIZE.h_mpi
  - option/size headers: AUTODIFF_OPTIONS.h, COST_OPTIONS.h, CPP_OPTIONS.h, CTRL_OPTIONS.h, DIC_OPTIONS.h, GMREDI_OPTIONS.h, PTRACERS_SIZE.h
  - modified/extra source: MDSIO_BUFF_WH.h, tamc.h

## Input variants (input*/)
- **input_ad**: data.pkg on: useGMRedi, usePTRACERS, useGCHEM, useMNC, useOFFLINE, useGrdchk
  - data: deltaTmom=900., deltaTtracer=43200., nTimeSteps=5, nIter0=0, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE., tempAdvScheme=2, saltAdvScheme=2
  - namelist files: data data.autodiff data.cost data.ctrl data.dic data.gchem data.gmredi data.grdchk data.mnc data.off data.optim data.pkg data.ptracers eedata prepare_run

## Reference results
`output_adm.txt` `output_tlm.txt.gz`

Run: `cd verification; ./testreport -of <optfile> -t tutorial_dic_adjoffline` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
