# verification/tutorial_tracer_adjsens

## README (first 40 lines)
```
Tutorial Example: "Centennial Time Scale Tracer Injection"
==========================================================
(formerly "carbon" verification ;
 also "tracer_adjoint_sensitivity" in release.1 branch)

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
  ln -s ../build/mitgcmuv_ad .
  ./mitgcmuv_ad > output_adm.txt
  cd ..

There is comparison output in the directory:
  results/output_adm.txt
grep for grdchk output:
  grep 'grdchk output' output_adm.txt

Comments:
  The input data is real*4

```

## Build variants (code*/)
- **code_ad**: packages.conf = `gfd cd_code gmredi kpp ptracers autodiff cost ctrl grdchk`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor cd_code gmredi kpp ptracers autodiff cost ctrl grdchk
  - SIZE.h: grid 90x40x20; sNx=45, sNy=20, OLx=3, OLy=3, nSx=2, nSy=2, nPx=1, nPy=1, Nr=20; has SIZE.h_mpi
  - option/size headers: AUTODIFF_OPTIONS.h, COST_OPTIONS.h, CPP_OPTIONS.h, CTRL_OPTIONS.h, CTRL_SIZE.h, GAD_OPTIONS.h, GMREDI_OPTIONS.h
  - modified/extra source: MDSIO_BUFF_WH.h, ptracers_forcing_surf.F, tamc.h
- **code_tap**: packages.conf = `gfd monitor cd_code gmredi kpp ptracers tapenade adjoint`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor cd_code gmredi kpp ptracers tapenade autodiff cost ctrl grdchk
  - SIZE.h: grid 90x40x20; sNx=45, sNy=20, OLx=3, OLy=3, nSx=2, nSy=2, nPx=1, nPy=1, Nr=20; has SIZE.h_mpi
  - option/size headers: AUTODIFF_OPTIONS.h, COST_OPTIONS.h, CPP_OPTIONS.h, CTRL_OPTIONS.h, CTRL_SIZE.h, GAD_OPTIONS.h, GMREDI_OPTIONS.h
  - modified/extra source: ptracers_forcing_surf.F

## Input variants (input*/)
- **input_ad**: data.pkg on: usePtracers, useGMRedi, useGrdchk
  - data: deltaTmom=2400., endTime=345600., startTime=0., nonlinFreeSurf=3, select_rStar=1, eosType='LINEAR', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE.
  - namelist files: data data.autodiff data.cost data.ctrl data.gmredi data.grdchk data.optim data.pkg data.ptracers eedata prepare_run
- **input_ad.som81**: data.pkg on: -
  - data: deltaTmom=2400., endTime=345600., startTime=0., nonlinFreeSurf=3, select_rStar=1, eosType='LINEAR', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE., tempAdvScheme=80, saltAdvScheme=81
  - namelist files: data
- **input_tap**: data.pkg on: -
  - data: deltaTmom=2400., endTime=345600., startTime=0., nonlinFreeSurf=3, select_rStar=1, eosType='LINEAR', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE.
  - namelist files: data prepare_run

## Reference results
`output_adm.som81.txt` `output_adm.txt` `output_tap_adj.txt` `output_tap_tlm.txt` `output_tlm.som81.txt.gz` `output_tlm.txt.gz`

Run: `cd verification; ./testreport -of <optfile> -t tutorial_tracer_adjsens` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
