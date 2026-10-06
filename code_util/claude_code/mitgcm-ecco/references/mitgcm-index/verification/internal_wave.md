# verification/internal_wave

## README (first 40 lines)
```
Example: "Internal Wave Forced by Open Boundary"
================================================

This uses a simply EOS (rho = - rho_o alpha T')
 and treats salt as a passive tracer.

Configure and compile the code:
  cd build
  ../../../tools/genmake2 -mods ../code [-of my_platform_optionFile]
  make depend
  make
  cd ..

To run:
  cd run
  ln -s ../input/* .
  ln -s ../build/mitgcmuv .
  ./mitgcmuv > output.txt
  cd ..

There is comparison output in:
  results/output.txt

Comments:
  The input data is real*8 and generated using the MATLAB script
  input/gendata.m

```

## Build variants (code*/)
- **code**: packages.conf = `gfd obcs kl10 mnc`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor obcs kl10 mnc
  - SIZE.h: grid 60x1x20; sNx=30, sNy=1, OLx=2, OLy=2, nSx=2, nSy=1, nPx=1, nPy=1, Nr=20; has SIZE.h_mpi
  - option/size headers: CPP_OPTIONS.h
  - modified/extra source: obcs_calc.F

## Input variants (input*/)
- **input**: data.pkg on: useOBCS, useMNC
  - data: deltaT=500., nTimeSteps=100, nIter0=0, nonlinFreeSurf=3, eosType='LINEAR', implicitFreeSurface=.TRUE., nonHydrostatic=.FALSE., usingCartesianGrid=.TRUE.
  - namelist files: data data.mnc data.obcs data.pkg eedata
- **input.kl10**: data.pkg on: useOBCS, useKL10
  - data: deltaT=450., nTimeSteps=300, nIter0=0, nonlinFreeSurf=3, eosType='LINEAR', implicitFreeSurface=.TRUE., nonHydrostatic=.FALSE., usingCartesianGrid=.TRUE.
  - namelist files: data data.kl10 data.pkg

## Reference results
`output.kl10.txt` `output.txt`

Run: `cd verification; ./testreport -of <optfile> -t internal_wave` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
