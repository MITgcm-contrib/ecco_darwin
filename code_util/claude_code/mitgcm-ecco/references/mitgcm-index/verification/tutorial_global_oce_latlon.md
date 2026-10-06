# verification/tutorial_global_oce_latlon

## README (first 40 lines)
```
Tutorial Example: "Global ocean"
(Global Ocean Simulation at 4^o Resolution)
===========================================

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

There is comparison output in the directory:
  results/output.txt

Comments:

```

## Build variants (code*/)
- **code**: packages.conf = `gfd cd_code gmredi ptracers mnc`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor cd_code gmredi ptracers mnc
  - SIZE.h: grid 90x40x15; sNx=45, sNy=40, OLx=2, OLy=2, nSx=2, nSy=1, nPx=1, nPy=1, Nr=15; has SIZE.h_mpi
  - modified/extra source: ptracers_apply_forcing.F, ptracers_forcing_surf.F

## Input variants (input*/)
- **input**: data.pkg on: useGMRedi, usePTRACERS
  - data: deltaTmom=1800., deltaTtracer=86400., nTimeSteps=20, nIter0=0, eosType='JMD95Z', usingSphericalPolarGrid=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.gmredi data.pkg data.ptracers eedata

## Reference results
`output.txt`

Run: `cd verification; ./testreport -of <optfile> -t tutorial_global_oce_latlon` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
