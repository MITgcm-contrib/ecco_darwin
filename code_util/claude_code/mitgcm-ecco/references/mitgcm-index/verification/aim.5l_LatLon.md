# verification/aim.5l_LatLon

## README (first 40 lines)
```
Intermediate Atmospheric physics, 5 layers Molteni Physics package.
Global spherical-grid configuration, 128x64x5 resolution.
====================================================================

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

---------------------------
Note:
Originally, this set up was very close to the one in
/development/adcroft/atmos/verification/molteni.128x64x5
with few modifications taken from run on hyades.

Others modifications have been added to improve the stability of the
model (fixed some bugs) and to get a less diffuse Q distribution:
o 3rd order scheme for the Horizontal advection of Q
o changes in the mapping between C-grid and A-grid  for surface stress
---------------------------
```

## Build variants (code*/)
- **code**: packages.conf = `atmospheric aim_v23 zonal_filt`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor shap_filt aim_v23 zonal_filt
  - SIZE.h: grid 128x64x5; sNx=128, sNy=16, OLx=3, OLy=3, nSx=1, nSy=4, nPx=1, nPy=1, Nr=5; has SIZE.h_mpi

## Input variants (input*/)
- **input**: data.pkg on: useAIM, useSHAP_FILT, useZONAL_FILT
  - data: deltaT=450.0, nTimeSteps=10, nIter0=69120, eosType='IDEALG', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE., buoyancyRelation='ATMOSPHERIC', saltAdvScheme=3
  - namelist files: data data.aimphys data.pkg data.shap data.zonfilt eedata

## Reference results
`output.txt`

Run: `cd verification; ./testreport -of <optfile> -t aim.5l_LatLon` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
