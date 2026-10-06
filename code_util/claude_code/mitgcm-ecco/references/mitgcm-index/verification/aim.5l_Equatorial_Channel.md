# verification/aim.5l_Equatorial_Channel

## README (first 40 lines)
```
Five level intermediate atmospheric physics example.
 Experimental configuration is an equatorial channel
 of 65 degrees wide cented on the Equator. A warm SST anomaly
 with a gaussian shape is defined at the center of the domain
 with an extension of ~30 degrees meridionally and ~60 degrees
 zonally.

The local copy of S/R aim_surf_bc.F (in code directory) provide
a) the SST field (hard coded).
b) the time in the year is fixed and corresponds to the spring equinox;
 this ensures a symetric forcing between N & S hemispheres.
The model uses the standard 5-level vertical resolution of the Speedy_v23
 code; the horizontal resolution is 2.8125 with 128x23 grid points.

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
Notes:
* To produce a short test that is relevant for the Atmosphere physics part,
  a restart file (pickup.0000051840) is included. The model reaches
  a statistical equilibrium after 1 year.
* the file aim_SST.* contains the SST field.
* Since aim pkg uses arrays with MAX_NO_THREADS as dimension, the maximum
... (truncated)
```

## Build variants (code*/)
- **code**: packages.conf = `atmospheric aim_v23`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor shap_filt aim_v23
  - SIZE.h: grid 128x23x5; sNx=32, sNy=23, OLx=2, OLy=2, nSx=4, nSy=1, nPx=1, nPy=1, Nr=5; has SIZE.h_mpi
  - modified/extra source: aim_surf_bc.F, ini_depths.F

## Input variants (input*/)
- **input**: data.pkg on: useAIM, useSHAP_FILT
  - data: deltaT=600., nTimeSteps=10, nIter0=51840, eosType='IDEALG', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE., buoyancyRelation='ATMOSPHERIC', saltAdvScheme=3
  - namelist files: data data.aimphys data.pkg data.shap eedata

## Reference results
`output.txt`

Run: `cd verification; ./testreport -of <optfile> -t aim.5l_Equatorial_Channel` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
