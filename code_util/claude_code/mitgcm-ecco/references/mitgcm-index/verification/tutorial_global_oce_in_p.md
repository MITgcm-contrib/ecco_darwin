# verification/tutorial_global_oce_in_p

## README (first 40 lines)
```
Tutorial Example: "P coordinate Global Ocean"
(Global Ocean Simulation at 4o Resolution in Pressure Coordinates)
==================================================================
(formerly "global_ocean_pressure")

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
o the set up is similar to that of tutorial_global_oce_latlon
o the code directory contains calc_phi_hyd.F, where the potential is computed
  according to the more natural finite volume discretization. Finite difference
  discretization is energy conserving, but the representation of the "fixed"
  surface (interface ocean-atmosphere) is less consistent.
o the code directory also contains dynamics.F which calls
  remove_mean_rl.F, a generic routine, to remove the mean from the
  diagnostic variable phiHydLow (sea surface height/gravity in pressure
  coordinates)

changes: 07 Feb. 2003 (jmc):
o find difficult to maintain the local version of dynamics.F up to date.
  therefore, has been remove from the code directory.
  One can recover the same version (but up to date) simply
  by activating the commented lines [between lines Cml( and Cml) ],
  at the end of the standard version of dynamics.F
o finite volume form of calc_phi_hyd.F is now a standard option.
... (truncated)
```

## Build variants (code*/)
- **code**: packages.conf = `gfd -mom_vecinv`
  - expanded: mom_common mom_fluxform generic_advdiff debug mdsio rw monitor
  - SIZE.h: grid 90x40x15; sNx=90, sNy=20, OLx=2, OLy=2, nSx=1, nSy=2, nPx=1, nPy=1, Nr=15; has SIZE.h_mpi
  - option/size headers: CPP_OPTIONS.h

## Input variants (input*/)
- **input**: data.pkg on: -
  - data: deltaTmom=1200.0, deltaTtracer=172800.0, endTime=3456000., startTime=0., nonlinFreeSurf=4, eosType='JMD95P', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE., buoyancyRelation='OCEANICP'
  - namelist files: data data.pkg eedata

## Reference results
`output.txt`

Run: `cd verification; ./testreport -of <optfile> -t tutorial_global_oce_in_p` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
