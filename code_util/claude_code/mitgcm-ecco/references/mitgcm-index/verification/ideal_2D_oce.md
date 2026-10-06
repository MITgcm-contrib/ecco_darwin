# verification/ideal_2D_oce

## README (first 40 lines)
```
test - experiment : ideal_2D_oce

purpose: Shows that the residual mean circulation is becoming small
 at the Equilibrium (Bolus cancel the Euleurian circ.) and is function
 of the diapycnal mixing.

Set up:
 Idealized 2D global ocean with flat bathymetry and no continent,
 symetric relative to the Eq.
 Forcing: zonal wind stress and surface temp. relaxation toward a
 "realistic" SST (function of Latitude).

To reduce diapycnal mixing, "exotic" parameters are used (and tested)
in this test-experiment:
a) vertical discretization: interface at the middle (use delRc)
b) GM advect form
c) use Visbeck
d) GM advect(Euler+Bolus) and Flux Limit Adv scheme.
e) 3 different time-steps (MOM,FS,Tracer)
f) oceanic exp using stagger time stepping
g) oceanic exp using cg2dTargetResWunit

```

## Build variants (code*/)
- **code**: packages.conf = `gfd cd_code gmredi diagnostics`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor cd_code gmredi diagnostics
  - SIZE.h: grid 1x56x15; sNx=1, sNy=14, OLx=3, OLy=3, nSx=1, nSy=4, nPx=1, nPy=1, Nr=15; has SIZE.h_mpi
  - option/size headers: DIAGNOSTICS_SIZE.h, GMREDI_OPTIONS.h

## Input variants (input*/)
- **input**: data.pkg on: useGMRedi, useDiagnostics
  - data: deltaTmom=1200.0, deltaTtracer=86400.0, nTimeSteps=20, nIter0=36000, eosType='LINEAR', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE., tempAdvScheme=77
  - namelist files: data data.diagnostics data.gmredi data.pkg eedata
- **input.geom**: data.pkg on: -
  - namelist files: data.diagnostics data.gmredi

## Reference results
`output.geom.txt` `output.txt`

Run: `cd verification; ./testreport -of <optfile> -t ideal_2D_oce` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
