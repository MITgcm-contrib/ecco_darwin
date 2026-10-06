# verification/exp4

## README (first 40 lines)
```
Example: "Flow over a bump with Open Boundaries and passive tracers"
====================================================================

This experiment is a 400 x 210 x 4.5 km channel (80 x 42 x 8 grid-points)
on an f-plane with a tall seamount in the middle.
It is intended to illustrate the interaction of a zonal flow over and around
a steep bump as well as the use of Open Boundary (`pkg/obcs`) and Floats
packages (`pkg/flt`).

### Overview:
This experiment contains 4 set-ups (with corresponding `input[.*]/` dir) that
can be run with the same executable (built from `build/` dir using customized
code from `code/`);
binary input files (all `real*8`) have been generated using matlab script
`gendata.m` from `input` dir.
All four set-ups use a simple EOS ( $\rho' = -\rho_0 ~ \alpha_T ~ \theta'$ )
and treat salt as a passive tracer ;

The **primary** test, using input files from `input/` dir, use four open
boundaries with simple specifications (`useOBCSprescribe`) from open boundary
parameter file `data.obcs`.
This is a non-hydrostatic set-up using the flux-form momentum equations.

Different kinds of open boundary values are used:
zonal (x-)velocity U is prescribed at all open boundaries with values that are
read from data files (specified in data.obcs);
meridional (y-)velocity V is set to zero on all boundaries, and temperature to
`tRef(z)`, both `in obcs_calc.F`, this is the default behavior;
at the western boundary, salinity values are used for salinity and one passive
tracer in the same way.
Salinity is set to sLev at all other boundaries, while a (nearly) homogeneous
Neumann condition is applied to the passive tracer (the latter is the default
in `obcs_calc.F`), with a relaxation (using pkg rbcs) in the Eastern part of
the channel.

The **secondary** test, using `input.nlfs/` dir, is similar to the primary test except
it is hydrostatic with the vector-invariant momentum formulation and using the
`z*` coordinate. A time-varying small imbalance between the Western boundary
inflow and Eastern boundary outflow generates sea-level fluctuations.

... (truncated)
```

## Build variants (code*/)
- **code**: packages.conf = `gfd obcs rbcs ptracers layers flt`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor obcs rbcs ptracers layers flt
  - SIZE.h: grid 80x42x8; sNx=40, sNy=21, OLx=3, OLy=3, nSx=2, nSy=2, nPx=1, nPy=1, Nr=8; has SIZE.h_mpi
  - option/size headers: CPP_OPTIONS.h, FLT_OPTIONS.h, OBCS_OPTIONS.h
  - modified/extra source: CPP_EEOPTIONS.h_mpi

## Input variants (input*/)
- **input**: data.pkg on: useOBCS, usePtracers, useRBCS
  - data: deltaT=600.0, nTimeSteps=10, nIter0=0, eosType='LINEAR', implicitFreeSurface=.TRUE., nonHydrostatic=.TRUE., usingCartesianGrid=.TRUE., saltAdvScheme=4
  - namelist files: data data.obcs data.pkg data.ptracers data.rbcs eedata
- **input.nlfs**: data.pkg on: useOBCS, usePtracers, useRBCS
  - data: deltaT=600.0, nTimeSteps=10, nIter0=0, nonlinFreeSurf=4, select_rStar=2, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., saltAdvScheme=4
  - namelist files: data data.obcs data.pkg eedata
- **input.stevens**: data.pkg on: useOBCS
  - data: deltaT=600.0, nTimeSteps=10, nIter0=0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., saltAdvScheme=4
  - namelist files: data data.obcs data.pkg eedata
- **input.with_flt**: data.pkg on: useFLT
  - data: deltaT=600.0, nTimeSteps=18, nIter0=0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE.
  - namelist files: data data.flt data.pkg eedata

## Reference results
`output.nlfs.txt` `output.stevens.txt` `output.txt` `output.with_flt.txt`

Run: `cd verification; ./testreport -of <optfile> -t exp4` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
