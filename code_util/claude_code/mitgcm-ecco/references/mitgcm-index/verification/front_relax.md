# verification/front_relax

## README (first 40 lines)
```
# Relaxation of a front in a channel : simplest example that uses GM-Redi parameterization


A 2-D, y-z set-up is used to mimic a zonally symmetric, reentrant channel with
a baroclinicly unstable initial density front.<br>
As meso-scale eddies are not resolved in this 2-D set-up, the GM-Redi
parameterization is used to represent their effects.

### Overview:
This experiment contains 5 set-ups (with corresponding `input[.*]/` dir) that
can be run with the same executable (built from `build/` dir using customized
code from `code/`); binary input files have been generated using matlab script
`gendata.m` from the corresponding `input` dir. All five set-ups use a
simple EOS ( $\rho' = -\rho_0 ~ \alpha_T ~ \theta'$ ) and treat salt as a
passive tracer ; without any surface forcing, the density front is expected to
flatten (GM effect) while salinity spread along isopycnal (Redi diffusion).

The **primary** test, using input files from `input/` dir, is the simplest
one, with flat bottom, non-uniform resolution in both direction (15 levels
from 50 m to 400 m thick near the bottom and, in Y-direction, 32 grid-points
with about 10 km spacing) and stratified every-where (background
$N = 2\times 10^{-3} ~s^{-1}$, see matlab script `input/gendata.m`),
avoiding the need for tapering or clipping.<br>
It uses the skew-flux formulation of GM with same Redi and GM diffusivity (
`GM_background_K` = 1000 $m^2/s$, see: `input/data.gmredi`). Note that 10 dead
levels were added (below the bottom) to allow to use the same executable
(compiled with `Nr = 25`) for all 5 set-ups.

The **secondary** test `input.in_p/` dir is the same as the primary test but
converted to use P-coordinates instead of height coordinates. For the purpose
of comparing P and Z coordinates, gravity and reference density `rhoNil` are
set to round number (respectively 10 and 1000) to facilitate conversions. It
uses the advective form of GM with same Redi and GM diffusivity (see:
`input.in_p/data.gmredi`).

The **next two secondary** set-ups, `input.mxl/` and `input.bvp/` are very
similar, sharing the same binary input files from `input.mxl/` dir ; they use
the full 25 level model to represent a 10 level, 200 m thick mixed layer on
top of a stratified warm bowl of water.  The `input.mxl/` illustrates the use
of the transition-layer tapering scheme 'fm07' with the skew-flux formulation
... (truncated)
```

## Build variants (code*/)
- **code**: packages.conf = `gfd gmredi diagnostics`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor gmredi diagnostics
  - SIZE.h: grid 1x32x25; sNx=1, sNy=16, OLx=3, OLy=3, nSx=1, nSy=2, nPx=1, nPy=1, Nr=25; has SIZE.h_mpi
  - option/size headers: DIAGNOSTICS_SIZE.h

## Input variants (input*/)
- **input**: data.pkg on: useGMRedi
  - data: deltaT=1800., nTimeSteps=20, startTime=0., eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., buoyancyRelation='OCEANIC', tempAdvScheme=20, saltAdvScheme=20
  - namelist files: data data.gmredi data.mpi data.pkg eedata
- **input.bvp**: data.pkg on: useGMRedi, useDiagnostics
  - data: deltaT=3600., nTimeSteps=25, startTime=0., eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE.
  - namelist files: data data.diagnostics data.gmredi data.pkg prepare_run
- **input.in_p**: data.pkg on: useGMRedi, useDiagnostics
  - data: deltaT=1800., nTimeSteps=20, startTime=0., eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., buoyancyRelation='OCEANICP', tempAdvScheme=20, saltAdvScheme=20
  - namelist files: data data.diagnostics data.gmredi data.pkg eedata
- **input.mxl**: data.pkg on: useGMRedi, useDiagnostics
  - data: deltaT=3600., nTimeSteps=25, startTime=0., eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE.
  - namelist files: data data.diagnostics data.gmredi data.pkg
- **input.top**: data.pkg on: useGMRedi, useDiagnostics
  - data: deltaT=3600., nTimeSteps=25, startTime=0., eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE.
  - namelist files: data data.diagnostics data.gmredi data.pkg

## Reference results
`output.bvp.txt` `output.in_p.txt` `output.mxl.txt` `output.top.txt` `output.txt`

Run: `cd verification; ./testreport -of <optfile> -t front_relax` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
