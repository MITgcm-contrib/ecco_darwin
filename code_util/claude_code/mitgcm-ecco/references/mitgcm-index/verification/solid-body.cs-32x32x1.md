# verification/solid-body.cs-32x32x1

## README (first 40 lines)
```
Simple solid-body rotation test on cubed-sphere grid
========================================================

### Overview:
This is a single level, steady-state atmospheric example
(`bouyancyRelation='ATMOSPHERIC'`) on cubed-sphere (cs-32) grid with initial
zonal wind field $U(\phi)$ and surface pressure anomaly $\eta(\phi)$,
both dependent on latitude $\phi$ only, that corresponds to an additional
relative rotation ($\omega\'$) on top of the solid-planet rotation ($\Omega$)
and around the same axis:

$$ U(\phi) = U_{eq} ~ \cos( \phi ) ~~~ \mathrm{with:} ~~~ U_{eq} = \omega' \times R $$

$$ \eta(\phi) = \rho_{const} ~ U_{eq} ~ ( \Omega R + U_{eq} / 2 ) ~~ ( \cos^{2}(\phi) - 2/3 ) $$

The parameters used here are slightly different from Earth (an opportunity to
test this capability) with a smaller planet radius (`rSphere`) $R = 5500 km$,
a slower rotation ( 30 h period, `rotationPeriod=108000.`) and an equatorial
zonal wind $U_{eq} = 80 m/s$ which corresponds to a 5 day revolution time.

The set-up uses linear free-surface with uniform density $\rho_{const} = 1$,
no viscosity and no bottom friction so that the solution is expected to remain
unchanged over time.
A bell-shape patch of passive tracer centered at mid-latitude
($\phi_{0} = 45^{o}$) is advected with the simulated wind field.

## Instructions
Configure and compile the code:

```
  cd build
  ../../../tools/genmake2 -mods ../code [-of my_platform_optionFile]
  make depend
  make
  cd ..
```

To run:

```
... (truncated)
```

## Build variants (code*/)
- **code**: packages.conf = `exch2 gfd -mom_fluxform diagnostics`
  - expanded: exch2 mom_common mom_vecinv generic_advdiff debug mdsio rw monitor diagnostics
  - SIZE.h: grid 192x32x1; sNx=32, sNy=32, OLx=2, OLy=2, nSx=6, nSy=1, nPx=1, nPy=1, Nr=1; has SIZE.h_mpi
  - option/size headers: DIAGNOSTICS_SIZE.h
  - modified/extra source: ini_psurf.F, ini_vel.F

## Input variants (input*/)
- **input**: data.pkg on: useDiagnostics
  - data: deltaT=450., nTimeSteps=25, nIter0=0, eosType='IDEALG', implicitFreeSurface=.TRUE., usingCurvilinearGrid=.TRUE., buoyancyRelation='ATMOSPHERIC'
  - namelist files: data data.diagnostics data.exch2 data.pkg eedata

## Reference results
`output.txt`

Run: `cd verification; ./testreport -of <optfile> -t solid-body.cs-32x32x1` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
