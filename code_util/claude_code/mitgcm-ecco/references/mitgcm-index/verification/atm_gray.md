# verification/atm_gray

## README (first 40 lines)
```
Gray atmosphere physics example on Cubed-Sphere grid
============================================================

Use gray atmospheric physics (O'Gorman and Schneider, JCl, 2008)
from package `atm_phys` inside  MITgcm dynamical core, in a global
cubed-sphere grid set-up (6 faces 32x32, 26 levels, non uniform deltaP).

### Overview:
This experiment contains 2 aqua-planet like set-ups (with corresponding `input[.*]/` dir)
that can be run with the same executable (built from `build/` dir using customized
code from `code/`); binary input files have been generated using matlab script
`gendata.m` from the `input` dir.
Both test experiments start from a spin-up state using pickup files written after 1 year.

The **primary** test, using input files from `input/` dir,
has an interactive SST with a 10m mixed layer depth and a prescribed,
time-invariant Q-flux. It also includes a weak damping of stratospheric winds.

The **secondary** test, using files from `input.ape/` dir, is the same as the primary
test but without stratospheric wind damping, and uses prescribed idealized SST from
Neale and Hoskins, 2001, Aqua-Planet Experiment (APE) project.

### Instructions:
Configure and compile the code:

```
  cd build
  ../../../tools/genmake2 -mods ../code [-of my_platform_optionFile]
  make depend
  make
  cd ..
```

To run primary test:

```
  cd run
  ln -s ../input/* .
  ./prepare_run
  ../build/mitgcmuv > output.txt
... (truncated)
```

## Build variants (code*/)
- **code**: packages.conf = `exch2 gfd shap_filt atm_phys diagnostics`
  - expanded: exch2 mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor shap_filt atm_phys diagnostics
  - SIZE.h: grid 192x32x26; sNx=32, sNy=32, OLx=4, OLy=4, nSx=6, nSy=1, nPx=1, nPy=1, Nr=26; has SIZE.h_mpi
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h

## Input variants (input*/)
- **input**: data.pkg on: useSHAP_FILT, useDiagnostics, useAtm_Phys
  - data: deltaT=384., nTimeSteps=10, nIter0=81000, nonlinFreeSurf=4, select_rStar=2, eosType='IDEALG', implicitFreeSurface=.TRUE., usingCurvilinearGrid=.TRUE., buoyancyRelation='ATMOSPHERIC', saltAdvScheme=77
  - namelist files: data data.atm_gray data.atm_phys data.diagnostics data.pkg data.shap eedata prepare_run
- **input.ape**: data.pkg on: -
  - namelist files: data.atm_gray data.atm_phys data.diagnostics eedata

## Reference results
`output.ape.txt` `output.txt`

Run: `cd verification; ./testreport -of <optfile> -t atm_gray` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
