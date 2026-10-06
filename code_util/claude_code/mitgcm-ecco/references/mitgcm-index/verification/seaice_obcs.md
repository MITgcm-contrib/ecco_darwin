# verification/seaice_obcs

## README (first 40 lines)
```
Test set-up for seaice pkg with Open-Boundary Conditions
========================================================

This verification experiment is used to test `pkg/seaice` with `pkg/obcs`
and the set-up itself is carved out from `../lab_sea/input.salt_plume/`.

The **primary** test uses input files from `input/` dir which have been generated using
the matlab script [`input/mk_input.m`](input/mk_input.m) together with the set of output files from running
test experiment `../lab_sea/input.salt_plume`.

The **secondary** test `input.seaiceSponge/` uses OBCS sponge-layer for seaice fields
(`useSeaiceSponge=.TRUE.`) in addition to prescribed OBCS from the primary test.

The **secondary** test `input.tides/` adds 4 tidal components to the barotropic velocity
at the open-boundaries (`useOBCStides=.TRUE.`) in addition to prescribed OBCS from the
primary test. The additional tidal component binary input files have been generated using
the matlab script `mk_tides.m` (see comments inside).

Note: naming of tidal input files ("OB\*File") has been changed and augmented to include
Open-Boundary tangential flows in PR [#752](https://github.com/MITgcm/MITgcm/pull/752)
 (see [`input.tides/update_TideFileName.sed`](input.tides/update_TideFileName.sed) on how
to update `data.obcs`)
```

## Build variants (code*/)
- **code**: packages.conf = `oceanic obcs exf seaice salt_plume`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor gmredi kpp obcs exf seaice salt_plume
  - SIZE.h: grid 10x8x23; sNx=5, sNy=8, OLx=4, OLy=4, nSx=2, nSy=1, nPx=1, nPy=1, Nr=23; has SIZE.h_mpi
  - option/size headers: CPP_OPTIONS.h, GMREDI_OPTIONS.h, OBCS_OPTIONS.h, SEAICE_OPTIONS.h

## Input variants (input*/)
- **input**: data.pkg on: useGMRedi, useKPP, useEXF, useCAL, useSEAICE, useSALT_PLUME, useOBCS
  - data: deltaTmom=3600.0, deltaTtracer=3600.0, endTime=21600., startTime=3600.0, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE., tempAdvScheme=7, saltAdvScheme=7, momStepping=.TRUE.
  - namelist files: data data.cal data.exf data.gmredi data.kpp data.obcs data.pkg data.salt_plume data.seaice eedata prepare_run
- **input.regDenom**: data.pkg on: -
  - namelist files: data.seaice
- **input.seaiceSponge**: data.pkg on: -
  - namelist files: data.obcs
- **input.tides**: data.pkg on: -
  - data: deltaT=3600.0, endTime=21600., startTime=3600.0, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE., tempAdvScheme=7, saltAdvScheme=7
  - namelist files: data data.obcs

## Reference results
`output.regDenom.txt` `output.seaiceSponge.txt` `output.tides.txt` `output.txt`

Run: `cd verification; ./testreport -of <optfile> -t seaice_obcs` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
