# verification/tutorial_plume_on_slope

## README (first 40 lines)
```
Tutorial Example: "Gravity plume on a continental slope"
========================================================
(formerly "plume_on_slope" verification ;
 also "nonhydrostatic_plume_on_slope" in release.1 branch)

### Overview:
This is a 2D set-up with (variable) high-resolution and non-hydrostatic dynamics, where dense water is produced on a shelf that then flows as a gravity current down the slope. 

The **primary** test uses a no-slip bottom boundary condition (`no_slip_bottom=.TRUE.`) and no explicit drag. 

The **secondary** test `rough.Bot` uses the logarithmic law of the wall to compute the drag coefficient for quadratic bottom drag as a function of distance from the bottom (i.e. cell thickness) and a prescribed roughness length `zRoughBot = 0.01` (in meters). For this configuration (i.e. vertical grid spacing) this value of `zRoughBot` corresponds to approximately `bottomDragQuadratic=5.E-2`. For consistency, the bottom boundary conditions is set to free slip (`no_slip_bottom=.FALSE.`).

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
  cd run
  ln -s ../input/* .
  ln -s ../build/mitgcmuv .
  ./mitgcmuv > output.txt
  cd ..
```

To run the **secondary** test `roughBot`:

```
  cd run
  rm *
  ln -s ../input.roughBot/* .
  ln -s ../input/* .
... (truncated)
```

## Build variants (code*/)
- **code**: packages.conf = `gfd obcs`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor obcs
  - SIZE.h: grid 320x1x60; sNx=80, sNy=1, OLx=3, OLy=3, nSx=4, nSy=1, nPx=1, nPy=1, Nr=60; has SIZE.h_mpi
  - option/size headers: CPP_OPTIONS.h, MOM_COMMON_OPTIONS.h, OBCS_OPTIONS.h

## Input variants (input*/)
- **input**: data.pkg on: useOBCS
  - data: deltaT=20.0, nTimeSteps=20, nIter0=0, eosType='LINEAR', implicitFreeSurface=.TRUE., nonHydrostatic=.TRUE., usingCartesianGrid=.TRUE., tempAdvScheme=33
  - namelist files: data data.obcs data.pkg eedata
- **input.roughBot**: data.pkg on: -
  - data: deltaT=20.0, nTimeSteps=20, nIter0=0, eosType='LINEAR', implicitFreeSurface=.TRUE., nonHydrostatic=.TRUE., usingCartesianGrid=.TRUE., tempAdvScheme=33
  - namelist files: data

## Reference results
`output.roughBot.txt` `output.txt`

Run: `cd verification; ./testreport -of <optfile> -t tutorial_plume_on_slope` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
