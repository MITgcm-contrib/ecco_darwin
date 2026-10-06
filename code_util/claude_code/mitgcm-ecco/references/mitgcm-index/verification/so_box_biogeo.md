# verification/so_box_biogeo

## README (first 40 lines)
```
Southern-Ocean box with Biochemistry, using Open-Boundary Conditions
 (pkg/obcs) at Northern, Eastern and Western edges of the domain.
======================================================================

This experiment illustrates and tests the use of package DIC with OBCS.

The configuration (e.g., resolution), model parameters and forcing are
almost identical to tutorial_global_oce_biogeo expect that the horizontal
domain is limited to a sub-domain around Drake passage with open-boundary
conditions coming from the last year of a 2 yrs global simulation and
initial conditions taken from t=1.yr of this same simulation.
This enable to compare directly the results of the first year simulation
of this regional set-up with the 2nd year results of the global set-up run.

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

----------------------------------------------------------------------
To generate inital and open boundary conditions:
a) The global set-up (using executable from: tutorial_global_oce_biogeo/build)
   was run for 2 years using model-parameter from inp_global/.
 The only differences with the ones used in tutorial_global_oce_biogeo are:
 - CD-Scheme is turned off (useCDscheme=F), since it is not implemented for OBCS ;
 - as a consequence, horizontal viscosity (viscAh) is increased from 2.E5 to 3.E5;
 - convective-adjustment diffusivity (ivdc_kappa) is reduced from 100. to 10.m^2/s;
 - implicit vertical viscosity is turned off (not needed);
... (truncated)
```

## Build variants (code*/)
- **code**: packages.conf = `gfd cd_code obcs gmredi ptracers gchem dic diagnostics`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor cd_code obcs gmredi ptracers gchem dic diagnostics
  - SIZE.h: grid 42x20x15; sNx=14, sNy=10, OLx=3, OLy=3, nSx=3, nSy=2, nPx=1, nPy=1, Nr=15; has SIZE.h_mpi
  - option/size headers: DIAGNOSTICS_SIZE.h, DIC_OPTIONS.h, PTRACERS_SIZE.h

## Input variants (input*/)
- **input**: data.pkg on: usePTRACERS, useGCHEM, useGMRedi, useOBCS, useDiagnostics
  - data: deltaTmom=900., deltaTtracer=43200., nTimeSteps=10, nIter0=0, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE.
  - namelist files: data data.diagnostics data.dic data.gchem data.gmredi data.obcs data.pkg data.ptracers eedata
- **input.caSat0**: data.pkg on: -
  - namelist files: data.diagnostics data.dic eedata
- **input.caSat3**: data.pkg on: -
  - namelist files: data.diagnostics data.dic eedata
- **input.saphe**: data.pkg on: -
  - namelist files: data.dic eedata

## Reference results
`output.caSat0.txt` `output.caSat3.txt` `output.saphe.txt` `output.txt`

Run: `cd verification; ./testreport -of <optfile> -t so_box_biogeo` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
