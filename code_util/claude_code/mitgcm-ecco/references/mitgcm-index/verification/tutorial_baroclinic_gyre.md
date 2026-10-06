# verification/tutorial_baroclinic_gyre

## README (first 40 lines)
```
Tutorial Example: "Baroclinic gyre"
(Baroclinic Ocean Gyre In Spherical Coordinates)
============================================================
(formerly "exp1" verification ;
 also "baroclinic_gyre_on_a_sphere" in release.1 branch)

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
  ../build/mitgcmuv > output.txt
```

There is comparison output in the directory:
  results/output.txt

Comments:
  The input data is real*4 and generated using the MATLAB script
  gendata.m.
```

## Build variants (code*/)
- **code**: packages.conf = `gfd diagnostics mnc`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor diagnostics mnc
  - SIZE.h: grid 62x62x15; sNx=31, sNy=31, OLx=2, OLy=2, nSx=2, nSy=2, nPx=1, nPy=1, Nr=15; has SIZE.h_mpi
  - option/size headers: DIAGNOSTICS_SIZE.h

## Input variants (input*/)
- **input**: data.pkg on: useMNC, useDiagnostics
  - data: deltaT=1200., endTime=12000., startTime=0., eosType='LINEAR', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE.
  - namelist files: data data.diagnostics data.mnc data.pkg eedata

## Reference results
`output.txt`

Run: `cd verification; ./testreport -of <optfile> -t tutorial_baroclinic_gyre` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
