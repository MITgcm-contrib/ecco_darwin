# verification/tutorial_barotropic_gyre

## README (first 40 lines)
```
Tutorial Example: "Barotropic gyre"
(Barotropic Ocean Gyre In Cartesian Coordinates)
================================================
(replaces old "exp0" verification ;
 also "barotropic_gyre_in_a_box" in release.1 branch)

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

There is comparison output in the directory:
  results/output.txt

Comments:
  The input data is real*4 and generated using the MATLAB script
  gendata.m.

```

## Build variants (code*/)
- **code**: packages.conf = `(inherits/none)`
  - SIZE.h: grid 62x62x1; sNx=62, sNy=62, OLx=2, OLy=2, nSx=1, nSy=1, nPx=1, nPy=1, Nr=1; has SIZE.h_mpi

## Input variants (input*/)
- **input**: data.pkg on: -
  - data: deltaT=1200.0, nTimeSteps=10, nIter0=0, implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE.
  - namelist files: data data.pkg eedata

## Reference results
`output.txt`

Run: `cd verification; ./testreport -of <optfile> -t tutorial_barotropic_gyre` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
