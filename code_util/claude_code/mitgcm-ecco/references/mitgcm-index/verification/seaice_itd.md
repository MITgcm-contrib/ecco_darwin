# verification/seaice_itd

## README (first 40 lines)
```
Seaice-only verification experiment in idealized periodic channel with
ice thickness distribution (otherwise very similar to
offline_exf_seaice): CPP-flag SEAICE_ITD is defined
-----------------------------------------------------------------

1) main forward experiment (code, input)

  Re-entrant zonally periodic channel (80x42 grid points) with just level (Nr=1)
   uniform resolution (5.km, 10m), solid Southern boundary with triangular shape
   coastline ("bathy_3c.bin")

  Use seaice (dynamics & thermodynamics from pkg/seaice) with EXF (see data.pkg)
   with initial ice thickness ranging from nearly 0 m in the "south"
   to over 7 m in the "north"(but no snow)
   (HeffFile  = 'heff_quartic.bin', in "input/data.seaice")
  Initial seaice concentration is 100 % everywhere
   (AreaFile='const100.bin', in "input/data.seaice")
  and seaice is initially at rest.

  Ridging is computed according to Thorndyke et al (1975) and Hibler
  (1980). Ice strength P is computed following Rothrock (1975)

  At runtime turn off time-stepping in 'data', PARM01, using:
    momStepping  = .FALSE.,
    saltStepping = .FALSE.,
    tempAdvection=.FALSE.,

 Forcing:
  None of the forcing vary with time; the input files have been
   generated using the python script "input/gendata.py".
  SST relaxation field is uniform in X, parabolic function of Y with
   maximum close to Southern boundary.

  Atmospheric air temp is uniform in Y, and only vary with X (~sin(2.pi.x/Lx))
   with an amplitude of 4.K ('tair_4x.bin');
  Uses constant Relative Humidity (70%, file 'qa70_4x.bin')
  constant and uniform downward shortwave (100.W/m2, 'dsw_100.bin'),
                       downward longwave (250.W/m^2, 'dlw_250.bin'),
                       zonal wind (10.m/s, 'windx.bin'),
  no meridional wind, no precip.
... (truncated)
```

## Build variants (code*/)
- **code**: packages.conf = `gfd exf seaice diagnostics`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor exf seaice diagnostics
  - SIZE.h: grid 80x42x1; sNx=40, sNy=21, OLx=3, OLy=3, nSx=2, nSy=2, nPx=1, nPy=1, Nr=1; has SIZE.h_mpi
  - option/size headers: DIAGNOSTICS_SIZE.h, SEAICE_OPTIONS.h, SEAICE_SIZE.h
  - modified/extra source: EXCH.h

## Input variants (input*/)
- **input**: data.pkg on: useEXF, useSEAICE, useDiagnostics
  - data: deltaT=1800.0, nTimeSteps=12, startTime=0.0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., momStepping=.FALSE.
  - namelist files: data data.diagnostics data.exf data.pkg data.seaice eedata
- **input.lipscomb07**: data.pkg on: -
  - namelist files: data.exf data.seaice
- **input.thermo**: data.pkg on: -
  - namelist files: data.seaice

## Reference results
`output.lipscomb07.txt` `output.thermo.txt` `output.txt`

Run: `cd verification; ./testreport -of <optfile> -t seaice_itd` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
