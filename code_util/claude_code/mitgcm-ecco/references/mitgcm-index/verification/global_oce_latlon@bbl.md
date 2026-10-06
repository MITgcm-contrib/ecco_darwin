# verification/global_oce_latlon@bbl

**vs its upstream base (merge-base):** added: input.bblphys/data, input.bblphys/data.bbl, input.bblphys/data.cal, input.bblphys/data.diagnostics, input.bblphys/data.exf, input.bblphys/data.pkg, input.bblphys/eedata.mth, input.bblphys/prepare_run, results/output.bblphys.txt, run/ADJdiffkr.0000000000.001.001.data, run/ADJdiffkr.0000000000.001.001.meta, run/ADJdiffkr.0000000000.001.002.data, run/ADJdiffkr.0000000000.001.002.meta, run/ADJdiffkr.0000000000.002.001.data, run/ADJdiffkr.0000000000.002.001.meta, run/ADJdiffkr.0000000000.002.002.data, run/ADJdiffkr.0000000000.002.002.meta, run/ADJempr.0000000000.001.001.data, run/ADJempr.0000000000.001.001.meta, run/ADJempr.0000000000.001.002.data, run/ADJempr.0000000000.001.002.meta, run/ADJempr.0000000000.002.001.data, run/ADJempr.0000000000.002.001.meta, run/ADJempr.0000000000.002.002.data, run/ADJempr.0000000000.002.002.meta, run/ADJetan.0000000000.001.001.data, run/ADJetan.0000000000.001.001.meta, run/ADJetan.0000000000.001.002.data, run/ADJetan.0000000000.001.002.meta, run/ADJetan.0000000000.002.001.data, run/ADJetan.0000000000.002.001.meta, run/ADJetan.0000000000.002.002.data, run/ADJetan.0000000000.002.002.meta, run/ADJggl90tke.0000000000.001.001.data, run/ADJggl90tke.0000000000.001.001.meta, run/ADJggl90tke.0000000000.001.002.data, run/ADJggl90tke.0000000000.001.002.meta, run/ADJggl90tke.0000000000.002.001.data, run/ADJggl90tke.0000000000.002.001.meta, run/ADJggl90tke.0000000000.002.002.data; changed: README.md, results/output.yearly.txt; removed: build/.gitignore

## README (first 40 lines)
```
# Global Ocean Simulation at 4 degree Resolution, including Adjoint Set-Up
First configuration to use OpenAD, started on 2005-08-19
 by heimbach@mit.edu, utke@mcs.anl.gov, cnh@mit.edu
***Note*** this experiment was previously named "OpenAD".

### Overview:
This experiment is derived from `tutorial_global_oce_latlon` (see also `global_ocean.90x40x15`)
with surface forcing provided by specific pkgs, either
`pkg/exf` or `pkg/ebm`, instead of relying on the main model surface forcing capability.
It contains 4 forward set-up, all using the same executable built from `code`
config but with specific input files from `input/` (primary test) and,
as secondary tests, from `input.yearly/`, `input.bblphys/` and `input.ebm/`.

It provides also adjoint settings for 2 AD compilers, TAF and Tapenade, with primary
test input files in `input_ad/` and `input_tap/` respectively, but also
several secondary test setting with each of the AD compilers.

## Part 1, Forward only tests:

The **primary** forward test, using input files from `input/`,
uses prescribed monthly-mean air-sea surface fluxes from `pkg/exf`.<br>

The **secondary** test, using input files from `input.yearly/`,
is very similar except for the specification of yearly input fields to`pkg/exf`.<br>

**Note:**

1.  These 2 set-up have been moved (in PR #830) from `verification/global_with_exf/`
    where a "README" still provides some details related to `pkg/exf` specific
    features used here.
2.  The ability, using `pkg/exf`, to compute surface fluxes from near surface
    atmospheric state and downward radiation as shown, e.g., in experiment
    `global_ocean.cs32x15` (secondary test `input.seaice` or `input.icedyn` or `input.in_p`)
    is not used here (`#undef ALLOW_BULKFORMULAE`).

The **secondary** test, using input files from `input.bblphys/`,
is the same as `input.yearly/` (which uses `pkg/bbl` with a constant lateral
BBL speed) but switches on the optional `pkg/bbl` physics: a lateral BBL
speed computed from the density contrast, slope and Coriolis parameter
(`bbl_useNofSpeed`) and Richardson-number entrainment into the BBL
... (truncated)
```

## Build variants (code*/)
- **code**: packages.conf = `gfd cd_code gmredi bbl ebm exf frazil profiles diagnostics`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor cd_code gmredi bbl ebm exf frazil profiles diagnostics
  - SIZE.h: grid 90x40x15; sNx=30, sNy=20, OLx=2, OLy=2, nSx=3, nSy=2, nPx=1, nPy=1, Nr=15; has SIZE.h_mpi
  - option/size headers: DIAGNOSTICS_SIZE.h, EXF_OPTIONS.h, PROFILES_OPTIONS.h
- **code_ad**: packages.conf = `gfd cd_code gmredi ggl90 exf mnc ebm adjoint`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor cd_code gmredi ggl90 exf mnc ebm autodiff cost ctrl grdchk
  - SIZE.h: grid 90x40x15; sNx=45, sNy=20, OLx=2, OLy=2, nSx=2, nSy=2, nPx=1, nPy=1, Nr=15; has SIZE.h_mpi
  - option/size headers: AUTODIFF_OPTIONS.h, COST_OPTIONS.h, CPP_OPTIONS.h, CTRL_OPTIONS.h, CTRL_SIZE.h, EXF_OPTIONS.h, GGL90_OPTIONS.h, GMREDI_OPTIONS.h
  - modified/extra source: tamc.h
- **code_tap**: packages.conf = `gfd cd_code gmredi ggl90 exf tapenade adjoint`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor cd_code gmredi ggl90 exf tapenade autodiff cost ctrl grdchk
  - SIZE.h: grid 90x40x15; sNx=45, sNy=20, OLx=2, OLy=2, nSx=2, nSy=2, nPx=1, nPy=1, Nr=15; has SIZE.h_mpi
  - option/size headers: AUTODIFF_OPTIONS.h, COST_OPTIONS.h, CPP_OPTIONS.h, CTRL_OPTIONS.h, CTRL_SIZE.h, EXF_OPTIONS.h, GGL90_OPTIONS.h, GMREDI_OPTIONS.h

## Input variants (input*/)
- **input**: data.pkg on: useEXF, useCAL, useGMRedi, useDiagnostics, usePROFILES
  - data: deltaTmom=1200., deltaTtracer=43200., nTimeSteps=20, nIter0=0, eosType='POLY3', usingSphericalPolarGrid=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.cal data.diagnostics data.exf data.gmredi data.pkg data.profiles eedata prepare_run
- **input.bblphys**: data.pkg on: useEXF, useCAL, useGMRedi, useFRAZIL, useBBL, useDiagnostics
  - data: deltaTmom=1200., deltaTtracer=43200., nTimeSteps=20, nIter0=0, eosType='POLY3', usingSphericalPolarGrid=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.bbl data.cal data.diagnostics data.exf data.pkg prepare_run
- **input.ebm**: data.pkg on: useGMRedi, useEBM
  - data: deltaTmom=1200., deltaTtracer=43200., nTimeSteps=20, nIter0=0, eosType='POLY3', usingSphericalPolarGrid=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.ebm data.gmredi data.pkg eedata prepare_run
- **input.yearly**: data.pkg on: useEXF, useCAL, useGMRedi, useFRAZIL, useBBL, useDiagnostics
  - data: deltaTmom=1200., deltaTtracer=43200., nTimeSteps=20, nIter0=0, eosType='POLY3', usingSphericalPolarGrid=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.bbl data.cal data.diagnostics data.exf data.pkg prepare_run
- **input_ad**: data.pkg on: useGMRedi, useGrdchk
  - data: deltaTmom=1200., deltaTtracer=43200., nTimeSteps=4, nIter0=0, eosType='JMD95Z', usingSphericalPolarGrid=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.autodiff data.cost data.ctrl data.gmredi data.grdchk data.optim data.pkg eedata prepare_run
- **input_ad.ebm**: data.pkg on: useGMRedi, useEBM, useGrdchk
  - data: deltaTmom=1200., deltaTtracer=43200., nTimeSteps=20, nIter0=0, eosType='POLY3', usingSphericalPolarGrid=.TRUE.
  - namelist files: data data.autodiff data.cost data.ctrl data.ebm data.gmredi data.grdchk data.optim data.pkg eedata prepare_run
- **input_ad.ggl90**: data.pkg on: useGMRedi, useGGL90, useGrdchk
  - data: deltaTmom=1200., deltaTtracer=43200., nTimeSteps=4, nIter0=0, eosType='JMD95Z', usingSphericalPolarGrid=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.ggl90 data.grdchk data.pkg
- **input_ad.w_exf**: data.pkg on: useEXF, useCAL, useGMRedi, useGrdchk
  - data: deltaTmom=1200., deltaTtracer=43200., nTimeSteps=20, nIter0=0, eosType='JMD95Z', usingSphericalPolarGrid=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.autodiff data.cal data.cost data.ctrl data.exf data.gmredi data.grdchk data.optim data.pkg eedata
- **input_tap**: data.pkg on: -
  - data: deltaTmom=1200., deltaTtracer=43200., nTimeSteps=4, nIter0=0, eosType='JMD95Z', usingSphericalPolarGrid=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.cost data.ctrl prepare_run
- **input_tap.w_exf**: data.pkg on: -
  - data: deltaTmom=1200., deltaTtracer=43200., nTimeSteps=20, nIter0=0, eosType='JMD95Z', usingSphericalPolarGrid=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.cost data.ctrl prepare_run

## Reference results
`README` `output.bblphys.txt` `output.ebm.txt` `output.txt` `output.yearly.txt` `output_adm.ebm.txt` `output_adm.ggl90.txt` `output_adm.txt` `output_adm.w_exf.txt` `output_tap_adj.txt` `output_tap_adj.w_exf.txt` `output_tap_tlm.txt` `output_tap_tlm.w_exf.txt` `output_tlm.ggl90.txt.gz` `output_tlm.txt.gz`

Run: `cd verification; ./testreport -of <optfile> -t global_oce_latlon` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
