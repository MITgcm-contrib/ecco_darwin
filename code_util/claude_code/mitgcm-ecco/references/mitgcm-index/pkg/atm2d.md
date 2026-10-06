# pkg/atm2d

2-D (zonally averaged) atmosphere coupled to an MITgcm ocean (IGSM-style climate coupling).

**runtime switch:** `useATM2D`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.atm2d`

## Namelist parameters
### COUPLE_PARM
- `dtcouple`
- `dtatm`
- `dtocn`
- `startYear`
- `endYear`
- `taveDump`
### PARM01_ATM2D
- `atmosTauuFile`
- `atmosTauvFile`
- `atmosWindFile`
### PARM02_ATM2D
- `tauuFile`
- `tauvFile`
- `windFile`
- `qnetFile`
- `evapFile`
- `precipFile`
### PARM03_ATM2D
- `thetaRelaxFile`
- `saltRelaxFile`
- `tauThetaRelax`
- `tauSaltRelax`
- `nttyperelax`
- `nstyperelax`
### PARM04_ATM2D
- `runoffFile`
- `runoffMapFile`
- `numbands`
- `rband`
### PARM05_ATM2D
- `useObsEmP`
- `useObsRunoff`
- `useAltDeriv`

## CPP options (defaults as shipped)
- `ATM2D_MPI_ON` (undef, ATM2D_OPTIONS.h) — turn on MPI or not
- `ALLOW_DBUG_ATM2D` (define, ATM2D_OPTIONS.h) — - allow single grid-point debugging write to standard-output
- `JBUGI` (define, ATM2D_OPTIONS.h)
- `JBUGJ` (define, ATM2D_OPTIONS.h)
- `CLM` (undef, ATM2D_OPTIONS.h) — - undocumented options:
- `CLM35` (undef, ATM2D_OPTIONS.h)
- `CPL_CHEM` (undef, ATM2D_OPTIONS.h)
- `CPL_NEM` (undef, ATM2D_OPTIONS.h)
- `CPL_OCEANCO2` (undef, ATM2D_OPTIONS.h)
- `CPL_TEM` (undef, ATM2D_OPTIONS.h)
- `DATA4TEM` (undef, ATM2D_OPTIONS.h)
- `IPCC_EMI` (undef, ATM2D_OPTIONS.h)
- `ML_2D` (undef, ATM2D_OPTIONS.h)
- `NCEPWIND` (undef, ATM2D_OPTIONS.h)
- `OCEAN_3D` (undef, ATM2D_OPTIONS.h)
- `ORBITAL_FOR` (undef, ATM2D_OPTIONS.h)
- `PREDICTED_AEROSOL` (undef, ATM2D_OPTIONS.h)

## Headers
- `AGRID.h` — 
- `ATM2D_OPTIONS.h` — Package-specific Options & Macros go here
- `ATM2D_VARS.h` — Files: mean 2D atmos fields used for wind anomaly coupling
- `ATMSIZE.h` — **** ATM 2D size definitions. Also declared in BD2G04.COM constants N_LAT and N_LEV declared in ctrparam.h
- `CPLIDS.h` — /==========================================================\ are used to identify this component and the fields it exchanges with other components. No
- `DRIVER.h` — 
- `OCNIDS.h` — are used to identify this component and the fields it exchanges with other components.
- `OCNSIZE.h` — /==========================================================\ OCN_SIZE.h Declare size of underlying computational grid for ocean component. \==========
- `OCNVARS.h` — grid. Arrays may need adding or removing different couplings.
- `ctrparam.h` — Purpose:     A header file contains cpp control parameters for the model Usage:       1. (un)comment #define line; or 2. #define/#undef x [number] to 

## Routines (31)
`atm2d_finish.F`, `atm2d_init_fixed.F`, `atm2d_init_vars.F`, `atm2d_read_pickup.F`, `atm2d_readparms.F`, `atm2d_write_pickup.F`, `atm2ocn_main.F`, `calc_1dto2d.F`, `calc_fileload.F`, `calc_zonal_means.F`, `fixed_flux_add.F`, `forward_step_atm2d.F`, `get_ocnvars.F`, `init_atm2d.F`, `init_sumvars.F`, `month_end_diags.F`, `norm_ocn_fluxes.F`, `pass_thsice_fluxes.F`, `put_ocnvars.F`, `read_atmos.F`, `relax_add.F`, `subtract_means.F`, `sum_ocn_fluxes.F`, `sum_thsice_out.F`, `sum_yr_end_diags.F`, `tave_end_diags.F`, `yr_end_diags.F`

## Called from outside the package
- `FORWARD_STEP_ATM2D` ← `model/src/main_do_loop.F:217`
- `ATM2D_INIT_FIXED` ← `model/src/packages_init_fixed.F:551`
- `ATM2D_INIT_VARS` ← `model/src/packages_init_variables.F:472`
