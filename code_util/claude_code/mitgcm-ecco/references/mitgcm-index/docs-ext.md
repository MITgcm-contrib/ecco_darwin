# External docs map (ECCO v4, ecco_darwin, darwin3 manual)

Section headings with line numbers (`L123`) and the names each section cites, for docs outside MITgcm
proper. Plain-text readmes show their first line and the MITgcm checkpoint / optfile / grid / git
commands they mention. Read the real file at the path shown (roots listed below).


# @ecco_v4_configs = `~/.cache/mitgcm-ecco-sources/ecco-docs/ECCO-v4-Configurations` (8393df4 2026-09-24 Merge pull request #246 from owang01/first_r7_commit)

## `README.md`
- L1 ECCOv4-Configurations
  - L8 ECCO Version 4
    - L12 Latest: Release 4: 1992-2017
    - L19 Release 3: 1992-2015
    - L23 Release 2: 1992-2011
    - L31 Release 1: 1992-2011
  - L38 Other Directories
    - L40 devel:
  - L44 Support
  - L53 References
- `ECCOv4 Release 1/README` (8 lines) — For the configurations for ECCOv4 Release 1 and Release 2, please see the ECCOv4 github page:
- `ECCOv4 Release 2/README` (8 lines) — For the configurations for ECCOv4 Release 1 and Release 2, please see the ECCOv4 github page:
- `ECCOv4 Release 3/README` (46 lines) — ECCO Version 4: Third Release [ECCO v4-r3] [ftp://ecco.jpl.nasa.gov/Version4/Release3/] — `checkpoint65u`
- `ECCOv4 Release 3/code/README` (3 lines) — Files originally stored on the MItgcm_contrib repository:
- `ECCOv4 Release 3/doc/README` (50 lines) — ECCO Version 4: Third Release [ECCO v4-r3] [ftp://ecco.jpl.nasa.gov/Version4/Release3/] — `checkpoint65p`, `V4r3`
- `ECCOv4 Release 3/namelist/README` (47 lines) — ECCO Version 4: Third Release [ECCO v4-r3] [ftp://ecco.jpl.nasa.gov/Version4/Release3/] — `checkpoint65u`, `LLC90`
- `ECCOv4 Release 4/README` (45 lines) — ECCO Version 4: Fourth Release [ECCO v4-r4] [https://ecco.jpl.nasa.gov/drive/files/Version4/Release4/] — `checkpoint66g`, `V4r4`

## `ECCOv4 Release 4/Docker/README.md`
  - L14 Directions
- `ECCOv4 Release 4/flux-forced/README` (39 lines) — ECCO Version 4: Fourth Release [ECCO v4-r4] [https://ecco.jpl.nasa.gov/drive/files/Version4/Release4/] — `checkpoint66g`

## `ECCOv4 Release 4/flux-forced/doc/README_fluxforced.md`
- L1 Configuration for flux-forced version of ECCO Version 4 Release 4
  - L5 Introduction
  - L35 Code
  - L49 Updated namelists
  - L63 Forcing
  - L100 Control variables
  - L129 Sample namelist to define a cost function

## `ECCOv4 Release 4/flux-forced/doc/README_offline_ptracer.md`
- L1 Offline passive tracer
  - L4 Overview
  - L13 Forward passive tracer — `code_offline_ptracer`, `namelist_offline_ptracer`, `state_weekly`
  - L24 Adjoint passive tracer — `uVeltave.0000227808.data`, `uVeltave.0000000000.data`, `reverseintime_all.sh`, `state_weekly`, `state_weekly_rev_time_227808`, `state_weekly_rev_time_8904`
  - L33 References:
- `ECCOv4 Release 4/optimization/README` (73 lines) — December 18, 2019
- `ECCOv4 Release 4/optimization/lsopt/README` (46 lines) — > Obtaining optimized BLAS routines for HPC platforms
- `ECCOv4 Release 4/optimization/optim/README` (70 lines) — c expid - experiment name
- `ECCOv4 Release 4/scripts/optim/README` (22 lines) — Example directory "optim" for the variable "optimdir" in
- `ECCOv4 Release 5/README` (38 lines) — ECCO Version 4: Fifth Release [ECCO V4r5] — `checkpoint68g`, `V4r5`

## `ECCOv4 Release 5/extension/code/slr_corr/README.md`
- L1 Overview of `slr_corr` — `slr_corr`
  - L4 slr_corr_readparms.F — `data.slr_corr`, `SLR_CORR_PARAM.h`
  - L7 slr_corr_init_fixed.F — `slrc_obs_timeseries`, `SLR_CORR_FIELDS.h`
  - L10 slr_corr_init_varia.F — `slr_corr_adjust_precip.F`, `slr_corr_init_varia.F`
  - L13 slr_corr_adjust_precip.F — `slrc_obs_timeseries`
- `ECCOv4 Release 6/README` (36 lines) — ECCO Version 4: Fifth Release [ECCO V4r6] — `checkpoint68g`, `V4r6`
- `ECCOv4 Release 7/README` (36 lines) — ECCO Version 4: Fifth Release [ECCO V4r6] — `checkpoint68g`, `V4r6`
- `devel/README` (19 lines) — This directory contains pieces of code that are useful, but not checked into the main MITgcm repository.
- `devel/V4r5_devel/README` (16 lines) — This directory contains the development version of the ECCO Version 4, Release 5 configuration. — `checkpoint68g`

## `devel/fluxforced/README.md`
- L1 Configurations for flux-forced runs
- `devel/iceshelf_icefront/README` (14 lines) — This directory contains the merged ice-shelf and ice-front code.

## `devel/reproduction/README.md`
- `devel/seaice_adjoint/README` (13 lines) — This directory contains code that has Ian Fenty's adjointable sea-ice code.
- `devel/seaice_adjoint/code/README` (24 lines) — /nobackupp7/owang/v4_release2/FROM_CVS_NEWEST_20170606/MITgcm/verification/release3/code06

# @ecco_v4_python = `~/.cache/mitgcm-ecco-sources/ecco-docs/ECCO-v4-Python-Tutorial` (3f0fcca 2026-05-15 Merge pull request #115 from andrewdelman/Steric_SSH_OBP)

## `README.md`
- L1 ECCO Version 4 Python Tutorial

## `Cloud_Setup/JPL_setup_instructions.md`
- L1 JPL setup for AWS EC2 instances
  - L7 Step 2: Start a JPL EC2 instance — `power_user`
    - L31 Step 3a: Enable ssh access — `sshd_enable.sh`
    - L62 Step 3b: Set up conda environment — `jupyter_env_setup.sh`

## `Docker/README.md`
- L1 Run tutorials on an AWS EC2 instance using a Docker container
  - L5 Getting started — `jupyter_env_setup.sh`
  - L9 Build the Docker image
  - L15 Run the Docker image — `ServerApp`, `LabApp`
  - L29 Open Jupyter lab in your browser — `Tutorials_as_Jupyter_Notebooks`
  - L45 Re-connect to Jupyter lab in Docker container — `ServerApp`, `LabApp`

## `Docker/smce/README.md`
- L1 Run tutorials on an AWS EC2 instance using a SMCE Docker image
  - L6 Getting started
  - L14 Build the Docker image
  - L24 Run the Docker image — `ServerApp`, `LabApp`
  - L44 Open Jupyter lab in your browser — `Tutorials_as_Jupyter_Notebooks`
  - L60 Re-connect to Jupyter lab in Docker container — `ServerApp`, `LabApp`

## `Intro_to_PO_Tutorials/Intro_to_PO_start.rst`
- L6 What are the Intro to PO Tutorials?
- L13 Who are these tutorials for?
- L21 What do I need to get started?
- L40 Which concepts are covered?
- `Tutorials_as_Jupyter_Notebooks/README.txt` (49 lines) — ECCO Version 4 Tutorial Juypter Notebooks Filenames

## `api_one_day/api.rst`
- L4 ecco_v4_py API reference
  - L8 Reading General Binary Files
  - L14 Reading LLC MDS (binary) Files
  - L20 Reading ECCO netcdf fields
  - L26 Convert LLC arrays to/from compact, face, and tile
  - L32 Plotting LLC fields
  - L38 Regridding fields to Lat-Lon grids
  - L44 Plotting fields on different projections
  - L50 Rotating LLC vector fields
  - L56 Exchanging/expanding values along LLC tile boundaries
  - L62 Other ECCO utility functions
  - L68 ECCO netcdf generation from MITgcm output

## `api_one_day/source/ecco_v4_py.rst`
- L1 ecco\_v4\_py package
  - L4 Submodules
  - L7 ecco\_v4\_py.ecco\_utils module
  - L15 ecco\_v4\_py.llc\_array\_conversion module
  - L23 ecco\_v4\_py.netcdf\_product\_generation module
  - L31 ecco\_v4\_py.read\_bin\_gen module
  - L39 ecco\_v4\_py.read\_bin\_llc module
  - L47 ecco\_v4\_py.resample\_to\_latlon module
  - L55 ecco\_v4\_py.test\_llc\_array\_loading\_and\_conversion module
  - L63 ecco\_v4\_py.tile\_exchange module
  - L71 ecco\_v4\_py.tile\_io module
  - L79 ecco\_v4\_py.tile\_plot module
  - L87 ecco\_v4\_py.tile\_plot\_proj module
  - L95 ecco\_v4\_py.tile\_rotation module
  - L104 Module contents

## `api_one_day/source/modules.rst`
- L1 ECCOv4-py

## `api_one_day/source/setup.rst`
- L1 setup module

## `doc/Installing_Python_and_Python_Packages.rst`
- L2 Python and Python Packages
  - L9 Why Python?
  - L25 Installing Python
    - L30 Anaconda
  - L49 Downloading the *ecco_v4_py* Python Package
    - L62 Option 1: Clone into the repository using git (recommended)
    - L80 Option 2: Download the repository using git (less recommended)
    - L92 Option 3: Use the *conda* package manager (less recommended)
    - L109 Option 4: Use the *pip* package manager (not at all recommended)
  - L119 Installing Dependencies
  - L133 Using the *ecco_v4_py* in your programs

## `doc/Intro_to_PO_start.rst`
- L6 What are the Intro to PO Tutorials?
- L13 Who are these tutorials for?
- L21 What do I need to get started?
- L40 Which concepts are covered?

## `doc/Tutorial_Introduction.rst`
- L2 Tutorial Overview
  - L6 What is the format of the tutorials?
    - L18 What are Jupyter notebooks?
  - L40 Will I learn Python just from reading these tutorials?
  - L46 What Python should I review before getting started?
    - L51 NumPy
    - L62 Matplotlib
    - L72 xarray
  - L83 What if I don't like the way you do X?
  - L89 What if I find a mistake?
  - L95 What if I would like to contribute with a tutorial of my own?
  - L101 Bonus Tutorials

## `doc/Tutorial_wget_Command_Line_HTTPS_Downloading_ECCO_Datasets_from_PODAAC.md`
- L1 Using _wget_ to Download ECCO Datasets from PO.DAAC
  - L7 Step 1: Create an account with NASA Earthdata
  - L15 Step 2: Set up your ```netrc``` and ```urs_cookies``` files — `_netrc`, `urs_cookies`
  - L47 Step 3: Prepare a list of granules (files) to download
  - L79 Step 4: Download files in a batch with GNU *_wget_*

## `doc/fields.rst`
- L2 ECCO v4 state estimate ocean, sea-ice, and atmosphere fields
  - L10 Geographical layout
    - L18 13-tile *native* lat-lon-cap 90 grid
      - L33 Available fields on the llc90 grid
    - L43 *interpolated* 0.5° x 0.5° latitude-longitude grid
      - L48 Available fields on the 0.5° x 0.5° latitude-longitude grid
    - L56 Miscellaneous fields and data
  - L65 Temporal frequency of state estimate fields
  - L71 Custom output

## `doc/index.rst`
- L6 Welcome to the ECCO Version 4 Tutorial
  - L12 Additional Resources
- L102 Indices and tables

## `doc/intro.rst`
- L2 The ECCO Ocean and Sea-Ice State Estimate
  - L6 What is the ECCO Central Production State Estimate?
    - L23 Relation to other ocean reanalyses
    - L28 Conservation properties of ECCO state estimates
  - L34 How is the ECCO Central Production State Estimate Made?

## `doc/support.rst`
- L1 Getting Help
  - L4 The ECCO Support Mailing List
  - L10 Problems installing Python libraries?

## `ecco_access/Downloading_ECCO_datasets_from_PODAAC/README.md`
- L1 Instructions for Downloading ECCO Datasets hosted on PODAAC
  - L3 Using Python 3 & Jupyter Notebooks
  - L9 Using Command Line *_wget_* and NASA Earthdata

## `ecco_access/Downloading_ECCO_datasets_from_PODAAC/Tutorial_wget_Command_Line_HTTPS_Downloading_ECCO_Datasets_from_PODAAC.md`
- L1 Using _wget_ to Download ECCO Datasets from PO.DAAC
  - L7 Step 1: Create an account with NASA Earthdata
  - L15 Step 2: Set up your ```netrc``` and ```urs_cookies``` files — `_netrc`, `urs_cookies`
  - L47 Step 3: Prepare a list of granules (files) to download
  - L79 Step 4: Download files in a batch with GNU *_wget_*

# @eccov4_gael = `~/.cache/mitgcm-ecco-sources/ecco-docs/ECCOv4` (f2c08e6 2024-06-05 Merge pull request #21 from gaelforget/v1.12c)

## `README.md`
- L1 ECCO Version 4
    - L12 References

## `docs/ECCOv4r1_mods.md`

## `docs/ECCOv4r3_mods.md`
    - L5 install MITgcm + configuration
    - L14 compile the MITgcm config.
    - L25 download input files
    - L35 setup the run directory
    - L47 run the model

## `docs/analyses.rst`
- L4 Analyze
  - L12 Julia Toolbox
  - L17 Python Toolbox
  - L22 Matlab Toolbox
  - L41 Other Resources

## `docs/biblirefs.rst`

## `docs/downloads.rst`
- L4 Download
  - L12 Model Solution
  - L19 Model Setup

## `docs/eccov4r2_dirtree.rst`

## `docs/eccov4r2_output.rst`

## `docs/eccov4r2_setup.rst`

## `docs/index.rst`
- L6 Welcome to ECCO version 4's documentation!

## `docs/introduction.rst`
- L4 Introduction

## `docs/runs.rst`
- L4 Reproduce
  - L19 The Release 2 Solution — `linux_amd64_gfortran`
  - L104 Other Known Solutions — `mitgcmuv_ad`
  - L121 Short Forward Tests
  - L164 Other Short Tests

## `docs/example_scripts/README.md`
- L1 Using cfncluster to run ECCO v4 r2
    - L22 Important warning:
  - L26 Instructions
    - L28 Step 1:
    - L37 Step 2:
    - L42 Step 3:
    - L46 Step 4:
    - L53 Step 5:
    - L60 Step 6:
    - L65 Step 7:
  - L71 References

# @ecco_darwin = `~/Documents/GitHub/ecco_darwin` (be01c37 2026-10-01 Merge branch 'master' of https://github.com/MITgcm-contrib/ecco_darwin)
- `readme.txt` (93 lines) — ECCO-Darwin Github Repository — `checkpoint66o`, `V4r4`, `V4r5`, `llc270`, `llc4320`, `llc1080`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git #`, `git clone https://github.com/MITgcm-contrib/ecco_darwin #`
- `adjoint/README.txt` (43 lines) — Adjoint-based experiments and analyses built on MITgcm/darwin3. — `darwin3.`

## `adjoint/darwin_bling_comparison/README.md`
- L1 darwin_bling_comparison
  - L24 Get the code
  - L37 Relation to Kay's methodology — `run_sensitivity_sweep_fields.sh`, `compute_multipoint_comparison.py`, `ctrl_map_genarr3d.F`, `SpeciesAdj_CO2`
  - L94 How to reproduce — `adxx_ptr1`
  - L124 Bugs found and fixed along the way — `compute_cost_fd.py`, `compute_multipoint_comparison.py`, `DARWIN_CALC_PCO2_APPROX`, `CALC_PCO2_APPROX`
  - L155 Results — `DARWIN_CALC_PCO2_APPROX`, `bling_bio_nitrogen.F`, `global_oce_biogeo_bling_SOCAT`, `PO4_lim`
  - L215 Directory layout — `darwin_vs_bling_multipoint.png`, `darwin_vs_bling_fd_vs_adjoint.png`, `initial_conditions_maps.png`, `initial_conditions_profiles.png`

## `adjoint/darwin_dic_comparison/README.md`
- L1 darwin_dic_comparison — `darwin_bling_comparison`, `multipoint_comparison_results.csv`, `_minus`
  - L16 Retraction: there is no PO4/DOP adjoint corruption — `compute_multipoint_comparison_dic.py`
    - L44 Evidence that the flagged points are ordinary — `run_dic_C`, `run_dic_D`
  - L81 Results (all 38 points, no exclusions)
  - L134 Why the singularity hypothesis was checked and rejected — `DIC_OPTIONS.h`, `DIC_NO_NEG`, `ALLOW_FE`, `DIC_AD_SAFE`
  - L156 Open item — `xx_ptr4`, `wt_DOP.bin`
  - L174 Directory layout

## `adjoint/global_oce_biogeo_bling_SOCAT/README.md`
- L1 global_oce_biogeo_bling_SOCAT — `global_oce_biogeo_bling`
  - L20 Get the code
  - L30 Why this isn't a typical MITgcm verification experiment
  - L43 Directory layout — `global_oce_biogeo_bling`, `dic_init.bin`, `alk_init.bin`, `ones_32b.bin`, `sample_prof.nc`, `mah_flux_smooth.bin`, `socat_pco2_clim_month01..12.nc`, `SOCATv2026_tracks_gridded_monthly.nc`
  - L74 The model configuration — `xx_theta`, `xx_salt`, `xx_ptr1`, `xx_ptr2`, `xx_ptr3`, `xx_ptr4`, `xx_ptr5`, `xx_ptr6`, `socat_pco2_clim_month01..12.nc`, `prof_PCOweight`
  - L110 Building — `mitgcmuv_ad`, `run_optim_loop.sh`
  - L143 Step 1: validate the adjoint (gradient check) — `xx_ptr2`, `xx_ptr1`, `run_optim_loop.sh`
  - L165 Step 2: run the M1QN3 optimization campaign — `mitgcmuv_ad`
    - L186 Continuing past an iteration-cap stop — `m1qn3_output.txt`, `optim_readparms.F`, `run_optim_loop.sh`
  - L219 Performance — `ALLOW_USE_MPI`
  - L246 Known limitations — `xx_genarr3d_weight`, `ones_32b.bin`, `mult_genarr3d`, `SOCATv2026_tracks_gridded_monthly.nc`, `prof_date`, `prof_PCOweight`, `fco2_std_weighted`
  - L310 Reproducing the figures — `BLING_RUN_OPT_DIR`, `BATHY_BIN`, `FIGURES_OUT_DIR`

## `adjoint/global_oce_biogeo_darwin/README.md`
- L1 global_oce_biogeo_darwin — `global_oce_biogeo_bling_SOCAT`
  - L17 Get the code
  - L26 Directory layout
  - L37 The ecosystem configuration: 6+4+0 — `DARWIN_SIZE.h`
  - L54 Reproducing the initial conditions
  - L73 Physics/forcing setup notes — `darwin_check.F`, `EXTERNAL_FIELDS_LOAD`, `tauThetaClimRelax`, `tauSaltClimRelax`, `climsstTauRelax`, `climsssTauRelax`, `CAL_FULLDATE`, `mah_flux_smooth.bin`, `DIAGNOSTICS_SET_POINTERS`
  - L101 Known limitation: pkg/profiles' PCO/PH/CHL/POC support — `pCO2`

## `adjoint/hybrid_darwin_bling/README.md`
- L1 hybrid_darwin_bling
  - L20 Why a surrogate gradient at all
  - L31 Get the code — `mitgcmuv_ad`
  - L45 Mechanism — `mitgcmuv_ad`, `ctrl_map_genarr3d.F`, `wt_DIC.bin`, `wt_ALK.bin`
  - L79 Two configurations — `run_bling`, `run_darwin`, `run_bling_1yr`, `run_darwin_1yr`, `mult_profiles`, `compute_darwin_cost.py`, `patch_ecco_cost.py`, `compute_darwin_cost_1yr.py`, `patch_ecco_cost_1yr.py`, `ctrl_map_ini_genarr.F`, `xx_genarr3d_file`, `darwinPCO2Snap01`, `prof_YYYYMMDD`, `prof_HHMMSS`
  - L105 Running — `pChkptFreq`, `chkptFreq`, `run_hybrid_loop_1yr.sh`, `nTimeSteps`, `m1qn3_output.txt`
  - L131 Results
  - L163 Known limitations / open items

## `adjoint/hybrid_darwin_dic/README.md`
- L1 hybrid_darwin_dic — `pilot_1month`
  - L17 Two attempts
    - L19 1. M1QN3-driven (`scripts/run_hybrid_loop_dic.sh`) -- stalled — `ecco_cost_MIT_CE_000.optNNNN`, `xx_ptr1`
    - L33 2. Manual steepest descent (`scripts/run_hybrid_loop_dic_manual.sh`) -- converges noisily — `adxx_ptr1`, `adxx_ptr2`, `dic_init.bin`, `alk_init.bin`
  - L46 Results — `dic_fc`
  - L78 Slot order in `data.ctrl` — `adxx_salt`, `MDS_READ_FIELD`, `run_dic_A`, `run_dic_B`, `adxx_ptr1`, `adxx_ptr2`, `run_dic`
  - L103 Layout
  - L122 Caveats — `run_dic_A`, `run_dic_B`, `dic_biotic_forcing.F`, `apply_control_to_darwin_ic_dic.py`, `xx_ptr1`, `xx_ptr2`
- `adjoint/tutorial_global_oce_biogeo/readme.txt` (43 lines) — Instructions for building and running adjoint of tutorial_global_oce_biogeo — `linux_amd64_ifort+mpi_ice_nas`, `git clone https://github.com/MITgcm-contrib/ecco_darwin.git git`
- `code_util/LOAC/C_GEM/readme.txt` (138 lines) — Python Idealized Estuary C-GEM model
- `code_util/LOAC/C_GEM/Guayas_estuary/readme.txt` (38 lines) — Config file and boundary conditions for C-GEM config of the Guayas estuary (Ecuador). — `git clone https://github.com/MITgcm-contrib/ecco_darwin.git mkdir`

## `code_util/LOAC/C_GEM/NS_RAD/CLAUDE.md`
- L1 CLAUDE.md
  - L5 What this is
  - L17 Running — `make_report.py`, `build_all.sh`, `make_validation_pdf.py`, `singlechannel_archive`, `netCDF4`, `CGEM_SITE`, `CGEM_MAXT_DAYS`, `CGEM_WARMUP_DAYS`, `CGEM_OUTPUT`, `CGEM_TS`, `CGEM_ICE`, `CGEM_MULTICHANNEL`, `CGEM_N_CHAN_UP`, `CGEM_WATERTEMP_FILE`, `CGEM_DISTANCE`, `CGEM_FONT_SCALE`, `ns_rad`, `tridag_module`, `uphyd_module`, `hyd_module._hyd_iterate`, `disp_sch`, `fun_module`
  - L155 Idealized verification site — `config.IS_IDEALIZED`, `schemes_module.openbound`, `WIND_FILE`, `SOLAR_FILE`, `AIRTEMP_FILE`, `RELHUM_FILE`, `PCO2_FILE`, `SEATEMP_FILE`, `BOUNDARY_FORCING`
  - L185 Multi-site configuration — `CGEM_SITE`, `init_module`, `initialize_substance`
    - L196 Geometry — now observation-based — `B_lb`, `B_ub`, `L_FLARE`, `build_all.sh`, `hyd_module`, `schemes_module.openbound`, `DEPTH_ub`, `n_chan`, `reach_id`, `CGEM_DISTANCE`, `SWORD_MOUTH_SUM`, `n_chan_mod`
    - L305 Boundary conditions (`BOUNDARIES`, per site) — `_baseline.py`, `bin_average`, `DEPTH_ub`
  - L358 Interannual forcing — `build_all.sh`, `freshwater_m3_d`, `DOC_g_d`, `C_DOC`, `biogeo_module.py`, `sed_module.py`, `TOC_cub`, `file_module.py`, `_load`, `repeatYear`, `river_watertemp_2022_degC.csv`, `q_ref`, `heat_module`, `file_module.exfread`, `ns_rad_report.pdf`, `idealized_verification.pdf`, `_paginated_grid`, `ns_rad_validation.pdf`, `VAR_PANELS`, `AREA_VARS`, `ns_rad_diagnostics.pdf`, `rdoc_ox`, `ch4_ox`, `ch4_ex`, `n2o_prod`, `n2o_ex`, `config.ARCTIC_BGC`, `BREAKUP_Q_FACTOR`, `_linregress_np`, `aer_deg` (+4)
    - L621 Full-forcing interannual (2005-2023) — `colville_interannual.py`, `SEATEMP_FILE`, `heat_module.py`, `build_river_temp.py`, `airtemp_2022_degC.csv`, `ice_module.py`, `T_air`, `U_wind`, `I_sw`
    - L706 Sub-daily (diurnal) forcing — `file_module._load`, `file_module.py`, `row_interval_sec`, `WIND_FREQ_SEC`, `AIRTEMP_FREQ_SEC`, `RELHUM_FREQ_SEC`, `SOLAR_FREQ_SEC`, `repeatYear`, `hourly_height`, `_fetch_year`, `qcrad_v3`, `WindSp`, `WindDr`
  - L792 Width and dispersion no longer use the Savenije estuary formulation
    - L796 `WIDTH_MODEL = "flare"` (was: whole-domain exponential) — `L_FLARE`, `B_ub`, `B_lb`
    - L818 `DISPERSION_MODEL = "seo"` (was: Van der Burgh / Savenije) — `DISP_MAX`, `disp_sch`
    - L861 `B` is the TOTAL conveyance width; dispersion uses a separate per-thread width — `CGEM_MULTICHANNEL`, `B_lb`, `B_ub`, `n_chan`, `B_UB_TOTAL`, `N_CHAN_LB`, `L_FLARE`, `B_thread`, `fun_module._piston_velocity_loop`, `aer_deg`, `pCO2`, `_pbar_rho`
  - L964 Forcings — `pCO2_barrow_2022.csv`, `pCO2_Barrow_2022.csv`, `daily_average_weather.csv`
    - L982 Tides — per-river harmonic reconstruction — `fun_module.Tide`
    - L995 Wind-driven storm surge — the real saltwater-intrusion driver (all four rivers) — `fun_module.Tide`, `config.SURGE_FILE`, `surge_prudhoe_2022_m.csv`
    - L1015 Water temperature is a TRANSPORTED field, not a scalar — `water_temp`, `p_bar`, `WATERTEMP_FILE`, `heat_module.py`, `ice_module.py`
    - L1036 Surface heat budget (`heat_module.py`) — `config.HEAT_BUDGET`, `relhum_2022_frac.csv`, `heat_budget`, `rel_hum`, `ICE_MODEL`, `heat_module.ice_energy_deficit`, `ice_module`
    - L1070 Prognostic ice model (`ice_module.py`) — `heat_budget`, `ice_frac`, `variables.ice_thickness`, `config.ICE_MODEL`, `ICE_MODEL`, `heat_module`, `T_FREEZE`, `biogeo_module`, `k_ice_PAR`, `transport_module`, `sed_module`, `file_module.icewrite`, `ice_thickness`
      - L1107 Bottom-fast ice was a one-way ratchet — FIXED — `ice_module`, `ice_frac`, `ICE_FORM_THRESH`
    - L1144 Provenance of the temperature forcing — do not use `watertemp.csv` as a river input — `river_watertemp_2022_degC.csv`, `airtemp_2022_degC.csv`
    - L1183 Discharge
  - L1206 Architecture — `BOUNDARY_FORCING`, `SURGE_FILE`, `Uw_sal`, `Uw_tid`, `fun_module.piston_velocity`, `file_module._FORCING_CACHE`, `init_module`, `lateral_module.py`, `config.ARCTIC_BGC`, `LATERAL_INFLOW`, `ICE_MODEL`, `ice_module`, `ice_frac`, `file_module.exfread`, `nbday_ice`, `P_WATERTEMP`, `P_SEATEMP`, `P_SURGE`, `fun_module`, `_h_solve_kg`, `biogeo_module`, `pCO2`, `_pbar_rho`, `ns_rad_diagnostics.pdf`, `config.OUTPUT_FORMAT`, `file_module.transwrite`
  - L1318 Known defects in the vendored code — `fun_module.pH`, `fun_module._pbar_rho`, `c_pH`, `init_module`, `sed_module`, `tau_dep`, `wMAX`, `kISS`, `sed_module.py`, `tau_dep_lb`, `tau_dep_ub`, `MG_TO_G`, `c_SPM`, `_pbar_rho`, `biogeo_module`, `schemes_module.tvd`, `_tvd`, `RIGID_UPSTREAM_LOWFLOW`, `SEO_CHEONG_FLOOR`, `schemes_module.py`, `fun_module.py`, `transport_module.py`, `init_module.py`, `Chezy_lb`, `Chezy_ub`, `tau_ero`, `piston_velocity`, `Uw_sal`, `Uw_tid`, `Mero_lb` (+7)
  - L1520 Configuration provenance — `config_OLD.py`, `code_python_FCO2`
  - L1535 Where the other documentation lives — `ice_model_plan.md`, `arctic_biogeochemistry.md`, `idealized_verification.md`, `CGEM_CARB_UNITS`, `na_sword_v17b.nc`

## `code_util/LOAC/C_GEM/NS_RAD/README.md`
- L1 NS-RAD — North Slope River-Aquatic-Delta Model
  - L15 What it does
  - L37 Requirements — `netCDF4`
  - L46 Quick start
  - L73 Repository layout — `build_all.sh`
  - L84 Documentation — `_validation.pdf`
  - L95 Status & caveats — `ARCTIC_BGC`, `config.ARCTIC_BGC`
  - L116 Citation
- `code_util/LOAC/C_GEM/NS_RAD/readme.txt` (32 lines) — NS-RAD -- the North Slope River-Aquatic-Delta Model.

## `code_util/LOAC/C_GEM/NS_RAD/docs/FEATURES.md`
- L1 NS-RAD — features added since the original C-GEM code
  - L14 1. Multi-river configuration & run infrastructure — `CGEM_SITE`, `CGEM_MAXT_DAYS`, `CGEM_WARMUP_DAYS`, `CGEM_TS`, `CGEM_OUTPUT`, `CGEM_ICE`, `CGEM_MULTICHANNEL`, `CGEM_N_CHAN_UP`, `CGEM_DISTANCE`, `CGEM_FONT_SCALE`
  - L31 2. Channel geometry & hydrodynamics — `CGEM_MULTICHANNEL`, `n_chan`, `multichannel_test.md`
  - L63 3. Observed 2022 forcings
  - L92 4. Temperature, heat & ice — new physics — `heat_module.py`, `ice_module.py`, `ice_model_plan.md`
  - L108 5. Carbonate chemistry — unit-correct solve + fixes — `CARBONATE_UNITS`
  - L118 6. Arctic biogeochemistry extension *(opt-in)* — `config.ARCTIC_BGC`, `arctic_biogeochemistry.md`, `lateral_module.py`, `LATERAL_INFLOW`, `LATERAL_CONC`
  - L138 7. General time-varying boundary conditions — `BOUNDARY_FORCING`
  - L147 8. Verification & testing — `idealized_verification.md`
  - L157 9. Performance
  - L166 10. Tooling & reproducibility — `fetch_discharge.py`, `build_tides.py`, `build_surge.py`, `build_river_temp.py`, `build_humidity.py`, `lter_boundary.py`, `build_boundary_chem.py`, `build_idealized_forcings.py`, `make_diagnostics_pdf`, `make_validation_pdf`, `make_geometry_pdf`, `make_summary_figures`, `make_river_maps`, `make_schematic_pdf`, `make_idealized_verification_pdf`, `make_movies`, `make_report`, `nsrad_style.py`
  - L179 11. Documentation — `arctic_biogeochemistry.md`, `idealized_verification.md`, `ice_model_plan.md`
  - L187 12. Fixed defects in the vendored code — `fun_module.pH`, `pCO2_barrow`, `pCO2_Barrow`

## `code_util/LOAC/C_GEM/NS_RAD/docs/README.md`
- L1 `docs/` — index
  - L15 Design docs (hand-written) — `model_description.md`, `ice_model_plan.md`, `arctic_biogeochemistry.md`, `config.ARCTIC_BGC`, `idealized_verification.md`, `BOUNDARY_FORCING`, `multichannel_test.md`
  - L30 Generated reports (regenerable) — `ns_rad_report.pdf`, `idealized_verification.pdf`, `ns_rad_interannual.pdf`, `idealized_forcings_preview.png`, `ns_rad_model_summary`, `ns_rad_model_schematic`, `ns_rad_geometry`, `ns_rad_river_networks`, `ns_rad_diagnostics`, `ns_rad_validation`, `make_report.py`
  - L47 Data files (tool inputs / outputs, regenerable) — `performance_metrics.json`, `validation_obs.json`, `usgs_velocity_obs.json`, `sword_widths.json`, `interannual_discharge_obs.json`, `ns_rad_interannual.pdf`
  - L57 Regenerating — `build_all.sh`

## `code_util/LOAC/C_GEM/NS_RAD/docs/arctic_biogeochemistry.md`
- L1 Arctic biogeochemistry extension — `LATERAL_INFLOW`
  - L20 New state variables — `rdoc_ox`, `ch4_ox`, `ch4_ex`, `n2o_prod`, `n2o_ex`
  - L34 1. Refractory DOC + CDOM photomineralisation — `f_photo_lab`, `I_use`
  - L56 2. Methane and nitrous oxide
  - L76 3. Benthic (sediment) efflux
  - L94 4. Distributed lateral loading — `LATERAL_INFLOW`, `variables.q_lat`, `LATERAL_CONC`
  - L110 Verification — `BOUNDARY_FORCING`
  - L127 What is still simplified — `PHOTO_EFF`

## `code_util/LOAC/C_GEM/NS_RAD/docs/ice_model_plan.md`
- L1 River / estuary ice model for NS-RAD — implementation plan — `ice_module.py`, `heat_module.py`, `heat_module`, `k_ice_PAR`, `ICE_MODEL`, `ice_frac`, `ice_thickness`
  - L35 0. The framing problem: the model's "ice" is SEA ice, not river ice
    - L50 Consequence 1: the spring freshet is discarded
    - L66 Consequence 2: every temperature-dependent rate is biased cold
    - L84 Why river ice is a different model from sea ice
    - L98 What this means for the forcing set
  - L105 1. What the model does today — `nbday_ice`, `file_module.exfread`, `transport_module`, `biogeo_module`, `sed_module`, `O2_ex`
  - L141 2. Forcing data
  - L183 3. Tiered implementation
    - L188 Tier 0 — correct the existing defects (no new physics) — `O2_ex`, `biogeo_module`
    - L226 Tier 1 — ice as a diagnostic scalar — `fetch_discharge.py`, `ice_module.py`, `stefan_alpha`, `melt_rate`, `h_open`, `h_full`, `k_ice_PAR`, `albedo_ice`, `fun_module.piston_velocity`, `O2_ex`, `biogeo_module`, `ice_frac`, `ice_h`, `file_module.Rates`
    - L272 Tier 2 — conserve state under ice (the real physics change) — `transport_module`, `biogeo_module`, `sed_module`
    - L295 Tier 3 — hydraulic coupling
    - L309 Tier 4 — optional refinements
  - L317 4. Validation
  - L333 4b. Coupling ice to river temperature (prognostic thermal model)
    - L341 Why the current setup cannot express the feedback — `water_temp`
    - L355 The three changes required — `water_temp`, `fun_module`, `biogeo_module`, `_fhet`, `_o2sat`, `piston_velocity`, `Q_sw`, `Q_sens`, `Q_lw`, `Q_lat`
    - L424 What this would buy
    - L437 Effort and risk, honestly — `fun_module`, `biogeo_module`
  - L453 5. Decisions needed before implementing

## `code_util/LOAC/C_GEM/NS_RAD/docs/idealized_verification.md`
- L1 Idealized verification experiment — `config.IS_IDEALIZED`
  - L13 What makes it "idealized"
  - L28 The feature under test: time-varying boundaries — `schemes_module.openbound`, `BOUNDARY_FORCING`
  - L62 Running it
  - L84 Verifying — `BOUNDARY_FORCING`, `IS_IDEALIZED`
  - L114 Diagnostics PDF and movies — `verify_idealized`, `BOUNDARY_FORCING`
  - L138 Files — `BOUNDARY_FORCING`

## `code_util/LOAC/C_GEM/NS_RAD/docs/model_description.md`
- L1 NS-RAD model description
  - L16 Contents
  - L37 1. Overview and governing equations
  - L63 2. Model domain and grid
  - L76 3. Channel geometry
  - L128 4. Hydrodynamics — `hyd_module`
  - L149 5. Longitudinal dispersion
  - L176 6. Transport scheme
  - L188 7. Water temperature and surface heat budget — `heat_module`
  - L207 8. Prognostic river ice — `ice_module`, `BREAKUP_Q_FACTOR`
  - L234 9. Aquatic carbonate system
  - L264 10. Pelagic biogeochemistry
    - L272 10.1 Primary production
    - L299 10.2 Organic-matter degradation and nitrogen cycling
    - L312 10.3 Gas exchange (piston velocity)
    - L332 10.4 Reaction stoichiometry
  - L354 11. Suspended sediment
  - L362 12. Arctic biogeochemistry extension — `config.ARCTIC_BGC`, `arctic_biogeochemistry.md`
    - L370 12.1 Refractory DOC and CDOM photomineralisation
    - L386 12.2 Methane and nitrous oxide
    - L405 12.3 Benthic (sediment) exchange
    - L420 12.4 Distributed lateral loading
  - L434 13. Boundary conditions and forcing — `BOUNDARY_FORCING`
  - L449 14. Numerical implementation
  - L465 15. Parameter tables
    - L467 15.1 Domain and numerics
    - L478 15.2 Hydrodynamics and sediment
    - L489 15.3 Phytoplankton and light
    - L504 15.4 Nutrient/oxygen half-saturations
    - L518 15.5 Reaction rates and Redfield ratios
    - L529 15.6 Arctic biogeochemistry extension
  - L549 16. References

## `code_util/LOAC/C_GEM/NS_RAD/docs/multichannel_test.md`
- L1 Multi-channel geometry — `__pycache__`, `colville_minor.py`
  - L19 The problem — `fun_module.river_dispersion`, `B_lb`, `B_ub`
  - L37 The two changes
    - L39 1. Width-role separation (`CGEM_MULTICHANNEL=on`) — `B_UB_TOTAL`, `n_chan`, `N_CHAN_LB`, `L_FLARE`, `B_thread`, `B_ub`
    - L57 2. Distributary resolution (`sites/colville_main.py`, `colville_minor.py`) — `Q_FRACTION`
  - L81 The data: read the braided total, don't reconstruct it — `B_ub`, `n_chan_mod`, `B_UB_TOTAL`
  - L106 Results
    - L120 Sensitivity to the thread count — the load-bearing result — `S_max`
    - L143 Adoption: all four rivers (`tools/compare_adoption.py`) — `fun_module._piston_velocity_loop`, `aer_deg`, `n_chan`
    - L175 What this does and does not do for the salinity misfit
  - L195 Caveats — `DISP_MAX`
  - L208 Verifying a run actually used the geometry you think — `CGEM_MULTICHANNEL`

## `code_util/LOAC/C_GEM/NS_RAD/docs/performance.md`
- L1 Speeding up NS-RAD — what was done and how it was verified — `CGEM_MAXT_DAYS`, `CGEM_WARMUP_DAYS`
  - L20 1. Measure first — `cProfile`, `tridag_module.coeff_a`, `tridag_module.tridag`, `uphyd_module.update`, `tridag_module.conv`, `density._seck`, `biogeo_module.biogeo`, `density._dens0`, `fun_module.K1_CO2`, `fun_module.K2_CO2`, `fun_module.p_bar`, `K1_CO2`, `K2_CO2`
    - L59 An aside worth recording
  - L72 2. What was changed
    - L74 2a. Hydrodynamic kernels → numba (`tridag_module.py`, `uphyd_module.py`) — `coeff_a`, `new_uh`
    - L102 2b. Density stack → numba (`density.py`) — `_dens0`, `_seck`, `fun_module.p_bar`
    - L115 2c. Hoist the carbonate constants (`biogeo_module.py`)
  - L127 3. Verification — bit-identity, not "looks right" — `__pycache__`
  - L176 4. What was deliberately NOT done
  - L214 5. What is left
  - L220 6. Environment — `fun_module.pH`
  - L228 7. Third pass — the Python orchestration layer (Tier 1 + Tier 2) — `disp_sch`, `cProfile`, `sed_module.sed`, `fun.river_dispersion`, `fun.piston_velocity`, `transport_module`, `sed_module`, `wISS`, `fun_module`, `piston_velocity`, `river_dispersion`, `_piston_velocity_loop`, `_river_dispersion_loop`, `_d_o2`, `hyd_module`, `BOUNDARY_FORCING`
  - L279 8. What is left now — `disp_sch`, `BOUNDARY_FORCING`, `file_module.write`, `file_module`
- `code_util/LOAC/C_GEM/North Slope/readme.txt` (30 lines) — Python Idealized Estuary C-GEM model for the North Slope (Alaska) with
- `code_util/LOAC/C_GEM/code_python_v2/readme.txt` (130 lines) — Python Idealized Estuary C-GEM model v2
- `code_util/LOAC/GloFas/README.txt` (104 lines) — Create daily freshwater runoff forcing files for ECCO from Global Flood Awareness System (GloFas), Copernicus, 3 min or 0.05x0.05 degree, — `git clone --depth 1`
- `code_util/LOAC/GlobalNews/readme.txt` (38 lines) — GlobalNEWS biogeochemical river exports (DIN, DON, DIP, DOP, DOC, PN, PP, POC, DSi)
- `code_util/bathy/readme.txt` (23 lines) — This folder contains MATLAB code written by Dustin Carroll for accurately re-gridding bathymetric products onto various LLC and lon-lat grids
- `code_util/transport/readme.txt` (1 lines) — ECCO-Darwin transport
- `idealized/1D_darwin/sea_ice_column/readme.txt` (27 lines) — v06 1D_ocean_ice_column_darwin idealized experiment — `darwin3-darwin`
- `offline/V4r6_darwin_offline/readme.txt` (394 lines) — Offline ECCO-Darwin: LLC90 biogeochemistry driven by pre-computed ECCOv4r6 physics — `checkpoint68g`, `darwin3`, `darwin3/pkg/darwin`, `darwin3/tools/genmake2`, `darwin3/tools/darwin/cogapp/cogapp.py`, `linux_amd64_ifort+mpi_ice_nas.electra_skylake_20251013`, `linux_amd64_ifort+mpi_ice_nas`, `V4r6`, `LLC90`, `LLC270`, `llc270`, `git clone https://github.com/MITgcm/MITgcm.git -b`, `git clone https://github.com/MITgcm-contrib/ecco_darwin.git cd`, `git clone https://github.com/ECCO-GROUP/ECCO-v4-Configurations.git mv`

## `pkg/wad/README.md`
- L1 pkg/wad — wetting and drying for MITgcm — `addMass`
  - L15 Contents
  - L32 Install — `WAD_PARM01`
  - L61 Required model settings — `wad_check`, `vectorInvariantMomentum`, `upwindShear`
  - L89 Parameters (`data.wad`, namelist `WAD_PARM01`) — `wadMinDepth`, `wadCritDepth`, `wadAdvDepth`, `wadGMDepth`, `wadKPPCapDepth`, `wadKPPDepth`, `wadDryForcing`, `wadUpwindFace`, `wadCarryVel`, `momImplVertAdv`, `wadDragDepth`, `wadManningN`, `ALLOW_BOTTOMDRAG_ROUGHNESS`, `wadDragMax`, `wadManningFile`, `wadMaxFroude`, `wadMaxSpeed`, `wadConserveVol`, `wadMonFreq`, `wadSmoothWidth`
  - L112 Diagnostics and monitor — `WADdryC`, `pickup_wad`, `pickup_wadface`
  - L121 Verification experiments — `wad_thacker_1d`, `wad_balzano`, `wad_flat_xz`, `wad_estuary_3d`, `wad_mudflat`
  - L138 Cost and time step
  - L145 Limitations — `wadSmoothWidth`
  - L156 References

## `pkg/wad/doc/wad_formulation.md`
- L1 pkg/wad: formulation
  - L6 1. What stock MITgcm does when a column empties — `CALC_R_STAR`, `hFacInf`, `Rmin_surf`, `hFacW`, `maskW`
  - L27 2. Definitions — `wadCritDepth`
  - L36 3. The face rule — `wadCritDepth`
  - L57 4. Closing a face in the geometry — `UPDATE_R_STAR`, `hFacW`
  - L70 5. Where pkg/wad sits in the time step
  - L88 6. Thickness of an open face — `wadUpwindFace`
  - L107 7. Outflow (positivity) limiter — `useRealFreshWaterFlux`
  - L136 8. Film floor and volume budget — `wadCritDepth`, `wadMonFreq`, `addMass`
  - L146 9. Momentum in thin water — `upwindShear`, `wadCritDepth`, `wadAdvDepth`, `wadCarryVel`
  - L162 10. Bottom drag
  - L177 11. Tapers on shallow columns — `wadCritDepth`, `wadKPPCapDepth`, `wadGMDepth`
  - L194 12. Verification
  - L207 13. Ripples, time step and cost

## `pkg/wad/verification/wad_balzano/README.md`
- L1 wad_balzano: Balzano (1998) tidal-flat tests: slope, step, pool — `wadCritDepth`
  - L9 Set-up
  - L18 Variants
  - L26 Inputs
  - L31 Build and run
  - L62 What to check — `wadMonFreq`, `minDepth`, `wadMinDepth`

## `pkg/wad/verification/wad_estuary_3d/README.md`
- L1 wad_estuary_3d: 3-D tidal estuary with drying flats
  - L9 Set-up
  - L22 Variants
  - L30 Inputs — `WAD_NX`, `WAD_NY`, `WAD_NR`
  - L36 Build and run
  - L67 What to check — `wadMonFreq`, `minDepth`, `wadMinDepth`

## `pkg/wad/verification/wad_flat_xz/README.md`
- L1 wad_flat_xz: Stratified beach at rest, and under a tide
  - L9 Set-up
  - L19 Variants
  - L25 Inputs
  - L29 Build and run
  - L60 What to check — `wadMonFreq`, `minDepth`, `wadMinDepth`

## `pkg/wad/verification/wad_mudflat/README.md`
- L1 wad_mudflat: Macrotidal mudflat with tidal creeks and a river
  - L8 Set-up — `addMass`
  - L19 Variants — `useRealFreshWaterFlux`
  - L26 Inputs — `addMass.bin`
  - L31 Build and run
  - L62 What to check — `wadMonFreq`, `minDepth`, `wadMinDepth`

## `pkg/wad/verification/wad_thacker_1d/README.md`
- L1 wad_thacker_1d: Thacker (1981) oscillating basin
  - L9 Set-up
  - L18 Variants
  - L25 Inputs — `bathy_thacker.bin`, `eta_thacker.bin`
  - L30 Build and run
  - L61 What to check — `wadMonFreq`, `minDepth`, `wadMinDepth`
- `regions/BayBengal/readme.txt` (51 lines) — Regional Bay of Bengal ECCO-Darwin cutout — `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `git clone https://github.com/darwinproject/darwin3 git`
- `regions/CCS/matlab/extract_cutout/readme.txt` (62 lines) — to fix volume drift, need to convert UVELMASS to UVEL
- `regions/CCS/v05/readme_macos.txt` (40 lines) — Instructions for building and running CCS regional simulation — `darwin3`, `darwin3/run`, `darwin_arm64_gfortran`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/CCS/v05/readme_pleiades.txt` (33 lines) — Instructions for building and running CCS regional simulation — `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/CCS/v05_kelp/readme_pleiades_kelp.txt` (66 lines) — Instructions for building and running CCS regional simulation — `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `llc270`, `git clone --depth 1`, `git clone https://github.com/darwinproject/darwin3 cd`, `git checkout 24885b71 mkdir`
- `regions/GoA/readme.txt` (113 lines) — https://github.com/MITgcm-contrib/ecco_darwin/tree/master/regions/GoA — `darwin3.readthedocs.io/en/latest/phys_pkgs/darwin.html`
- `regions/GoA/v05/readme_macos.txt` (39 lines) — Instructions for building and running Gulf of Alaska (GoA) regional simulation — `darwin3`, `darwin3/run`, `darwin_arm64_gfortran`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/GoA/v05/readme_pleiades.txt` (33 lines) — Instructions for building and running Gulf of Alaska (GoA) regional simulation — `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`

## `regions/GoM/GoM_1km/README.md`
- L1 Gulf of Mexico downscaled regional set-up
  - L4 Preliminary information
  - L9 1. Get code
  - L21 2. Build executable
  - L33 3. Instructions for running simulation
- `regions/GoM/high_res/readme.txt` (48 lines) — Guld of Mexico regional setup based on LLC4320 grid + LLC270 IC/BC — `linux_amd64_ifort+mpi_ice_nas`, `LLC4320`, `LLC270`, `git clone https://github.com/MITgcm/MITgcm.git git`
- `regions/GoM/high_res/readme_darwin.txt` (66 lines) — Guld of Mexico regional setup based on LLC4320 grid + LLC270/ECCO-Darwin IC/BC — `darwin3`, `darwin3/configurations/downscaled_greenland/L1/L1_GOM/exf`, `darwin3/configurations/downscaled_greenland/L1/L1_GOM/obcs`, `linux_amd64_ifort+mpi_ice_nas`, `LLC4320`, `LLC270`, `git clone --depth 1`, `git clone https://github.com/darwinproject/darwin3 cd`, `git checkout 24885b71 mkdir`
- `regions/GoM/llc270/readme.txt` (113 lines) — https://github.com/MITgcm-contrib/ecco_darwin/tree/master/regions/GoA — `darwin3.readthedocs.io/en/latest/phys_pkgs/darwin.html`
- `regions/GoM/llc270/v05/readme_macos.txt` (38 lines) — Instructions for building and running Gulf of Alaska (GoA) regional simulation — `darwin3`, `darwin3/run`, `darwin_arm64_gfortran`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/GoM/llc270/v05/readme_pleiades.txt` (33 lines) — Instructions for building and running Gulf of Alaska (GoA) regional simulation — `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/GoM/llc90/readme_macos.txt` (63 lines) — Instructions for building and running Gulf of Mexico (GoM) — `darwin3`, `llc90`, `git clone --depth 1`, `git clone https://github.com/darwinproject/darwin3 #`, `git checkout backport_ckpt68y mkdir`
- `regions/GoM/llc90/readme_pleiades.txt` (34 lines) — Instructions for building and running Gulf of Mexico (GoM) — `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout darwin_ckpt68y ==============`
- `regions/GoM/llc90/readme_ubuntu.txt` (41 lines) — Instructions for building and running GoML regional simulation — `darwin3`, `darwin3/run`, `darwin_arm64_gfortran`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/GoM/llc90/readme_windows.txt` (73 lines) — Instructions for building and running Gulf of Mexico (GoM) — `darwin3`, `darwin3/run`, `darwin_arm64_gfortran`, `llc90`, `git clone --depth 1`, `git clone https://github.com/darwinproject/darwin3 #`, `git checkout backport_ckpt68y mkdir`
- `regions/GulfGuinea/v05/readme_macos.txt` (39 lines) — Instructions for building and running Gulf of Guinea regional — `darwin3`, `darwin3/run`, `darwin_arm64_gfortran`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/LR17/80m/readme_macos.txt` (34 lines) — Instructions for building and running LR17 regional simulation — `git clone git@github.com:MITgcm/MITgcm.git git`
- `regions/LR17/v05/readme_macos.txt` (45 lines) — Instructions for building and running LR17 regional simulation — `darwin3`, `darwin_arm64_gfortran`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 cp`
- `regions/LR17/v05/readme_pleiades.txt` (33 lines) — Instructions for building and running LR17 regional simulation — `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/LR17/v05/readme_ubuntu.txt` (40 lines) — Instructions for building and running LR17 regional simulation — `darwin3`, `darwin3/run`, `darwin_arm64_gfortran`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/Med/readme.txt` (101 lines) — https://github.com/MITgcm-contrib/ecco_darwin/tree/master/regions/Med — `darwin3.readthedocs.io/en/latest/phys_pkgs/darwin.html`
- `regions/Med/v05/readme_macos.txt` (39 lines) — Instructions for building and running Mediterranean regional simulation — `darwin3`, `darwin3/run`, `darwin_arm64_gfortran`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/Med/v05/readme_pleiades.txt` (33 lines) — Instructions for building and running Mediterranean regional simulation — `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/Med/v05/readme_ubuntu.txt` (45 lines) — Instructions for building and running Mediterranean regional simulation — `darwin3`, `darwin3/run`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/Med/v05/readme_zorbas.txt` (65 lines) — Instructions for building and running Mediterranean regional simulation — `darwin3`, `darwin3/verification`, `darwin3/run`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/RedSea/readme.txt` (101 lines) — https://github.com/MITgcm-contrib/ecco_darwin/tree/master/regions/RedSea — `darwin3.readthedocs.io/en/latest/phys_pkgs/darwin.html`
- `regions/RedSea/Coral/readme.txt` (7 lines) — Coral bleaching model developed by Ioannis Chatzonikolakis
- `regions/RedSea/kaust_baby/readme_macos.txt` (62 lines) — Building RedSea/kaust_baby on arm64 macOS — `darwin3`, `darwin3/build`, `darwin3/run`, `darwin_arm64_gfortran`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/RedSea/kaust_v1s3/readme_darwin.txt` (31 lines) — Instructions for building and running a RedSea/kaust_v1s3 — `darwin3`, `darwin3/build`, `darwin3/run`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/RedSea/v05/readme_macos.txt` (43 lines) — Instructions for building and running a Red Sea regional simulation — `darwin3`, `darwin3/run`, `darwin_arm64_gfortran`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/RedSea/v05/readme_pleiades.txt` (33 lines) — Instructions for building and running a Red Sea regional simulation — `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/RedSea/v05/readme_ubuntu.txt` (39 lines) — Instructions for building and running Red Sea regional simulation — `darwin3`, `darwin3/run`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/RedSea/v05/readme_ubuntu_4cpu.txt` (45 lines) — Instructions for building and running Red Sea regional simulation with 4 proccesors — `darwin3`, `darwin3/run`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/RedSea/v05_coral/readme_macos.txt` (39 lines) — Instructions for building and running a Red Sea regional simulation — `darwin3`, `darwin3/run`, `darwin_arm64_gfortran`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/RedSea/v05_coral/readme_pleiades.txt` (33 lines) — Instructions for building and running a Red Sea regional simulation — `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/RedSea/v05_coral/readme_ubuntu.txt` (39 lines) — Instructions for building and running Red Sea regional simulation — `darwin3`, `darwin3/run`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/RedSea/v05_kpp/readme.txt` (32 lines) — This code replaces the single-column swfrac with swfrac2d to add
- `regions/RedSea/v05_kpp/readme_macos.txt` (55 lines) — Instructions for building and running a Red Sea regional simulation — `darwin3`, `darwin3/run`, `darwin_arm64_gfortran`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/RedSea/v05_kpp/readme_ubuntu.txt` (55 lines) — Instructions for building and running Red Sea regional simulation — `darwin3`, `darwin3/run`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/RedSea/v05_kpp/verification/readme.txt` (35 lines) — A set of 1-day ("endtime=86400," in data) verification experiments were — `git checkout 7d24893", that`
- `regions/RedSea/v05_swfrac/readme.txt` (13 lines) — This code replaces the single-column swfrac with swfrac2d to add
- `regions/RedSea/v05_swfrac/readme_macos.txt` (52 lines) — Instructions for building and running a Red Sea regional simulation — `darwin3`, `darwin3/run`, `darwin_arm64_gfortran`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/RedSea/v05_swfrac/readme_ubuntu.txt` (56 lines) — Instructions for building and running Red Sea regional simulation — `darwin3`, `darwin3/run`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/RedSea/v05_swfrac/verification/readme.txt` (37 lines) — A set of 1-day ("endtime=86400," in data) verification experiments were
- `regions/bering_strait/readme.txt` (113 lines) — https://github.com/MITgcm-contrib/ecco_darwin/tree/master/regions/GoA — `darwin3.readthedocs.io/en/latest/phys_pkgs/darwin.html`
- `regions/bering_strait/v05/readme_macos.txt` (39 lines) — Instructions for building and running Gulf of Alaska (GoA) regional simulation — `darwin3`, `darwin3/run`, `darwin_arm64_gfortran`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`
- `regions/bering_strait/v05/readme_pleiades.txt` (33 lines) — Instructions for building and running Gulf of Alaska (GoA) regional simulation — `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `llc270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`, `git checkout 24885b71 mkdir`

## `regions/downscaling/README.md`
- L1 How to generate an ECCO regional cut out
  - L6 General information
  - L12 Main steps — `diagnostic_vec`
  - L20 Getting Started
    - L26 1. Get ECCO-Darwin v5 set-up & merge ``diagnostic_vec`` with MITgcm — `diagnostic_vec`
    - L64 2. Create the python3 anaconda environment

## `regions/downscaling/STEP1.md`
- L1 Generate regional cut-out input files
  - L4 Preliminary information — `diagnostics_vec`
  - L12 I. Generate the ``.mitgrid`` file — `gen_mitgrid.py`
  - L43 II. Generate the bathymetry file — `gen_bathy.py`
  - L81 III. Generate tiles for multiprocessing
  - L105 IV. Generate grid file (netcdf)
    - L110 a. Set the cut-out configuration parameters
    - L142 b. Compile and run the regional model on 1 step — `regnm_bathymetry.bin`
    - L177 c. Stitch the grid tiles in a netcdf file — `stitch_ncgrid.py`
  - L210 V. Generate mask files for ``diagnostics_vec``
    - L212 a. Compile and run ECCO global state estimate on 1 step — `job_ECCO_darwin`
    - L258 c. Generate the masks files — `gen_dvmasks.py`

## `regions/downscaling/STEP2.md`
- L1 Extract regional model boundary information
  - L4 Preliminary information — `diagnostic_vec`
  - L13 I. Prepare the simulation for ``diagnostic_vec`` extraction
    - L15 a. Turn on ``diagnostic_vec`` package
    - L32 b. Generate a ``data.diagnostics_vec`` parameter file — `data.diagnostics_vec`, `diagnostic_vec`
    - L43 c. Set the compile time ``DIAGNOSTICS_VEC_SIZE.h`` file — `DIAGNOSTICS_VEC_SIZE.h`, `gen_dvmasks`, `data.diagnostics_vec`
  - L55 II. Compile and run the simulation — `data.diagnostics_vec`

## `regions/downscaling/STEP3.md`
- L1 Generate regional set-up, initial and boundary conditions
  - L5 Preliminary information
  - L15 I. Generate the initial conditions — `gen_pickups.py`
  - L64 II. Generate the boundary conditions — `gen_obcs.py`
  - L94 III. Verify the boundary conditions (recommended) — `gen_obcs.py`
    - L109 III.1 Boundary velocity sections — `plot_boundary_velocity_profile.py` — `gen_obcs.py`, `gen_obcs2.py`
    - L144 III.2 Boundary transport check — `check_obcs_transport.py` — `gen_obcs.py`, `diagnostics_vec`, `gen_obcs2.py`, `_ncgrid.nc`, `HFacS`, `gen_mod`

## `regions/downscaling/STEP4.md`
- L1 Build and run your downscaled regional set-up
  - L4 Preliminary information
  - L9 I. Prepare code directories
  - L50 II. Prepare namelists
    - L59 a. data namelist (ecco_darwin/regions/YOURSETUP/code/data)
    - L108 b. data.cal namelist (ecco_darwin/regions/YOURSETUP/inputs/data.cal)
    - L113 c. data.exf namelist (ecco_darwin/regions/YOURSETUP/inputs/data.exf)
    - L133 d. data.obcs namelist (ecco_darwin/regions/YOURSETUP/inputs/data.obcs)
    - L137 e. data.ggl90 namelist (ecco_darwin/regions/YOURSETUP/inputs/data.ggl90)
    - L145 f. data.diagnostics namelist (ecco_darwin/regions/YOURSETUP/inputs/data.diagnostics)
    - L149 g. data.pkg namelist (ecco_darwin/regions/YOURSETUP/inputs/data.pkg)
    - L165 h. Additional namelists for Darwin run (ecco_darwin/regions/YOURSETUP/inputs_darwin/)
  - L190 III. Prepare forcing files
    - L192 a. Forcing files for physics
    - L198 b. Additional forcings if running with Darwin
  - L222 IV. Compile and run
    - L226 a. Compile
    - L250 b. Run

## `regions/kerguelen/README.md`
- L1 Kerguelen Plateau ~1 km regional ECCO-Darwin configuration
  - L20 Contents
  - L34 What is **not** here — `kerguelen_ncgrid.nc`, `kerguelen_bathymetry.bin`, `diagnostics_vec`
  - L50 Two things that will silently misconfigure a rebuild — `code_darwin`, `OBCS_OPTIONS.h`, `DIAGNOSTICS_SIZE.h`, `ALLOW_OBCS_SPONGE`
  - L65 Provenance

## `regions/kerguelen/v05/PIPELINE.md`
- L1 How the Kerguelen inputs were generated
  - L14 STEP 1 — regional grid and bathymetry — `run_step1_mitgrid.sh`, `job_step1_bathy.sh`, `kerguelen_bathymetry.bin`, `job_step1_stitch.sh`, `kerguelen_ncgrid.nc`, `job_step1_dvmasks.sh`, `delX`, `delY`, `delYFile`
  - L37 STEP 2 — parent-side extraction — `nVEC_mask`, `nml_vecFiles`, `VEC_points`
  - L61 STEP 3 — regional boundary conditions and initial conditions
    - L63 OBCS files — `job_gen_obcs_fast.sh` -> `<config>/forcings/OBCS/` (164 files, ~100 GB) — `gen_obcs_fast.py`, `THETA_east.bin`
    - L70 Barotropic transport correction — `job_gen_obcs_transport_correction.sh` — `nonlinFreeSurf`, `exactConserv`
    - L91 Initial conditions — `job_gen_pickups.sh` -> `<config>/forcings/pickups/` (78 files, ~50 GB) — `gen_pickups.py`, `hydrogThetaFile`, `hydrogSaltFile`
  - L104 STEP 4 — regional run — `readme_pleiades.txt`, `code_darwin`
  - L112 Regeneration order
- `regions/kerguelen/v05/readme_pleiades.txt` (134 lines) — Instructions for building and running the Kerguelen Plateau ~1 km regional — `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `LLC270`, `git clone git@github.com:MITgcm-contrib/ecco_darwin.git git`

## `regions/kerguelen/v05/scripts/pipeline/README.md`
- L1 Pipeline scripts (STEP1–STEP3) — `run_step1_mitgrid.sh`, `run_step1_mitgrid_repad.sh`, `job_step1_bathy.sh`, `kerguelen_bathymetry.bin`, `job_step1_stitch.sh`, `kerguelen_ncgrid.nc`, `job_step1_dvmasks.sh`, `job_step2_parent_dv.sh`, `diagnostics_vec`, `job_gen_obcs_fast.sh`, `job_gen_obcs_transport_correction.sh`, `job_gen_pickups.sh`
- `regions/mac_delta/LatLon/readme.txt` (56 lines) — Mackenzie Delta regional setup based on LatLon — `linux_amd64_ifort+mpi_ice_nas`, `git clone https://github.com/MITgcm/MITgcm.git git`
- `regions/mac_delta/LatLon/readme_darwin.txt` (66 lines) — Mackenzie Delta regional setup based on LatLon w/ Darwin — `darwin3`, `darwin3/configurations/downscaled_greenland/L1/L1_mac_delta/exf`, `darwin3/configurations/downscaled_greenland/L1/L1_mac_delta/obcs`, `linux_amd64_ifort+mpi_ice_nas`, `git clone --depth 1`, `git clone https://github.com/darwinproject/darwin3 cd`, `git checkout 24885b71 mkdir`

## `regions/mac_delta/llc270/biogeochem_setup/README.md`
- L1 ECCO-Darwin: Mackenzie Delta configuration
  - L7 1. Get the code
  - L13 2. Build executable
  - L27 3. Run the setup
    - L33 Make links to forcing files
    - L50 Copy setup data files
    - L55 Run

## `regions/mac_delta/llc270/biogeochem_setup/carroll_2020_ecosystem/CDOM_setup/README.md`
- L1 Baseline Setup Content
  - L3 Introduction
  - L13 Ecosystem Characteristics
  - L26 Tracers (n=33)

## `regions/mac_delta/llc270/biogeochem_setup/carroll_2020_ecosystem/ED-SBS_Bertin_etal_2023/README.md`
- L1 ED-SBS: ECCO-Darwin Mackenzie Delta configuration
  - L7 1. Get the code
  - L12 1. Build executable
  - L27 2. Download the forcing files
  - L35 3. Prepare the simulation
  - L61 4. Run the code

## `regions/mac_delta/llc270/biogeochem_setup/carroll_2020_ecosystem/ED-SBS_Bertin_etal_2023/setup_files/README.md`
- L1 Baseline Setup Content
  - L3 Introduciton
  - L11 Ecosystem Characteristics
  - L24 Tracers (n=32)

## `regions/mac_delta/llc270/biogeochem_setup/carroll_2020_ecosystem/baseline_setup/README.md`
- L1 Baseline Setup Content
  - L3 Introduciton
  - L10 Ecosystem Characteristics
  - L23 Tracers (n=31)
- `regions/mac_delta/llc270/biogeochem_setup/oldAO_dev/readme.txt` (69 lines) — The new Darwin ecosystem has been created to simulate general plankton dynamics taking place in the Arctic Ocean. — `darwin3`, `git clone -b cdom-carbon`
- `regions/mac_delta/llc270/biogeochem_setup/oldAO_dev/Setup_files_v0/code_darwin/readme.txt` (39 lines) — From Steph (4/20/2021):
- `regions/mac_delta/llc270/biogeochem_setup/oldAO_dev/Setup_files_v0/input/readme.txt` (18 lines) — From Steph (4/20/2021)
- `regions/mac_delta/llc270/physics_setup/readme.txt` (70 lines) — Mackenzie Delta regional setup based on LLC270 — `linux_amd64_ifort+mpi_ice_nas`, `PATH_TO_MPI_ENVIRONMENT_VARIABLE`, `LLC270`, `llc270`, `git clone https://github.com/MITgcm/MITgcm.git #svn`, `git clone https://github.com/MITgcm-contrib/ecco_darwin ln`
- `regions/mac_delta/llc4320/readme.txt` (85 lines) — Mackenzie Delta regional setup based on LLC4320 — `linux_amd64_ifort+mpi_ice_nas`, `linux_amd64_ifort+gcc`, `LLC4320`, `llc4320`, `git clone https://github.com/MITgcm/MITgcm.git svn`
- `regions/mac_delta/llc4320/AO_ecosystem/readme.txt` (46 lines) — ARCTIC ECOSYSTEM README FILE:
- `regions/mac_delta/llc4320/AO_ecosystem/code_darwin/readme.txt` (39 lines) — From Steph (4/20/2021):
- `regions/mac_delta/llc4320/AO_ecosystem/input/readme.txt` (18 lines) — From Steph (4/20/2021)
- `regions/north_slope/readme.txt` (113 lines) — https://github.com/MITgcm-contrib/ecco_darwin/tree/master/regions/GoA — `darwin3.readthedocs.io/en/latest/phys_pkgs/darwin.html`
- `regions/totten/readme.txt` (56 lines) — %Contents will no longer updated here. Please check — `darwin3-dev`, `Darwin3`, `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `git clone git://gud.mit.edu/darwin3-dev darwin3`, `git checkout ecco_darwin_v4_llc270_darwin3 cd`
- `v02/cs510_Brix/readme.txt` (137 lines) — Build executable for ECCO-Darwin version 2 (ag4) — `checkpoint62w`, `linux_amd64_ifort+mpi_ice_nas`, `cs510`, `CS510`, `git clone --depth 1`
- `v02/cs510_Manizza_AO/readme.txt` (11 lines) — ECCO2-Darwin AO simulation used in Manizza et al. publications — `cs510`
- `v02/cs510_Manizza_AO/code/readme.txt` (4 lines) — originally obtained from ~hbrix/MITgcm/code_darwin_p6 — `linux_amd64_ifort+mpi_ice_nas`
- `v03/cs510_Brix/readme.txt` (55 lines) — Build executable for ECCO-Darwin version 3 — `linux_amd64_ifort+mpi_ice_nas`, `cs510`, `git clone --depth 1`
- `v03/cs510_latest/readme.txt` (46 lines) — Build executable for ECCO-Darwin version 3 — `linux_amd64_ifort+mpi_ice_nas`, `cs510`, `git clone --depth 1`
- `v04/3deg/readme budget.txt` (63 lines) — Verification experiment, initially based on — `h_mpi`, `data_mpi`, `git clone --depth 1`
- `v04/3deg/readme.txt` (53 lines) — Verification experiment, initially based on — `h_mpi`, `data_mpi`, `git clone --depth 1`
- `v04/llc270_JAMES_budget/readme_darwin.txt` (34 lines) — Instructions for llc270 ECCO-Darwin simulation — `linux_amd64_ifort+mpi_ice_nas`, `llc270`, `git clone --depth 1`
- `v04/llc270_JAMES_paper/readme/readme_85_92.txt` (66 lines) — Instructions for llc270 physical simulation without Darwin, — `linux_amd64_ifort+mpi_ice_nas`, `llc270`, `git clone https://github.com/MITgcm-contrib/ecco_darwin cd`
- `v04/llc270_JAMES_paper/readme/readme_ecco_darwin.txt` (42 lines) — Instructions for llc270 ECCO-Darwin simulation — `linux_amd64_ifort+mpi_ice_nas`, `llc270`, `git clone --depth 1`
- `v04/llc270_JAMES_paper/readme/readme_ecco_darwin_verification.txt` (39 lines) — Instructions for llc270 verification experiment — `linux_amd64_ifort+mpi_ice_nas`, `llc270`, `git clone --depth 1`
- `v04/llc270_JAMES_paper/readme/readme_physics.txt` (63 lines) — Instructions for llc270 physical simulation without Darwin — `linux_amd64_ifort+mpi_ice_nas`, `llc270`, `git clone https://github.com/MITgcm-contrib/ecco_darwin cd`
- `v04/llc270_OAE_ship_track_paper/readme_ecco_darwin.txt` (49 lines) — Instructions for llc270 ECCO-Darwin simulation — `linux_amd64_ifort+mpi_ice_nas`, `llc270`, `git clone --depth 1`
- `v04/llc270_devel/readme_darwin.txt` (39 lines) — development version of ECCO-Darwin with nonlinear water-column dissolution and — `linux_amd64_ifort+mpi_ice_nas`, `git clone --depth 1`
- `v05/1deg/readme_darwin_v4r4.txt` (51 lines) — v05 1deg Darwin3 simulation based on ECCOV4r4 set-up — `Darwin3`, `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `V4r4`, `git clone --branch darwin_ckpt68d_at_c66g`, `git clone --depth 1`
- `v05/1deg/readme_darwin_v4r5.txt` (54 lines) — v05 1deg Darwin3 simulation based on ECCOV4r5 set-up — `Darwin3`, `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `v4r4`, `V4r5`, `LLC90`, `git clone --branch backport_ckpt68g`, `git clone --depth 1`
- `v05/1deg/readme_v4r4.txt` (38 lines) — v05 1deg Darwin3 simulation based on ECCOV4r4 set-up — `Darwin3`, `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `git clone --branch backport-c66g`, `git clone --depth 1`
- `v05/1deg/readme_v4r5.txt` (46 lines) — v05 1deg Darwin3 simulation based on ECCOV4r5 set-up — `checkpoint68g`, `Darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `v4r4`, `LLC90`, `git clone https://github.com/MITgcm/MITgcm.git -b`, `git clone --depth 1`
- `v05/1deg_CDR/readme_darwin.txt` (58 lines) — v05 1deg Darwin3 simulation based on ECCOV4r4 set-up — `Darwin3`, `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `V4r4`, `git clone --branch darwin_ckpt68d_at_c66g`, `git clone --depth 1`
- `v05/1deg_RADIv2/readme_darwin_v4r5.txt` (65 lines) — v05 1deg Darwin3 simulation based on ECCOV4r5 set-up with RADIv2 metamodel — `Darwin3`, `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `v4r4`, `V4r5`, `LLC90`, `git clone --branch backport_ckpt68g`, `git clone --depth 1`
- `v05/1deg_RADIv2_Nutrient_Kplnkt/readme_darwin_v4r5.txt` (63 lines) — v05 1deg Darwin3 simulation based on ECCOV4r5 set-up with RADIv2 metamodel — `Darwin3`, `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `v4r4`, `V4r5`, `LLC90`, `git clone --branch backport_ckpt68g`, `git clone --depth 1`
- `v05/1deg_oaemip/README.txt` (63 lines) — v05 1deg Darwin3 simulation based on ECCOV4r5 set-up with OAE experiments — `Darwin3`, `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `v4r4`, `V4r5`, `LLC90`, `git clone --branch backport_ckpt68g`, `git clone --depth 1`
- `v05/1deg_runoff/readme_darwin_v4r5.txt` (60 lines) — v05 1deg Darwin3 simulation based on ECCOV4r5 set-up with daily point source — `Darwin3`, `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `v4r4`, `V4r5`, `LLC90`, `git clone --branch backport_ckpt68g`, `git clone --depth 1`
- `v05/3deg/readme.txt` (68 lines) — v05 3deg darwin3 verification experiment with volume, salt, salinity, DIC, and FeT budget — `darwin3`, `darwin3/pkg/darwin`, `h_mpi`, `linux_amd64_ifort+mpi_ice_nas`, `data_mpi`, `llc270`, `git clone https://github.com/darwinproject/darwin3 cd`, `git checkout 24885b71 cd`, `git clone --depth 1`
- `v05/3deg/readme_sed.txt` (70 lines) — v05 3deg darwin3 verification experiment with volume, salt, salinity, DIC, and FeT budget — `darwin3`, `darwin3/pkg/darwin`, `h_mpi`, `linux_amd64_ifort+mpi_ice_nas`, `data_mpi`, `llc270`, `git clone https://github.com/darwinproject/darwin3 cd`, `git checkout 24885b71 cd`, `git clone --depth 1`
- `v05/3deg/useful/macOS_catalina/readme.txt` (5 lines) — For compiling on local OSX machine w/ Catalina or later use — `darwin3/eesupp/src`, `darwin_amd64_gfortran`
- `v05/3deg_CDR/readme_3deg_CDR.txt` (75 lines) — v05 3deg Darwin3 setup for Carbon Dioxide Removal (CDR) simulations — `Darwin3`, `darwin3`, `darwin3/pkg/darwin`, `linux_amd64_ifort+mpi_ice_nas`, `data_mpi`, `llc270`, `git clone https://github.com/darwinproject/darwin3 cd`, `git checkout 24885b71 cd`, `git clone --depth 1`
- `v05/3deg_CDR_MacroA/readme_3deg_CDR_MacroA.txt` (76 lines) — v05 3deg Darwin3 setup for Carbon Dioxide Removal (CDR) simulations — `Darwin3`, `darwin3`, `darwin3/pkg/darwin`, `linux_amd64_ifort+mpi_ice_nas`, `data_mpi`, `llc270`, `git clone https://github.com/darwinproject/darwin3 cd`, `git checkout 24885b71 cd`, `git clone --depth 1`
- `v05/llc270/readme.txt` (109 lines) — for regular ECCO-Darwin (w/ pkg/ctrl) — `checkpoint67d`, `darwin3`, `darwin3-dev`, `linux_amd64_ifort+mpi_ice_nas`, `LLC270`, `llc270`, `git clone --depth 1`, `git clone https://github.com/darwinproject/darwin3 cd`, `git checkout 24885b71 mkdir`, `git clone git://gud.mit.edu/darwin3-dev darwin3`
- `v05/llc270/readme2.txt` (46 lines) — for slim ECCO-Darwin (w/o pkg/ctrl) — `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `LLC270`, `llc270`, `git clone --depth 1`, `git clone https://github.com/darwinproject/darwin3 cd`, `git checkout 24885b71 mkdir`
- `v05/llc270/readme_1985.txt` (121 lines) — Instructions for building and running ECCO-Darwin v05 with Darwin 3 — `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `llc270`, `git clone --depth 1`, `git clone https://github.com/darwinproject/darwin3 cd`, `git checkout 24885b71 mkdir`
- `v05/llc270/readme_v5r1.txt` (48 lines) — Instructions for building and running ECCO-Darwin v05 with Darwin 3 — `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `llc270`, `git clone --depth 1`, `git clone https://github.com/darwinproject/darwin3 cd`, `git checkout darwin_ckpt68g mkdir`
- `v05/llc270_RADIv1/readme.txt` (112 lines) — iter42/input is available at https://data.nas.nasa.gov/ecco/data.php?dir=/eccodata/llc_270/iter42/input — `checkpoint67d`, `Darwin3`, `darwin3`, `darwin3-dev`, `linux_amd64_ifort+mpi_ice_nas`, `llc270`, `git clone --depth 1`, `git clone https://github.com/darwinproject/darwin3 cd`, `git checkout 24885b71 mkdir`, `git clone git://gud.mit.edu/darwin3-dev darwin3`
- `v05/llc270_jra55do/readme.txt` (85 lines) — iter42/input is available at https://data.nas.nasa.gov/ecco/data.php?dir=/eccodata/llc_270/iter42/input — `checkpoint67d`, `darwin3`, `darwin3-dev`, `linux_amd64_ifort+mpi_ice_nas`, `LLC270`, `llc270`, `git clone --depth 1`, `git clone https://github.com/darwinproject/darwin3 cd`, `git checkout 24885b71 mkdir`, `git clone git://gud.mit.edu/darwin3-dev darwin3`
- `v05/llc270_jra55do_mangroves/readme.txt` (95 lines) — Experimental set-up for mangrove carbon export based on LLC270 jra55-do physical solution. — `checkpoint67d`, `darwin3`, `darwin3-dev`, `linux_amd64_ifort+mpi_ice_nas`, `LLC270`, `llc270`, `git clone --depth 1`, `git clone https://github.com/darwinproject/darwin3 cd`, `git checkout 24885b71 mkdir`, `git clone git://gud.mit.edu/darwin3-dev darwin3`
- `v05/llc270_jra55do_nutrients/readme.txt` (89 lines) — iter42/input is available at https://data.nas.nasa.gov/ecco/data.php?dir=/eccodata/llc_270/iter42/input — `checkpoint67d`, `darwin3`, `darwin3-dev`, `linux_amd64_ifort+mpi_ice_nas`, `LLC270`, `llc270`, `git clone --depth 1`, `git clone https://github.com/darwinproject/darwin3 cd`, `git checkout 24885b71 mkdir`, `git clone git://gud.mit.edu/darwin3-dev darwin3`
- `v05/llc270_oaemip/README.txt` (63 lines) — OAEMIP v05 LLC270 Darwin3 simulation based on ECCOV5r1 set-up with OAE experiments — `Darwin3`, `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `LLC270`, `llc270`, `git clone --branch backport_ckpt68g`, `git clone --depth 1`
- `v06/1deg/readme_darwin_v4r5.txt` (50 lines) — v06 1deg Darwin3 simulation based on ECCOV4r5 set-up — `Darwin3`, `darwin3`, `darwin3/pkg/darwin`, `linux_amd64_ifort+mpi_ice_nas`, `V4r5`, `llc270`, `git clone https://github.com/darwinproject/darwin3 git`, `git checkout backport_ckpt68y cd`
- `v06/1deg/readme_darwin_v4r5_AO.txt` (52 lines) — v06 1deg Darwin3 simulation based on ECCOV4r5 set-up — `Darwin3`, `darwin3`, `darwin3/pkg/darwin`, `linux_amd64_ifort+mpi_ice_nas`, `V4r5`, `llc270`, `git clone https://github.com/darwinproject/darwin3 git`, `git checkout backport_ckpt68y cd`
- `v06/1deg/readme_v4r5_v2.txt` (40 lines) — ECCOV4r5 set-up — `checkpoint68y`, `linux_amd64_ifort+mpi_ice_nas`, `git clone https://github.com/MITgcm-contrib/ecco_darwin.git git`, `git checkout checkpoint68y #`
- `v06/3deg/readme.txt` (90 lines) — v06 3deg darwin3 verification experiment — `darwin3`, `darwin3/pkg/darwin`, `darwin3/run`, `h_mpi`, `linux_amd64_ifort+mpi_ice_nas`, `data_mpi`, `llc270`, `git clone https://github.com/darwinproject/darwin3 git`, `git checkout backport_ckpt68y cd`
- `v06/llc270/readme_darwin.txt` (41 lines) — Instructions for building and running v06 ECCO-Darwin — `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `llc270`, `git clone https://github.com/darwinproject/darwin3 git`
- `v06/llc270/readme_physics.txt` (35 lines) — Instructions for building and running v06 ECCO-Darwin llc270 physics — `darwin3`, `linux_amd64_ifort+mpi_ice_nas`, `llc270`, `git clone https://github.com/darwinproject/darwin3 git`

# @darwin3_manual = `~/Documents/GitHub/darwin3/doc` (28ac947a3 2024-09-24 Revert scavenging defaults to POP-based equivalents)

## `phys_pkgs/darwin.rst`
- L5 DARWIN package
  - L113 Compiling and Running
    - L116 Compiling — `GCHEM_SEPARATE_FORCING`, `ALLOW_CLIMSST_RELAXATION`, `ALLOW_CLIMSSS_RELAXATION`, `DARWIN_ALLOW_NQUOTA`, `DARWIN_ALLOW_PQUOTA`, `DARWIN_ALLOW_FEQUOTA`, `DARWIN_ALLOW_SIQUOTA`, `DARWIN_ALLOW_CHLQUOTA`, `DARWIN_ALLOW_CDOM`, `DARWIN_CDOM_UNITS_CARBON`, `DARWIN_ALLOW_CSTORE`, `DARWIN_ALLOW_CSTORE_DIAGS`, `DARWIN_ALLOW_CARBON`, `DARWIN_SOLVESAPHE`, `DARWIN_TOTALPHSCALE`, `DARWIN_USE_PLOAD`, `DARWIN_ALLOW_RADI`, `DARWIN_ALLOW_DENIT`, `DARWIN_ALLOW_EXUDE`, `ALLOW_OLD_VIRTUALFLUX`, `DARWIN_NITRATE_FELIMIT`, `DARWIN_BOTTOM_SINK`, `DARWIN_NUTRIENT_RUNOFF`, `DARWIN_AVPAR`, `DARWIN_ALLOW_GEIDER`, `DARWIN_GEIDER_RHO_SYNTH`, `DARWIN_CHL_INIT_LEGACY`, `DARWIN_SCATTER_CHL`, `DARWIN_DIAG_IOP`, `DARWIN_GRAZING_SWITCH` (+18)
    - L210 Running
      - L216 Runtime Parameters — `DARWIN_FORCING_PARAMS`, `DARWIN_INTERP_PARAMS`, `USE_EXF_INTERPOLATION`, `DARWIN_PARAMS`, `DARWIN_CDOM_PARAMS`, `DARWIN_RADTRANS_PARAMS`, `DARWIN_RANDOM_PARAMS`, `DARWIN_TRAIT_PARAMS`, `darwin_useQsw`, `darwin_useSEAICE`, `darwin_useEXFwind`, `alpfe`, `PARfile`, `PARconst`, `PARperiod`, `PARRepCycle`, `repeatPeriod`, `useExfYearlyFields`, `PARStartTime`, `PARstartdate1`, `PARstartdate2`, `PAR_exfremo_intercept`, `PAR_exfremo_slope`, `PARmask`, `darwin_inscal_PAR`, `R_ALK_DIC_runoff`, `R_NO3_DIN_runoff`, `R_NO2_DIN_runoff`, `R_NH4_DIN_runoff`, `R_DIP_IP_runoff` (+147)
      - L536 Traits — `isPhoto`, `bactType`, `isAerobic`, `isDenit`, `hasSi`, `hasPIC`, `diazo`, `useNH4`, `useNO2`, `useNO3`, `combNO`, `isPrey`, `isPred`, `tempMort`, `tempMort2`, `tempGraz`, `Xmin`, `amminhib`, `acclimtimescl`, `mort`, `mort2`, `ExportFracMort`, `ExportFracMort2`, `ExportFracExude`, `phytoTempCoeff`, `phytoTempExp1`, `phytoTempAe`, `phytoTempExp2`, `phytoTempOptimum`, `phytoDecayPower` (+70)
      - L670 Allometric trait generation — `isPhoto`, `grp_photo`, `bactType`, `grp_bacttype`, `isAerobic`, `grp_aerobic`, `isDenit`, `grp_denit`, `isPred`, `grp_pred`, `isPrey`, `grp_prey`, `hasSi`, `grp_hasSi`, `hasPIC`, `grp_hasPIC`, `diazo`, `grp_diazo`, `useNH4`, `grp_useNH4`, `useNO2`, `grp_useNO2`, `useNO3`, `grp_useNO3`, `combNO`, `grp_combNO`, `aptype`, `grp_aptype`, `tempMort`, `grp_tempMort` (+208)
  - L838 Diagnostics — `SMR_____MR`, `SMRP____MR`, `SM_P____MR`, `SM_P____LR`, `SM_P____L1`, `SM_P____U1`, `SM______L1`
  - L958 Call Tree

## `phys_pkgs/darwin_airsea.rst`
- L5 Air-sea exchanges — `DARWIN_ALLOW_CARBON`, `icefile`, `useSEAICE`, `m3perkg`, `darwin_useSEAICE`

## `phys_pkgs/darwin_bacteria.rst`
- L5 Bacteria — `bactType`, `isAerobic`, `isDenit`, `grp_bacttype`, `grp_aerobic`, `grp_denit`, `KDOC`
  - L38 Growth and energy sources
  - L74 Generic particle-associated
  - L118 Generic free-living
  - L155 Bacteria parameters — `bactType`, `grp_bacttype`, `isAerobic`, `grp_aerobic`, `isDenit`, `grp_denit`, `pcoefO2`, `pmaxDIN`, `ksatDIN`, `alpha_hydrol`, `PCmax`, `a`, `b_PCmax`, `yield`, `yod`, `ynd`, `yieldO2`, `yoe`, `yieldNO3`, `yne`, `ksatPON`, `a_ksatPON`, `ksatPOC`, `ksatPOP`, `ksatPOFe`, `ksatDON`, `a_ksatDON`, `ksatDOC`, `ksatDOP`, `ksatDOFe`

## `phys_pkgs/darwin_carbon.rst`
- L5 Carbon chemistry
  - L8 Carbon chemistry options — `DARWIN_ALLOW_CARBON`, `DARWIN_SOLVESAPHE`, `DARWIN_TOTALPHSCALE`, `DARWIN_ALLOW_RADI`, `selectPHsolver`, `surfSaltMinInit`, `surfSaltMaxInit`, `surfTempMinInit`, `surfTempMaxInit`, `surfDICMinInit`, `surfDICMaxInit`, `surfALKMinInit`, `surfALKMaxInit`, `surfPO4MinInit`, `surfPO4MaxInit`, `surfSiMinInit`, `surfSiMaxInit`, `surfSaltMin`, `surfSaltMax`, `surfTempMin`, `surfTempMax`, `surfDICMin`, `surfDICMax`, `surfALKMin`, `surfALKMax`, `surfPO4Min`, `surfPO4Max`, `surfSiMin`, `surfSiMax`
  - L134 Calcite dissolution — `darwin_disscSelect`, `DARWIN_SOLVESAPHE`, `Kdissc`, `darwin_KeirCoeff`, `darwin_KeirExp`
  - L202 Diagnostics — `DARWIN_ALLOW_CARBON`, `DARWIN_ALLOW_RADI`, `SM_P____L1`, `SMR_____MR`, `SMRP____MR`, `SM______L1`, `SM______U1`, `SM_P____U1`, `SM_P____M1`

## `phys_pkgs/darwin_cdom.rst`
- L5 Dynamic CDOM — `DARWIN_ALLOW_CDOM`, `DARWIN_CDOM_UNITS_CARBON`, `DARWIN_ALLOW_DENIT`, `fracCDOM`, `CDOMdegrd`, `CDOMbleach`, `PARCDOM`, `R_NP_CDOM`, `R_FeP_CDOM`, `R_CP_CDOM`, `R_NC_CDOM`, `R_PC_CDOM`, `R_FeC_CDOM`

## `phys_pkgs/darwin_changes.rst`
- L3 Change Log

## `phys_pkgs/darwin_chl.rst`
- L5 Chlorophyll synthesis
  - L8 With N quota:
  - L27 Without N quota: — `DARWIN_GEIDER_RHO_SYNTH`
  - L60 Without Chl quota, — `chl2nmax`, `acclimtimescl`, `a_acclimtimescl`

## `phys_pkgs/darwin_cons.rst`
- L5 Conservation of chemical elements — `DARWIN_ALLOW_CONS`, `DARWIN_BOTTOM_SINK`, `darwin_linFSConserv`, `DARWIN_ALLOW_CARBON`, `ALLOW_OLD_VIRTUALFLUX`, `ironFile`, `DARWIN_MINFE`, `freefemax`, `DARWIN_PART_SCAV`, `fesedflux`, `fesedflux_pcm`, `DARWIN_IRON_SED_SOURCE_VARIABLE`, `DARWIN_ALLOW_DENIT`, `DARWIN_ALLOW_NQUOTA`
  - L72 Linear free surface — `implicitFreeSurface`, `DARWIN_linFSConserve`

## `phys_pkgs/darwin_cstore.rst`
- L5 Internal carbon store and exudation — `DARWIN_ALLOW_NQUOTA`, `DARWIN_ALLOW_PQUOTA`, `DARWIN_ALLOW_FEQUOTA`, `DARWIN_ALLOW_SIQUOTA`, `DARWIN_ALLOW_CSTORE`, `FracExudeC`, `a_FracExudeC`

## `phys_pkgs/darwin_denit.rst`
- L5 Denitrification — `DARWIN_ALLOW_DENIT`, `DARWIN_ALLOW_CDOM`, `denit_NP`, `denit_NO3`, `O2crit`, `NO3crit`

## `phys_pkgs/darwin_equations.rst`
- L3 Model equations — `DARWIN_ALLOW_CDOM`, `DARWIN_ALLOW_CARBON`, `R_PICPOC`, `a_R_PICPOC`, `R_OP`

## `phys_pkgs/darwin_exude.rst`
- L3 Exudation — `DARWIN_ALLOW_EXUDE`, `kexcc`, `a_kexcC`, `kexcn`, `a_kexcN`, `kexcp`, `a_kexcP`, `kexcfe`, `a_kexcFe`, `kexcsi`, `a_kexcSi`, `ExportFracExude`, `a_ExportFracExude`

## `phys_pkgs/darwin_grazing.rst`
- L5 Grazing — `DARWIN_GRAZING_SWITCH`, `hillnumGraz`
  - L114 Implementation
  - L191 Runtime Parameters — `grazemax`, `a_grazemax`, `kgrazesat`, `a_kgrazesat`, `tempGraz`, `grp_tempGraz`, `inhib_graz`, `inhib_graz_exp`, `hillnumGraz`, `hollexp`, `phygrazmin`, `palat`, `asseff`, `grp_ass_eff`, `ExportFracPreyPred`, `grp_ExportFracPreyPred`, `DARWIN_ALLOMETRIC_PALAT`, `grp_pred`, `grp_prey`, `a`, `b_ppOpt`, `a_ppSig`, `palat_min`

## `phys_pkgs/darwin_growth.rst`
- L5 Growth — `DARWIN_ALLOW_GEIDER`, `aphy_chl_ave`, `DARWIN_ALLOW_CHLQUOTA`, `ksatPAR`, `a_ksatPAR`, `kinhPAR`, `a_kinhPAR`, `PCmax`, `a`, `b_PCmax`, `PARmin`, `mQyield`, `a_mQyield`, `chl2cmax`, `a_chl2cmax`, `inhibGeider`, `a_inhibGeider`, `aphy_chl_ps`, `aphy_chl_ps_type`, `grp_aptype`, `darwin_phytoAbsorbFile`

## `phys_pkgs/darwin_iron.rst`
- L5 Iron chemistry
  - L24 Dust deposition — `ironfile`, `alpfe`, `darwin_inscal_iron`
  - L34 Sedimentation — `depthfesed`, `fesedflux`, `DARWIN_IRON_SED_SOURCE_VARIABLE`, `DARWIN_IRON_SED_SOURCE_POP`
  - L62 Hydrothermal vents — `DARWIN_ALLOW_HYDROTHERMAL_VENTS`, `depthFeVent`, `ventHe3file`, `solFeVent`, `R_FeHe3_vent`
  - L82 Scavenging — `DARWIN_PART_SCAV`, `DARWIN_PART_SCAV_POP`, `DARWIN_MINFE`, `alpfe`, `depthfesed`, `fesedflux`, `fesedflux_pcm`, `DARWIN_IRON_SED_SOURCE_VARIABLE`, `fesedflux_min`, `R_CP_fesed`, `depthFeVent`, `solFeVent`, `R_FeHe3_vent`, `scav`, `scav_tau`, `scav_inter`, `scav_exp`, `scav_POC_wgt`, `scav_PSi_wgt`, `scav_PIC_wgt`, `scav_rPOM`, `ligand_tot`, `ligand_stab`, `freefemax`, `scav_rat`, `scav_R_POPPOC`

## `phys_pkgs/darwin_light.rst`
- L5 Non-spectral Light — `PARfile`, `PARconst`, `darwin_useQsw`, `icefile`, `DARWIN_AVPAR`, `DARWIN_ALLOW_CHLQUOTA`, `DARWIN_ALLOW_GEIDER`, `R_ChlC`, `parfrac`, `parconv`, `katten_w`, `katten_Chl`, `aphy_chl_ave`

## `phys_pkgs/darwin_mort.rst`
- L5 Mortality
  - L31 Parameters — `mort`, `a_mort`, `mort2`, `a_mort2`, `Xmin`, `a_Xmin`, `tempMort`, `grp_tempMort`, `tempMort2`, `grp_tempMort2`, `ExportFracMort`, `a_ExportFracMort`, `ExportFracMort2`, `a_ExportFracMort2`

## `phys_pkgs/darwin_remin.rst`
- L5 Remineralization and Nitrification — `DARWIN_ALLOW_DENIT`, `KDOC`, `KDOP`, `KDON`, `KDOFe`, `KPOC`, `KPOP`, `KPON`, `KPOFe`, `KPOSi`, `Knita`, `Knitb`, `PAR_oxi`

## `phys_pkgs/darwin_resp.rst`
- L5 Respiration
  - L46 Parameters — `respRate`, `a`, `b_respRate_c`, `qcarbon`, `b_qcarbon`, `Xmin`, `a_Xmin`, `a_respRate_c`

## `phys_pkgs/darwin_sink.rst`
- L5 Sinking and Swimming — `DARWIN_BOTTOM_SINK`, `wPIC_sink`, `wC_sink`, `wN_sink`, `wP_sink`, `wSi_sink`, `wFe_sink`, `biosink`, `a`, `b_biosink`, `bioswim`, `b_bioswim`

## `phys_pkgs/darwin_spectral.rst`
- L5 Spectral Light — `darwin_waterAbsorbFile`, `grp_aptype`, `darwin_phytoAbsorbFile`, `DARWIN_SCATTER_CHL`, `darwin_particleAbsorbFile`, `DARWIN_ALLOW_CDOM`, `darwin_bbmin`, `darwin_bbw`, `darwin_RPOC`, `darwin_rCDOM`, `CDOMcoeff`, `darwin_lambda_aCDOM`, `darwin_Sdom`, `darwin_aCDOM_fac`, `darwin_part_size_P`, `aphy_chl`, `aphy_chl_ps`, `aphy_mgC`, `bphy_mgC`, `bbphy_mgC`
  - L121 Format of optical spectra files — `grp_aptype`, `darwin_waterAbsorbFile`, `darwin_particleAbsorbFile`, `darwin_phytoAbsorbFile`
  - L157 Allometric scaling of absorption and scattering spectra — `darwin_allomSpectra`
    - L169 Absorption
    - L192 Total scattering
    - L223 Backscattering — `darwin_allomSpectra`, `darwin_aCarCell`, `darwin_bCarCell`, `darwin_absorpSlope`, `darwin_bbbSlope`, `darwin_scatSwitchSizeLog`, `darwin_scatSlopeSmall`, `darwin_scatSlopeLarge`
  - L265 Photosynthetically Active Radation

## `phys_pkgs/darwin_temperature.rst`
- L5 Temperature dependence — `DARWIN_TEMP_VERSION`, `DARWIN_TEMP_RANGE`, `DARWIN_NOTEMP`, `tempMort`, `tempMort2`, `tempGraz`
  - L18 DARWIN_TEMP_VERSION 1 — `DARWIN_TEMP_RANGE`
  - L39 DARWIN_TEMP_VERSION 2 — `DARWIN_TEMP_RANGE`
  - L66 DARWIN_TEMP_VERSION 3
  - L85 DARWIN_TEMP_VERSION 4 — `DARWIN_TEMP_RANGE`, `phytoTempCoeff`, `a_phytoTempCoeff`, `phytoTempExp1`, `a_phytoTempExp1`, `tempnorm`, `TempCoeffArr`, `TempAeArr`, `TempRefArr`, `phytoTempAe`, `a_phytoTempAe`, `hetTempAe`, `a_hetTempAe`, `grazTempAe`, `a_grazTempAe`, `reminTempAe`, `mortTempAe`, `mort2TempAe`, `uptakeTempAe`, `phytoTempExp2`, `a_phytoTempExp2`, `phytoTempOptimum`, `a_phytoTempOptimum`, `phytoDecayPower`, `a_phytoDecayPower`, `hetTempExp2`, `a_hetTempExp2`, `hetTempOptimum`, `a_hetTempOptimum`, `hetDecayPower` (+9)

## `phys_pkgs/darwin_uptake.rst`
- L5 Nutrient uptake and limitation
  - L31 Without P quota:
  - L40 With P quota:
  - L66 Si: — `hasSi`
  - L79 Without N quota:
    - L82 diazotroph:
    - L91 not diazotroph: — `combNO`
  - L146 With N quota:
    - L185 diazotroph:
    - L204 not diazotroph:
  - L210 Without Fe quota:
  - L220 With Fe quota,
  - L253 Effective half saturation constants — `DARWIN_effective_ksat`, `darwin_select_kn_allom`, `a_ksatNO3`, `b_ksatNO3`
  - L287 Uptake and limitation parameters — `synthcost`, `hasSi`, `grp_hasSi`, `diazo`, `grp_diazo`, `useNH4`, `grp_useNH4`, `useNO2`, `grp_useNO2`, `useNO3`, `grp_useNO3`, `combNO`, `grp_combNO`, `Qnmin`, `a`, `b_Qnmin`, `Qnmax`, `b_Qnmax`, `Qpmin`, `b_Qpmin`, `Qpmax`, `b_Qpmax`, `Qsimin`, `b_Qsimin`, `Qsimax`, `b_Qsimax`, `Qfemin`, `b_Qfemin`, `Qfemax`, `b_Qfemax` (+41)

## `phys_pkgs/radtrans.rst`
- L5 RADTRANS package
  - L8 Introduction
  - L23 Compiling and Running
    - L26 Compiling — `RT_useMeanCosSolz`, `ALLOW_CLIMSST_RELAXATION`, `ALLOW_CLIMSSS_RELAXATION`, `RADTRANS_DIAG_SOLUTION`
    - L52 Running — `useRADTRANS`
      - L57 Runtime Parameters — `RADTRANS_FORCING_PARAMS`, `RADTRANS_INTERP_PARAMS`, `USE_EXF_INTERPOLATION`, `RADTRANS_PARAMS`, `RT_Edfile`, `RT_Esfile`, `RT_Ed_const`, `RT_Es_const`, `RT_E_period`, `RT_E_RepCycle`, `repeatCycle`, `useOasimYearlyFields`, `RT_E_StartTime`, `RT_E_startdate1`, `RT_E_startdate2`, `RT_E_mask`, `RT_Ed_exfremo_intercept`, `RT_Es_exfremo_intercept`, `RT_Ed_exfremo_slope`, `RT_Es_exfremo_slope`, `RT_inscal_Ed`, `RT_inscal_Es`, `RT_refract_water`, `RT_rmud_max`, `RT_wbEdges`, `RT_wbRefWLs`, `RT_kmax`, `RT_useMeanCosSolz`, `RT_sfcIrrThresh`
  - L130 Diagnostics — `SM_P____L1`, `SMRP____LR`, `SMRP____MR`, `SMR_____MR`, `SM_P____MR`
  - L177 Call tree
