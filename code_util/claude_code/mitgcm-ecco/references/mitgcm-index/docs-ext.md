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

# @ecco_darwin = `~/Documents/GitHub/ecco_darwin` (c0f1de2 2026-10-19 Update 1-D time series stations set-ups)
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

## `code_util/claude_code/mitgcm-ecco/SKILL.md`
- L6 MITgcm + ECCO + Darwin
  - L19 Ground rules (these come from real incidents) — `nTimeSteps`, `nEndIter`
  - L43 How to be the MITgcm expert in the room — `EXACT_CONSERV`, `backport_ckpt68y`
  - L80 Route to the right reference
  - L103 Fast facts worth always having in mind — `mitgcmuv_ad`, `mitgcmuv_ftl`, `SIZE.h`, `readBinaryPrec`, `writeBinaryPrec`, `available_diagnostics.log`

## `code_util/claude_code/mitgcm-ecco/references/analysis.md`
- L1 Reading and plotting MITgcm / ECCO output
  - L3 Formats — `dimList`, `fldList`
  - L14 Python (preferred) — `ecco_v4_py`, `scipy.spatial.cKDTree`
  - L33 MATLAB
  - L40 Analysis conventions the user wants — `dumpAtLast`, `cDIC`, `C_DIC`, `cDIC_PIC`, `C_DICPIC`, `freeFe`, `freeFeLs`
  - L54 PO.DAAC / NetCDF publishing — `standard_name`, `long_name`

## `code_util/claude_code/mitgcm-ecco/references/build.md`
- L1 Building MITgcm / darwin3
  - L3 Standard sequence — `PACKAGES_CONFIG.h`
  - L27 Local Mac optfiles (arm64 Mac) — `_timenow_`, `numDiags`, `numLevels`
  - L45 Pleiades (NAS) optfiles — `optfile_toss4`
  - L58 Adjoint / tangent linear — `mitgcmuv_ad`, `mitgcmuv_ftl`, `AD_TAF_FLAGS`, `FTL_TAF_FLAGS`, `MDSIO_BUFF_WH.h`, `CTRL_INIT_VARIABLES`, `mult_gencost`

## `code_util/claude_code/mitgcm-ecco/references/config-consistency.md`
- L1 Configuration sanity: code_* + run_* pairs
  - L7 Compile-time ↔ run-time couplings — `SIZE.h`, `delZ`, `tRef`, `sRef`, `useCDscheme`, `PTRACERS_num`, `PTRACERS_SIZE.h`, `PTRACERS_numInUse`, `nGroup`, `nPhoto`, `DARWIN_SIZE.h`, `grp_nplank`, `grp_photo`, `grp_names`, `DARWIN_ALLOW_CHLQUOTA`, `RADTRANS_SIZE.h`, `RT_EdFile`, `RT_EsFile`, `RT_wbRefWLs`, `RT_wbEdges`, `DIAGNOSTICS_SIZE.h`, `numDiags`, `MNC_SIZE.h`, `RW_MFLDS.h`, `CPP_OPTIONS.h`, `USE_EXF_INTERPOLATION`, `DARWIN_INDICES.h`, `PTRACERS_initialFile`
  - L26 Time and calendar — `startDate_1`, `calendarDumps`, `deltaToffline`, `offlineIter0`, `offlineForcingPeriod`, `offlineForcingCycle`, `nIter0`, `data_org`, `repeatPeriod`, `EXF_GetFFieldRec`, `USE_EXF_INTERPOLATION`
  - L49 Namelist traps — `viscAz`, `viscAr`, `diag_mnc`, `useMNC`, `offlineForcingCycle`, `momStepping`, `tempStepping`, `ggl90_check.F`, `GAD_SMOLARKIEWICZ_HACK`, `EXF_INTERP_READ`, `MDS_READ_FIELD`, `RT_icefile`
  - L67 Offline runs must mirror their forward run — `hFacMin`, `hFacMinDr`, `bathyFile`, `useShelfIce`, `diffKrFile`, `ALLOW_3D_DIFFKR`, `xx_kapredi`
  - L83 Paths and inputs — `readBinaryPrec`, `available_diagnostics.log`, `fileName`, `stat_freq`, `deltaT`
  - L93 Before a long run — `nTimeSteps`

## `code_util/claude_code/mitgcm-ecco/references/configs.md`
- L1 Where things live: configurations and projects
  - L6 Local MITgcm / darwin3 clones (keep separate; each has its own branch work) — `ALLOW_OBCS_SEAICE_SPONGE`, `backport_ckpt68y`
  - L20 Configurations — `run_v4r6_forward`, `backport_ckpt68y`, `ggl90_check.F`, `COMPILE_AND_RUN.md`, `EXTERNAL_DATA_NEEDED.md`, `global_oce_latlon`, `global_ocean.cs32x15`
  - L37 Data sources — `SEAICEnonLinIterMax`

## `code_util/claude_code/mitgcm-ecco/references/downscaling.md`
- L1 Regional downscaling (ECCO-Darwin parent → regional MITgcm/darwin3)
  - L9 Where things live — `readme_pleiades.txt`, `chain_watchdog_selftest.sh`
  - L24 Regions (snapshot 2026-10-06) — `kerguelen_4km`, `job_production_chain_r2688.sh`, `OBCS_1700_preext`, `pan_AO`, `west_AO`, `GoM_jra55do_nutrients`, `mac_delta`, `disko_bay`, `nares_strait`, `CCS_surf`, `pan_greenland`
  - L36 The pipeline (STEP1–4) and the one rule that matters — `run_step1_mitgrid.sh`, `job_step1_bathy.sh`, `job_step1_stitch.sh`, `DIAGNOSTICS_VEC_SIZE.h`, `nVEC_mask`, `nml_vecFiles`, `data.diagnostics_vec`, `VEC_points`, `gen_obcs_fast.py`, `gen_pickups.py`, `hydrogThetaFile`, `gen_derived_tracer_obcs.py`, `_pickup.py`, `gen_seaice_ic.py`, `gen_obcs_transport_correction.py`, `stage_regional_run.sh`, `setup_run.sh`, `stage_run_4km.sh`, `job_production_chain.sh`, `find_restart_iter.py`
  - L64 Working rules (on top of the SKILL.md ground rules) — `make_pan_AO_grid_nc.py`
  - L78 Failure modes seen (symptom → cause → fix) — `job_build.sh`, `code_darwin`, `rom_ait`, `bro_ele`, `nTimeSteps`, `numDiags`, `DIAGNOSTICS_SIZE.h`, `oasim_calc_solz`, `oasim_dTsolz`, `deltaT`, `make_4km.py`, `readBinaryPrec`, `useOBCSbalance`, `eco_tracers.SEED_FRAC`, `ZOO_SEED_FRAC`, `chain_watchdog.sh`, `time_tsnumber`
  - L113 Analysis recipes

## `code_util/claude_code/mitgcm-ecco/references/hpc.md`
- L1 Running: NASA Pleiades and the local Mac
  - L3 Pleiades access (from Claude Code)
  - L21 PBS job script conventions (csh) — `sky_ele`, `nIter0`
  - L53 Remote paths — `v05_1deg_V4r5`, `_bbl`, `sea_ice_BCs_latlon`, `from_steph`
  - L66 Local Mac runs

## `code_util/claude_code/mitgcm-ecco/references/literature.md`
- L1 Literature behind ECCO, ECCO-Darwin and Dustin's MITgcm work
  - L11 ECCO state estimation (framework, releases)
  - L28 ECCO-Darwin and BGC state estimation — `llc270_RADIv1`
  - L49 Regional / river-plume configurations (Mackenzie)
  - L57 Glacier fjords and subglacial plumes (MITgcm, pre-ECCO work)
  - L66 Adjoint tools
  - L72 Updating

## `code_util/claude_code/mitgcm-ecco/references/package-development.md`
- L1 Developing MITgcm and Darwin packages
  - L15 Anatomy of a package `pkg/<name>/` — `PACKAGES_CONFIG.h`, `CPP_OPTIONS.h`, `OPEN_COPY_DATA_FILE`, `_BEGIN_MASTER`, `_END_MASTER`, `_BARRIER`, `ALLOW_DIAGNOSTICS`, `DIAGNOSTICS_ADDTOLIST`, `forward_step.F`, `external_forcing_surf.F`, `do_oceanic_phys.F`, `packages_readparms.F`, `packages_check.F`, `PACKAGES_ERROR_MSG`, `packages_init_fixed.F`, `packages_init_variables.F`, `packages_write_pickup.F`, `phys_pkgs.rst`
  - L43 Fortran conventions (match them exactly) — `_EXCH_XYZ_RL`, `kLowC`, `kSurfC`, `select_rStar`, `READ_FLD_XY_RL`, `READ_REC_3D_RL`, `WRITE_FLD_XY_RL`, `GLOBAL_SUM_TILE_RL`
  - L72 Diagnostics — `available_diagnostics.log`
  - L86 Darwin (darwin3) specifics — `gchem_readparms`, `gchem_init_fixed`, `gchem_init_vari`, `gchem_tr_register`, `gchem_forcing_sep`, `DARWIN_FORCING`, `gchem_output`, `DARWIN_DIAGS`, `gchem_write_pickup`, `GCHEM_SEPARATE_FORCING`, `DARWIN_SIZE.h`, `DARWIN_INDICES.h`, `DARWIN_PARAMS.h`, `DARWIN_TRAITPARAMS.h`, `DARWIN_TRAITS.h`, `darwin_generate_allometric.F`, `darwin_generate_random.F`, `darwin_light_radtrans.F`, `darwin_light.F`, `darwin_carbon_chem.F`, `darwin_solvesaphe.F`, `backport_ckpt68y`, `backport_ckpt68g`, `DARWIN_ALLOW_DIAZ`, `darwin_tr_register.F`, `darwin_cons.F`
  - L114 Adjoint-readiness — `WAD_SMOOTH_MASK`, `MDSIO_BUFF_WH.h`
  - L125 Testing a new feature (minimum before calling it done) — `BBL_CHECK_BUDGET`

## `code_util/claude_code/mitgcm-ecco/references/verification-upstream.md`
- L1 Verification, testreport and upstreaming
  - L3 testreport — `internal_wave`, `global_ocean.cs32x15`, `global_oce_latlon`, `check_wad_refs.sh`, `check_wad_refs_x86.sh`
  - L24 A new verification experiment
  - L37 Correctness checks the user expects
  - L47 Repositories — `backport_ckpt68y`, `backport_ckpt68g`
  - L56 Upstreaming decisions so far
  - L68 MITgcm PR conventions

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/README.md`
- L1 MITgcm source index — `backport_ckpt68y`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/bibliography.md`
- L1 Bibliography: papers behind MITgcm / darwin3 schemes — `manual_references.bib`
  - L8 Part 1: citations by manual page
    - L10 darwin3 `doc/phys_pkgs/darwin_airsea.rst`
    - L16 darwin3 `doc/phys_pkgs/darwin_carbon.rst`
    - L37 darwin3 `doc/phys_pkgs/darwin_iron.rst` — `DARWIN_PART_SCAV`
    - L42 darwin3 `doc/phys_pkgs/darwin_spectral.rst`
    - L47 `doc/algorithm/adv-schemes.rst`
    - L52 `doc/algorithm/algorithm.rst` — `viscAhReMax`, `viscC2Smag`, `ALLOW_LEITH_QG`
    - L81 `doc/algorithm/c-grid.rst`
    - L84 `doc/algorithm/finitevol-meth.rst`
    - L87 `doc/algorithm/nonlinear-freesurf.rst`
    - L92 `doc/algorithm/vert-grid.rst`
    - L95 `doc/autodiff/autodiff.rst`
    - L111 `doc/contributing/contributing.rst`
    - L116 `doc/examples/advection_in_gyre/advection_in_gyre.rst`
    - L120 `doc/examples/baroclinic_gyre/baroclinic_gyre.rst` — `no_slip_sides`
    - L132 `doc/examples/barotropic_gyre/barotropic_gyre.rst`
    - L145 `doc/examples/deep_convection/deep_convection.rst`
    - L148 `doc/examples/examples.rst`
    - L172 `doc/examples/global_oce_biogeo/global_oce_biogeo.rst`
    - L188 `doc/examples/global_oce_in_p/global_oce_in_p.rst`
    - L205 `doc/examples/global_oce_latlon/global_oce_latlon.rst`
    - L221 `doc/examples/global_oce_optim/global_oce_optim.rst`
    - L230 `doc/examples/held_suarez_cs/held_suarez_cs.rst`
    - L248 `doc/examples/reentrant_channel/reentrant_channel.rst`
    - L283 `doc/examples/tracer_adjsens/tracer_adjsens.rst`
    - L289 `doc/getting_started/getting_started.rst` — `ALLOW_BL79_LAT_VARY`, `diffKrBL79surf`, `diffKrBL79deep`, `diffKrBL79scl`, `diffKrBL79Ho`
    - L308 `doc/ocean_state_est/ocean_state_est.rst`
    - L319 `doc/overview/finding_pressure.rst`
    - L327 `doc/overview/global_atmos_hs.rst`
    - L331 `doc/overview/hydrostatic.rst`
    - L338 `doc/overview/overview.rst`
    - L352 `doc/overview/soln_strategy.rst`
    - L355 `doc/phys_pkgs/fizhi.rst`
    - L387 `doc/phys_pkgs/generic_advdiff.rst` — `GAD_ALLOW_TS_SOM_ADV`, `GAD_SMOLARKIEWICZ_HACK`
    - L391 `doc/phys_pkgs/ggl90.rst`
    - L394 `doc/phys_pkgs/gmredi.rst`
    - L413 `doc/phys_pkgs/gridalt.rst`
    - L416 `doc/phys_pkgs/kl10.rst`
    - L424 `doc/phys_pkgs/kpp.rst`
    - L430 `doc/phys_pkgs/obcs.rst`
    - L435 `doc/phys_pkgs/opps.rst`
    - L438 `doc/phys_pkgs/rbcs.rst`
    - L441 `doc/phys_pkgs/remesh.rst`
    - L445 `doc/phys_pkgs/seaice.rst` — `SEAICEuseEVPstar`, `SEAICEuseEVPrev`, `SEAICEaEVPcStar`, `SEAICEaEVPalphaMin`, `SEAICEbasalDragK1`, `SEAICESideDrag`, `useHibler79IceStrength`, `SEAICEsimpleRidging`, `SEAICE_cf`, `SEAICEpartFunc`, `SEAICEredistFunc`, `SEAICEgStar`, `SEAICEhStar`, `SEAICEuseLinRemapITD`, `Hlimit_c3`, `SEAICE_cStar`, `SEAICE_tensilFac`, `SEAICE_ALLOW_MCS`
    - L549 `doc/phys_pkgs/shelfice.rst` — `SHI_ALLOW_GAMMAFRICT`, `SHELFICEuseGammaFrict`
    - L566 `doc/phys_pkgs/streamice.rst`
    - L578 `doc/software_arch/software_arch.rst`
  - L582 Part 2: full references (`key`: citation; cited N times)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/core-calltree.md`
- L1 Core call tree (static) — `THE_MODEL_MAIN`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/core-params.md`
- L1 Core namelists
  - L5 model/src/ini_parms.F  (reads data)
    - L6 PARM01 — `gravitySign`, `nh_Am2`, `gBaro`, `gravityFile`, `tAlpha`, `sBeta`, `selectCoriMap`, `fPrime`, `rotationPeriod`, `viscAh`, `viscAhW`, `viscAhMax`, `viscAhGrid`, `viscAhGridMax`, `viscAhGridMin`, `viscC2leith`, `viscC4leith`, `smag3D_coeff`, `useSmag3D`, `useFullLeith`, `useAnisotropicViscAgridMax`, `useStrainTensionVisc`, `useAreaViscLength`, `viscC2leithD`, `viscC4leithD`, `viscC2LeithQG`, `viscC2smag`, `viscC4smag`, `viscAhD`, `viscAhZ` (+178)
    - L221 PARM02 — `cg2dMaxIters`, `cg2dMinItersNSA`, `cg2dChkResFreq`, `cg2dUseMinResSol`, `cg2dTargetResidual`, `cg2dTargetResWunit`, `cg2dpcOffDFac`, `cg2dPreCondFreq`, `cg3dMaxIters`, `cg3dChkResFreq`, `cg3dTargetResidual`, `cg3dTargetResWunit`, `useNSACGSolver`, `useSRCGSolver`, `printResidualFreq`
    - L237 PARM03 — `nIter0`, `nTimeSteps`, `nTimeSteps_l2`, `nEndIter`, `baseTime`, `startTime`, `endTime`, `deltaT`, `deltaTClock`, `deltaTMom`, `deltaTtracer`, `dTtracerLev`, `deltaTFreeSurf`, `forcing_In_AB`, `momForcingOutAB`, `tracForcingOutAB`, `momDissip_In_AB`, `doAB_onGtGs`, `abEps`, `alph_AB`, `beta_AB`, `startFromPickupAB2`, `applyExchUV_early`, `tauCD`, `epsAB_CD`, `cAdjFreq`, `chkPtFreq`, `pChkPtFreq`, `pickupSuff`, `pickupStrictlyMatch` (+21)
    - L290 PARM04 — `usingCartesianGrid`, `usingCylindricalGrid`, `usingSphericalPolarGrid`, `usingCurvilinearGrid`, `xgOrigin`, `ygOrigin`, `dxSpacing`, `dySpacing`, `delX`, `delY`, `delXFile`, `delYFile`, `horizGridFile`, `phiEuler`, `thetaEuler`, `psiEuler`, `rSphere`, `radius_fromHorizGrid`, `deepAtmosphere`, `seaLev_Z`, `top_Pres`, `delZ`, `delP`, `delR`, `delRc`, `delRFile`, `delRcFile`, `useMin4hFacEdges`, `interViscAr_pCell`, `interDiffKr_pCell` (+13)
    - L334 PARM05 — `bathyFile`, `topoFile`, `addWwallFile`, `addSwallFile`, `shelfIceFile`, `diffKrFile`, `viscAhDfile`, `viscAhZfile`, `viscA4Dfile`, `viscA4Zfile`, `hydrogThetaFile`, `hydrogSaltFile`, `maskIniTemp`, `maskIniSalt`, `checkIniTemp`, `checkIniSalt`, `zonalWindFile`, `meridWindFile`, `thetaClimFile`, `saltClimFile`, `surfQfile`, `surfQnetFile`, `surfQswFile`, `EmPmRfile`, `saltFluxFile`, `uVelInitFile`, `vVelInitFile`, `pSurfInitFile`, `dQdTFile`, `ploadFile` (+13)
  - L378 model/src/packages_boot.F  (reads data.pkg)
    - L379 PACKAGES — `useOBCS`, `useSHAP_FILT`, `useZONAL_FILT`, `useOPPS`, `usePP81`, `useKL10`, `useMY82`, `useGGL90`, `useKPP`, `useGMRedi`, `useDOWN_SLOPE`, `useBBL`, `useCAL`, `useEXF`, `useBulkForce`, `useEBM`, `useCheapAML`, `useAUTODIFF`, `useGrdchk`, `useSMOOTH`, `usePROFILES`, `useOBSFIT`, `useECCO`, `useCTRL`, `useSBO`, `useFLT`, `usePTRACERS`, `useGCHEM`, `useRBCS`, `useOffLine` (+22)
  - L432 eesupp/src/eeset_parms.F  (reads ?)
    - L433 EEPARMS — `usingMPI`, `useCubedSphereExchange`, `useCoupler`, `useNEST_PARENT`, `useNEST_CHILD`, `useNest2W_parent`, `useNest2W_child`, `useOASIS`, `useSETRLSTK`, `useSIGREG`, `debugMode`, `printMapIncludesZeros`, `maxLengthPrt1D`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/docs-ext.md`
- L1 External docs map (ECCO v4, ecco_darwin, darwin3 manual)
- L8 @ecco_v4_configs = `~/.cache/mitgcm-ecco-sources/ecco-docs/ECCO-v4-Configurations` (8393df4 2026-09-24 Merge pull request #246 from owang01/first_r7_commit)
  - L10 `README.md`
  - L29 `ECCOv4 Release 4/Docker/README.md`
  - L33 `ECCOv4 Release 4/flux-forced/doc/README_fluxforced.md`
  - L42 `ECCOv4 Release 4/flux-forced/doc/README_offline_ptracer.md` — `code_offline_ptracer`, `namelist_offline_ptracer`, `state_weekly`, `uVeltave.0000227808.data`, `uVeltave.0000000000.data`, `reverseintime_all.sh`, `state_weekly_rev_time_227808`, `state_weekly_rev_time_8904`
  - L54 `ECCOv4 Release 5/extension/code/slr_corr/README.md` — `slr_corr`, `data.slr_corr`, `SLR_CORR_PARAM.h`, `slrc_obs_timeseries`, `SLR_CORR_FIELDS.h`, `slr_corr_adjust_precip.F`, `slr_corr_init_varia.F`
  - L65 `devel/fluxforced/README.md`
  - L69 `devel/reproduction/README.md`
- L73 @ecco_v4_python = `~/.cache/mitgcm-ecco-sources/ecco-docs/ECCO-v4-Python-Tutorial` (3f0fcca 2026-05-15 Merge pull request #115 from andrewdelman/Steric_SSH_OBP)
  - L75 `README.md`
  - L78 `Cloud_Setup/JPL_setup_instructions.md` — `power_user`, `sshd_enable.sh`, `jupyter_env_setup.sh`
  - L84 `Docker/README.md` — `jupyter_env_setup.sh`, `ServerApp`, `LabApp`, `Tutorials_as_Jupyter_Notebooks`
  - L92 `Docker/smce/README.md` — `ServerApp`, `LabApp`, `Tutorials_as_Jupyter_Notebooks`
  - L100 `Intro_to_PO_Tutorials/Intro_to_PO_start.rst`
  - L107 `api_one_day/api.rst`
  - L121 `api_one_day/source/ecco_v4_py.rst`
  - L138 `api_one_day/source/modules.rst`
  - L141 `api_one_day/source/setup.rst`
  - L144 `doc/Installing_Python_and_Python_Packages.rst`
  - L157 `doc/Intro_to_PO_start.rst`
  - L163 `doc/Tutorial_Introduction.rst`
  - L177 `doc/Tutorial_wget_Command_Line_HTTPS_Downloading_ECCO_Datasets_from_PODAAC.md` — `urs_cookies`, `_netrc`
  - L184 `doc/fields.rst`
  - L195 `doc/index.rst`
  - L200 `doc/intro.rst`
  - L207 `doc/support.rst`
  - L212 `ecco_access/Downloading_ECCO_datasets_from_PODAAC/README.md`
  - L217 `ecco_access/Downloading_ECCO_datasets_from_PODAAC/Tutorial_wget_Command_Line_HTTPS_Downloading_ECCO_Datasets_from_PODAAC.md` — `urs_cookies`, `_netrc`
- L224 @eccov4_gael = `~/.cache/mitgcm-ecco-sources/ecco-docs/ECCOv4` (f2c08e6 2024-06-05 Merge pull request #21 from gaelforget/v1.12c)
  - L226 `README.md`
  - L230 `docs/ECCOv4r1_mods.md`
  - L232 `docs/ECCOv4r3_mods.md`
  - L239 `docs/analyses.rst`
  - L246 `docs/biblirefs.rst`
  - L248 `docs/downloads.rst`
  - L253 `docs/eccov4r2_dirtree.rst`
  - L255 `docs/eccov4r2_output.rst`
  - L257 `docs/eccov4r2_setup.rst`
  - L259 `docs/index.rst`
  - L262 `docs/introduction.rst`
  - L265 `docs/runs.rst` — `linux_amd64_gfortran`, `mitgcmuv_ad`
  - L272 `docs/example_scripts/README.md`
- L285 @ecco_darwin = `~/Documents/GitHub/ecco_darwin` (be01c37 2026-10-01 Merge branch 'master' of https://github.com/MITgcm-contrib/ecco_darwin)
  - L289 `adjoint/darwin_bling_comparison/README.md` — `run_sensitivity_sweep_fields.sh`, `compute_multipoint_comparison.py`, `ctrl_map_genarr3d.F`, `SpeciesAdj_CO2`, `adxx_ptr1`, `compute_cost_fd.py`, `DARWIN_CALC_PCO2_APPROX`, `CALC_PCO2_APPROX`, `bling_bio_nitrogen.F`, `global_oce_biogeo_bling_SOCAT`, `PO4_lim`, `darwin_vs_bling_multipoint.png`, `darwin_vs_bling_fd_vs_adjoint.png`, `initial_conditions_maps.png`, `initial_conditions_profiles.png`
  - L298 `adjoint/darwin_dic_comparison/README.md` — `darwin_bling_comparison`, `multipoint_comparison_results.csv`, `_minus`, `compute_multipoint_comparison_dic.py`, `run_dic_C`, `run_dic_D`, `DIC_OPTIONS.h`, `DIC_NO_NEG`, `ALLOW_FE`, `DIC_AD_SAFE`, `xx_ptr4`, `wt_DOP.bin`
  - L307 `adjoint/global_oce_biogeo_bling_SOCAT/README.md` — `global_oce_biogeo_bling`, `dic_init.bin`, `alk_init.bin`, `ones_32b.bin`, `sample_prof.nc`, `mah_flux_smooth.bin`, `socat_pco2_clim_month01..12.nc`, `SOCATv2026_tracks_gridded_monthly.nc`, `xx_theta`, `xx_salt`, `xx_ptr1`, `xx_ptr2`, `xx_ptr3`, `xx_ptr4`, `xx_ptr5`, `xx_ptr6`, `prof_PCOweight`, `mitgcmuv_ad`, `run_optim_loop.sh`, `m1qn3_output.txt`, `optim_readparms.F`, `ALLOW_USE_MPI`, `xx_genarr3d_weight`, `mult_genarr3d`, `prof_date`, `fco2_std_weighted`, `BLING_RUN_OPT_DIR`, `BATHY_BIN`, `FIGURES_OUT_DIR`
  - L321 `adjoint/global_oce_biogeo_darwin/README.md` — `global_oce_biogeo_bling_SOCAT`, `DARWIN_SIZE.h`, `darwin_check.F`, `EXTERNAL_FIELDS_LOAD`, `tauThetaClimRelax`, `tauSaltClimRelax`, `climsstTauRelax`, `climsssTauRelax`, `CAL_FULLDATE`, `mah_flux_smooth.bin`, `DIAGNOSTICS_SET_POINTERS`, `pCO2`
  - L330 `adjoint/hybrid_darwin_bling/README.md` — `mitgcmuv_ad`, `ctrl_map_genarr3d.F`, `wt_DIC.bin`, `wt_ALK.bin`, `run_bling`, `run_darwin`, `run_bling_1yr`, `run_darwin_1yr`, `mult_profiles`, `compute_darwin_cost.py`, `patch_ecco_cost.py`, `compute_darwin_cost_1yr.py`, `patch_ecco_cost_1yr.py`, `ctrl_map_ini_genarr.F`, `xx_genarr3d_file`, `darwinPCO2Snap01`, `prof_YYYYMMDD`, `prof_HHMMSS`, `pChkptFreq`, `chkptFreq`, `run_hybrid_loop_1yr.sh`, `nTimeSteps`, `m1qn3_output.txt`
  - L340 `adjoint/hybrid_darwin_dic/README.md` — `pilot_1month`, `ecco_cost_MIT_CE_000.optNNNN`, `xx_ptr1`, `adxx_ptr1`, `adxx_ptr2`, `dic_init.bin`, `alk_init.bin`, `dic_fc`, `adxx_salt`, `MDS_READ_FIELD`, `run_dic_A`, `run_dic_B`, `run_dic`, `dic_biotic_forcing.F`, `apply_control_to_darwin_ic_dic.py`, `xx_ptr2`
  - L353 `code_util/LOAC/C_GEM/NS_RAD/CLAUDE.md` — `make_report.py`, `build_all.sh`, `make_validation_pdf.py`, `singlechannel_archive`, `netCDF4`, `CGEM_SITE`, `CGEM_MAXT_DAYS`, `CGEM_WARMUP_DAYS`, `CGEM_OUTPUT`, `CGEM_TS`, `CGEM_ICE`, `CGEM_MULTICHANNEL`, `CGEM_N_CHAN_UP`, `CGEM_WATERTEMP_FILE`, `CGEM_DISTANCE`, `CGEM_FONT_SCALE`, `ns_rad`, `tridag_module`, `uphyd_module`, `hyd_module._hyd_iterate`, `disp_sch`, `fun_module`, `config.IS_IDEALIZED`, `schemes_module.openbound`, `WIND_FILE`, `SOLAR_FILE`, `AIRTEMP_FILE`, `RELHUM_FILE`, `PCO2_FILE`, `SEATEMP_FILE` (+140)
  - L382 `code_util/LOAC/C_GEM/NS_RAD/README.md` — `netCDF4`, `build_all.sh`, `_validation.pdf`, `ARCTIC_BGC`, `config.ARCTIC_BGC`
  - L393 `code_util/LOAC/C_GEM/NS_RAD/docs/FEATURES.md` — `CGEM_SITE`, `CGEM_MAXT_DAYS`, `CGEM_WARMUP_DAYS`, `CGEM_TS`, `CGEM_OUTPUT`, `CGEM_ICE`, `CGEM_MULTICHANNEL`, `CGEM_N_CHAN_UP`, `CGEM_DISTANCE`, `CGEM_FONT_SCALE`, `n_chan`, `multichannel_test.md`, `heat_module.py`, `ice_module.py`, `ice_model_plan.md`, `CARBONATE_UNITS`, `config.ARCTIC_BGC`, `arctic_biogeochemistry.md`, `lateral_module.py`, `LATERAL_INFLOW`, `LATERAL_CONC`, `BOUNDARY_FORCING`, `idealized_verification.md`, `fetch_discharge.py`, `build_tides.py`, `build_surge.py`, `build_river_temp.py`, `build_humidity.py`, `lter_boundary.py`, `build_boundary_chem.py` (+14)
  - L408 `code_util/LOAC/C_GEM/NS_RAD/docs/README.md` — `model_description.md`, `ice_model_plan.md`, `arctic_biogeochemistry.md`, `config.ARCTIC_BGC`, `idealized_verification.md`, `BOUNDARY_FORCING`, `multichannel_test.md`, `ns_rad_report.pdf`, `idealized_verification.pdf`, `ns_rad_interannual.pdf`, `idealized_forcings_preview.png`, `ns_rad_model_summary`, `ns_rad_model_schematic`, `ns_rad_geometry`, `ns_rad_river_networks`, `ns_rad_diagnostics`, `ns_rad_validation`, `make_report.py`, `performance_metrics.json`, `validation_obs.json`, `usgs_velocity_obs.json`, `sword_widths.json`, `interannual_discharge_obs.json`, `build_all.sh`
  - L415 `code_util/LOAC/C_GEM/NS_RAD/docs/arctic_biogeochemistry.md` — `LATERAL_INFLOW`, `rdoc_ox`, `ch4_ox`, `ch4_ex`, `n2o_prod`, `n2o_ex`, `f_photo_lab`, `I_use`, `variables.q_lat`, `LATERAL_CONC`, `BOUNDARY_FORCING`, `PHOTO_EFF`
  - L425 `code_util/LOAC/C_GEM/NS_RAD/docs/ice_model_plan.md` — `ice_module.py`, `heat_module.py`, `heat_module`, `k_ice_PAR`, `ICE_MODEL`, `ice_frac`, `ice_thickness`, `nbday_ice`, `file_module.exfread`, `transport_module`, `biogeo_module`, `sed_module`, `O2_ex`, `fetch_discharge.py`, `stefan_alpha`, `melt_rate`, `h_open`, `h_full`, `albedo_ice`, `fun_module.piston_velocity`, `ice_h`, `file_module.Rates`, `water_temp`, `fun_module`, `_fhet`, `_o2sat`, `piston_velocity`, `Q_sw`, `Q_sens`, `Q_lw` (+1)
  - L448 `code_util/LOAC/C_GEM/NS_RAD/docs/idealized_verification.md` — `config.IS_IDEALIZED`, `schemes_module.openbound`, `BOUNDARY_FORCING`, `IS_IDEALIZED`, `verify_idealized`
  - L457 `code_util/LOAC/C_GEM/NS_RAD/docs/model_description.md` — `hyd_module`, `heat_module`, `ice_module`, `BREAKUP_Q_FACTOR`, `config.ARCTIC_BGC`, `arctic_biogeochemistry.md`, `BOUNDARY_FORCING`
  - L491 `code_util/LOAC/C_GEM/NS_RAD/docs/multichannel_test.md` — `__pycache__`, `colville_minor.py`, `fun_module.river_dispersion`, `B_lb`, `B_ub`, `B_UB_TOTAL`, `n_chan`, `N_CHAN_LB`, `L_FLARE`, `B_thread`, `Q_FRACTION`, `n_chan_mod`, `S_max`, `fun_module._piston_velocity_loop`, `aer_deg`, `DISP_MAX`, `CGEM_MULTICHANNEL`
  - L505 `code_util/LOAC/C_GEM/NS_RAD/docs/performance.md` — `CGEM_MAXT_DAYS`, `CGEM_WARMUP_DAYS`, `cProfile`, `tridag_module.coeff_a`, `tridag_module.tridag`, `uphyd_module.update`, `tridag_module.conv`, `density._seck`, `biogeo_module.biogeo`, `density._dens0`, `fun_module.K1_CO2`, `fun_module.K2_CO2`, `fun_module.p_bar`, `K1_CO2`, `K2_CO2`, `tridag_module.py`, `uphyd_module.py`, `coeff_a`, `new_uh`, `_dens0`, `_seck`, `biogeo_module.py`, `__pycache__`, `fun_module.pH`, `disp_sch`, `sed_module.sed`, `fun.river_dispersion`, `fun.piston_velocity`, `transport_module`, `sed_module` (+11)
  - L528 `pkg/wad/README.md` — `addMass`, `WAD_PARM01`, `wad_check`, `vectorInvariantMomentum`, `upwindShear`, `wadMinDepth`, `wadCritDepth`, `wadAdvDepth`, `wadGMDepth`, `wadKPPCapDepth`, `wadKPPDepth`, `wadDryForcing`, `wadUpwindFace`, `wadCarryVel`, `momImplVertAdv`, `wadDragDepth`, `wadManningN`, `ALLOW_BOTTOMDRAG_ROUGHNESS`, `wadDragMax`, `wadManningFile`, `wadMaxFroude`, `wadMaxSpeed`, `wadConserveVol`, `wadMonFreq`, `wadSmoothWidth`, `WADdryC`, `pickup_wad`, `pickup_wadface`, `wad_thacker_1d`, `wad_balzano` (+3)
  - L540 `pkg/wad/doc/wad_formulation.md` — `CALC_R_STAR`, `hFacInf`, `Rmin_surf`, `hFacW`, `maskW`, `wadCritDepth`, `UPDATE_R_STAR`, `wadUpwindFace`, `useRealFreshWaterFlux`, `wadMonFreq`, `addMass`, `upwindShear`, `wadAdvDepth`, `wadCarryVel`, `wadKPPCapDepth`, `wadGMDepth`
  - L556 `pkg/wad/verification/wad_balzano/README.md` — `wadCritDepth`, `wadMonFreq`, `minDepth`, `wadMinDepth`
  - L564 `pkg/wad/verification/wad_estuary_3d/README.md` — `WAD_NX`, `WAD_NY`, `WAD_NR`, `wadMonFreq`, `minDepth`, `wadMinDepth`
  - L572 `pkg/wad/verification/wad_flat_xz/README.md` — `wadMonFreq`, `minDepth`, `wadMinDepth`
  - L580 `pkg/wad/verification/wad_mudflat/README.md` — `addMass`, `useRealFreshWaterFlux`, `addMass.bin`, `wadMonFreq`, `minDepth`, `wadMinDepth`
  - L588 `pkg/wad/verification/wad_thacker_1d/README.md` — `bathy_thacker.bin`, `eta_thacker.bin`, `wadMonFreq`, `minDepth`, `wadMinDepth`, `darwin_arm64_gfortran`
  - L604 `regions/GoM/GoM_1km/README.md` — `darwin_arm64_gfortran`
  - L652 `regions/downscaling/README.md` — `diagnostic_vec`
  - L660 `regions/downscaling/STEP1.md` — `diagnostics_vec`, `gen_mitgrid.py`, `gen_bathy.py`, `regnm_bathymetry.bin`, `stitch_ncgrid.py`, `job_ECCO_darwin`, `gen_dvmasks.py`
  - L674 `regions/downscaling/STEP2.md` — `diagnostic_vec`, `data.diagnostics_vec`, `DIAGNOSTICS_VEC_SIZE.h`, `gen_dvmasks`
  - L683 `regions/downscaling/STEP3.md` — `gen_pickups.py`, `gen_obcs.py`, `plot_boundary_velocity_profile.py`, `gen_obcs2.py`, `check_obcs_transport.py`, `diagnostics_vec`, `_ncgrid.nc`, `HFacS`, `gen_mod`
  - L692 `regions/downscaling/STEP4.md`
  - L712 `regions/kerguelen/README.md` — `kerguelen_ncgrid.nc`, `kerguelen_bathymetry.bin`, `diagnostics_vec`, `code_darwin`, `OBCS_OPTIONS.h`, `DIAGNOSTICS_SIZE.h`, `ALLOW_OBCS_SPONGE`
  - L719 `regions/kerguelen/v05/PIPELINE.md` — `run_step1_mitgrid.sh`, `job_step1_bathy.sh`, `kerguelen_bathymetry.bin`, `job_step1_stitch.sh`, `kerguelen_ncgrid.nc`, `job_step1_dvmasks.sh`, `delX`, `delY`, `delYFile`, `nVEC_mask`, `nml_vecFiles`, `VEC_points`, `job_gen_obcs_fast.sh`, `gen_obcs_fast.py`, `THETA_east.bin`, `job_gen_obcs_transport_correction.sh`, `nonlinFreeSurf`, `exactConserv`, `job_gen_pickups.sh`, `gen_pickups.py`, `hydrogThetaFile`, `hydrogSaltFile`, `readme_pleiades.txt`, `code_darwin`
  - L731 `regions/kerguelen/v05/scripts/pipeline/README.md` — `run_step1_mitgrid.sh`, `run_step1_mitgrid_repad.sh`, `job_step1_bathy.sh`, `kerguelen_bathymetry.bin`, `job_step1_stitch.sh`, `kerguelen_ncgrid.nc`, `job_step1_dvmasks.sh`, `job_step2_parent_dv.sh`, `diagnostics_vec`, `job_gen_obcs_fast.sh`, `job_gen_obcs_transport_correction.sh`, `job_gen_pickups.sh`
  - L736 `regions/mac_delta/llc270/biogeochem_setup/README.md`
  - L745 `regions/mac_delta/llc270/biogeochem_setup/carroll_2020_ecosystem/CDOM_setup/README.md`
  - L751 `regions/mac_delta/llc270/biogeochem_setup/carroll_2020_ecosystem/ED-SBS_Bertin_etal_2023/README.md`
  - L759 `regions/mac_delta/llc270/biogeochem_setup/carroll_2020_ecosystem/ED-SBS_Bertin_etal_2023/setup_files/README.md`
  - L765 `regions/mac_delta/llc270/biogeochem_setup/carroll_2020_ecosystem/baseline_setup/README.md` — `PATH_TO_MPI_ENVIRONMENT_VARIABLE`, `h_mpi`, `data_mpi`, `darwin_amd64_gfortran`
- L824 @darwin3_manual = `~/Documents/GitHub/darwin3/doc` (28ac947a3 2024-09-24 Revert scavenging defaults to POP-based equivalents)
  - L826 `phys_pkgs/darwin.rst` — `GCHEM_SEPARATE_FORCING`, `ALLOW_CLIMSST_RELAXATION`, `ALLOW_CLIMSSS_RELAXATION`, `DARWIN_ALLOW_NQUOTA`, `DARWIN_ALLOW_PQUOTA`, `DARWIN_ALLOW_FEQUOTA`, `DARWIN_ALLOW_SIQUOTA`, `DARWIN_ALLOW_CHLQUOTA`, `DARWIN_ALLOW_CDOM`, `DARWIN_CDOM_UNITS_CARBON`, `DARWIN_ALLOW_CSTORE`, `DARWIN_ALLOW_CSTORE_DIAGS`, `DARWIN_ALLOW_CARBON`, `DARWIN_SOLVESAPHE`, `DARWIN_TOTALPHSCALE`, `DARWIN_USE_PLOAD`, `DARWIN_ALLOW_RADI`, `DARWIN_ALLOW_DENIT`, `DARWIN_ALLOW_EXUDE`, `ALLOW_OLD_VIRTUALFLUX`, `DARWIN_NITRATE_FELIMIT`, `DARWIN_BOTTOM_SINK`, `DARWIN_NUTRIENT_RUNOFF`, `DARWIN_AVPAR`, `DARWIN_ALLOW_GEIDER`, `DARWIN_GEIDER_RHO_SYNTH`, `DARWIN_CHL_INIT_LEGACY`, `DARWIN_SCATTER_CHL`, `DARWIN_DIAG_IOP`, `DARWIN_GRAZING_SWITCH` (+69)
  - L837 `phys_pkgs/darwin_airsea.rst` — `DARWIN_ALLOW_CARBON`, `useSEAICE`, `darwin_useSEAICE`
  - L840 `phys_pkgs/darwin_bacteria.rst` — `bactType`, `isAerobic`, `isDenit`, `grp_bacttype`, `grp_aerobic`, `grp_denit`, `pcoefO2`, `pmaxDIN`, `ksatDIN`, `alpha_hydrol`, `b_PCmax`, `yieldO2`, `yieldNO3`, `ksatPON`, `a_ksatPON`, `ksatPOC`, `ksatPOP`, `ksatPOFe`, `ksatDON`, `a_ksatDON`, `ksatDOC`, `ksatDOP`, `ksatDOFe`
  - L847 `phys_pkgs/darwin_carbon.rst` — `DARWIN_ALLOW_CARBON`, `DARWIN_SOLVESAPHE`, `DARWIN_TOTALPHSCALE`, `DARWIN_ALLOW_RADI`, `selectPHsolver`, `surfSaltMinInit`, `surfSaltMaxInit`, `surfTempMinInit`, `surfTempMaxInit`, `surfDICMinInit`, `surfDICMaxInit`, `surfALKMinInit`, `surfALKMaxInit`, `surfPO4MinInit`, `surfPO4MaxInit`, `surfSiMinInit`, `surfSiMaxInit`, `surfSaltMin`, `surfSaltMax`, `surfTempMin`, `surfTempMax`, `surfDICMin`, `surfDICMax`, `surfALKMin`, `surfALKMax`, `surfPO4Min`, `surfPO4Max`, `surfSiMin`, `surfSiMax`, `darwin_disscSelect` (+9)
  - L853 `phys_pkgs/darwin_cdom.rst` — `DARWIN_ALLOW_CDOM`, `DARWIN_CDOM_UNITS_CARBON`, `DARWIN_ALLOW_DENIT`, `fracCDOM`, `R_NP_CDOM`, `R_FeP_CDOM`, `R_CP_CDOM`, `R_NC_CDOM`, `R_PC_CDOM`, `R_FeC_CDOM`
  - L856 `phys_pkgs/darwin_changes.rst`
  - L859 `phys_pkgs/darwin_chl.rst` — `DARWIN_GEIDER_RHO_SYNTH`, `a_acclimtimescl`
  - L865 `phys_pkgs/darwin_cons.rst` — `DARWIN_ALLOW_CONS`, `DARWIN_BOTTOM_SINK`, `darwin_linFSConserv`, `DARWIN_ALLOW_CARBON`, `ALLOW_OLD_VIRTUALFLUX`, `ironFile`, `DARWIN_MINFE`, `DARWIN_PART_SCAV`, `fesedflux_pcm`, `DARWIN_IRON_SED_SOURCE_VARIABLE`, `DARWIN_ALLOW_DENIT`, `DARWIN_ALLOW_NQUOTA`, `implicitFreeSurface`, `DARWIN_linFSConserve`
  - L869 `phys_pkgs/darwin_cstore.rst` — `DARWIN_ALLOW_NQUOTA`, `DARWIN_ALLOW_PQUOTA`, `DARWIN_ALLOW_FEQUOTA`, `DARWIN_ALLOW_SIQUOTA`, `DARWIN_ALLOW_CSTORE`, `FracExudeC`, `a_FracExudeC`
  - L872 `phys_pkgs/darwin_denit.rst` — `DARWIN_ALLOW_DENIT`, `DARWIN_ALLOW_CDOM`, `denit_NP`, `denit_NO3`
  - L875 `phys_pkgs/darwin_equations.rst` — `DARWIN_ALLOW_CDOM`, `DARWIN_ALLOW_CARBON`, `R_PICPOC`, `a_R_PICPOC`, `R_OP`
  - L878 `phys_pkgs/darwin_exude.rst` — `DARWIN_ALLOW_EXUDE`, `a_kexcC`, `a_kexcN`, `a_kexcP`, `a_kexcFe`, `a_kexcSi`, `ExportFracExude`, `a_ExportFracExude`
  - L881 `phys_pkgs/darwin_grazing.rst` — `DARWIN_GRAZING_SWITCH`, `hillnumGraz`, `a_grazemax`, `a_kgrazesat`, `tempGraz`, `grp_tempGraz`, `inhib_graz`, `inhib_graz_exp`, `grp_ass_eff`, `ExportFracPreyPred`, `grp_ExportFracPreyPred`, `DARWIN_ALLOMETRIC_PALAT`, `grp_pred`, `grp_prey`, `b_ppOpt`, `a_ppSig`, `palat_min`
  - L886 `phys_pkgs/darwin_growth.rst` — `DARWIN_ALLOW_GEIDER`, `aphy_chl_ave`, `DARWIN_ALLOW_CHLQUOTA`, `ksatPAR`, `a_ksatPAR`, `kinhPAR`, `a_kinhPAR`, `b_PCmax`, `mQyield`, `a_mQyield`, `a_chl2cmax`, `inhibGeider`, `a_inhibGeider`, `aphy_chl_ps`, `aphy_chl_ps_type`, `grp_aptype`, `darwin_phytoAbsorbFile`
  - L889 `phys_pkgs/darwin_iron.rst` — `darwin_inscal_iron`, `DARWIN_IRON_SED_SOURCE_VARIABLE`, `DARWIN_IRON_SED_SOURCE_POP`, `DARWIN_ALLOW_HYDROTHERMAL_VENTS`, `depthFeVent`, `ventHe3file`, `solFeVent`, `R_FeHe3_vent`, `DARWIN_PART_SCAV`, `DARWIN_PART_SCAV_POP`, `DARWIN_MINFE`, `fesedflux_pcm`, `fesedflux_min`, `R_CP_fesed`, `scav_tau`, `scav_inter`, `scav_exp`, `scav_POC_wgt`, `scav_PSi_wgt`, `scav_PIC_wgt`, `scav_rPOM`, `ligand_tot`, `ligand_stab`, `scav_rat`, `scav_R_POPPOC`
  - L896 `phys_pkgs/darwin_light.rst` — `darwin_useQsw`, `DARWIN_AVPAR`, `DARWIN_ALLOW_CHLQUOTA`, `DARWIN_ALLOW_GEIDER`, `R_ChlC`, `katten_w`, `katten_Chl`, `aphy_chl_ave`
  - L899 `phys_pkgs/darwin_mort.rst` — `a_mort`, `a_mort2`, `a_Xmin`, `tempMort`, `grp_tempMort`, `tempMort2`, `grp_tempMort2`, `ExportFracMort`, `a_ExportFracMort`, `ExportFracMort2`, `a_ExportFracMort2`
  - L903 `phys_pkgs/darwin_remin.rst` — `DARWIN_ALLOW_DENIT`, `PAR_oxi`
  - L906 `phys_pkgs/darwin_resp.rst` — `respRate`, `b_respRate_c`, `b_qcarbon`, `a_Xmin`, `a_respRate_c`
  - L910 `phys_pkgs/darwin_sink.rst` — `DARWIN_BOTTOM_SINK`, `wPIC_sink`, `wC_sink`, `wN_sink`, `wP_sink`, `wSi_sink`, `wFe_sink`, `b_biosink`, `b_bioswim`
  - L913 `phys_pkgs/darwin_spectral.rst` — `darwin_waterAbsorbFile`, `grp_aptype`, `darwin_phytoAbsorbFile`, `DARWIN_SCATTER_CHL`, `darwin_particleAbsorbFile`, `DARWIN_ALLOW_CDOM`, `darwin_bbmin`, `darwin_bbw`, `darwin_RPOC`, `darwin_rCDOM`, `darwin_lambda_aCDOM`, `darwin_Sdom`, `darwin_aCDOM_fac`, `darwin_part_size_P`, `aphy_chl`, `aphy_chl_ps`, `aphy_mgC`, `bphy_mgC`, `bbphy_mgC`, `darwin_allomSpectra`, `darwin_aCarCell`, `darwin_bCarCell`, `darwin_absorpSlope`, `darwin_bbbSlope`, `darwin_scatSwitchSizeLog`, `darwin_scatSlopeSmall`, `darwin_scatSlopeLarge`
  - L922 `phys_pkgs/darwin_temperature.rst` — `DARWIN_TEMP_VERSION`, `DARWIN_TEMP_RANGE`, `DARWIN_NOTEMP`, `tempMort`, `tempMort2`, `tempGraz`, `phytoTempCoeff`, `a_phytoTempCoeff`, `phytoTempExp1`, `a_phytoTempExp1`, `TempCoeffArr`, `TempAeArr`, `TempRefArr`, `phytoTempAe`, `a_phytoTempAe`, `hetTempAe`, `a_hetTempAe`, `grazTempAe`, `a_grazTempAe`, `reminTempAe`, `mortTempAe`, `mort2TempAe`, `uptakeTempAe`, `phytoTempExp2`, `a_phytoTempExp2`, `phytoTempOptimum`, `a_phytoTempOptimum`, `phytoDecayPower`, `a_phytoDecayPower`, `hetTempExp2` (+4)
  - L929 `phys_pkgs/darwin_uptake.rst` — `hasSi`, `combNO`, `DARWIN_effective_ksat`, `darwin_select_kn_allom`, `a_ksatNO3`, `b_ksatNO3`, `grp_hasSi`, `grp_diazo`, `useNH4`, `grp_useNH4`, `useNO2`, `grp_useNO2`, `useNO3`, `grp_useNO3`, `grp_combNO`, `b_Qnmin`, `b_Qnmax`, `b_Qpmin`, `b_Qpmax`, `b_Qsimin`, `b_Qsimax`, `b_Qfemin`, `b_Qfemax`
  - L945 `phys_pkgs/radtrans.rst` — `RT_useMeanCosSolz`, `ALLOW_CLIMSST_RELAXATION`, `ALLOW_CLIMSSS_RELAXATION`, `RADTRANS_DIAG_SOLUTION`, `useRADTRANS`, `RADTRANS_FORCING_PARAMS`, `RADTRANS_INTERP_PARAMS`, `USE_EXF_INTERPOLATION`, `RADTRANS_PARAMS`, `RT_Edfile`, `RT_Esfile`, `RT_Ed_const`, `RT_Es_const`, `RT_E_period`, `RT_E_RepCycle`, `repeatCycle`, `useOasimYearlyFields`, `RT_E_StartTime`, `RT_E_startdate1`, `RT_E_startdate2`, `RT_E_mask`, `RT_Ed_exfremo_intercept`, `RT_Es_exfremo_intercept`, `RT_Ed_exfremo_slope`, `RT_Es_exfremo_slope`, `RT_inscal_Ed`, `RT_inscal_Es`, `RT_refract_water`, `RT_rmud_max`, `RT_wbEdges` (+8)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/docs.md`
- L1 MITgcm manual section map (doc/**/*.rst)
  - L9 `doc/algorithm/adv-schemes.rst` — Linear advection schemes — `uTrans`, `vTrans`, `rTrans`, `gad_c2_adv_x.F`, `gad_c2_adv_y.F`, `gad_c2_adv_r.F`, `gad_u3_adv_x.F`, `gad_u3_adv_y.F`, `gad_u3_adv_r.F`, `gad_c4_adv_x.F`, `gad_c4_adv_y.F`, `gad_c4_adv_r.F`, `gad_fluxlimit_adv_x.F`, `gad_fluxlimit_adv_y.F`, `gad_fluxlimit_adv_r.F`, `gad_dst3_adv_x.F`, `gad_dst3_adv_y.F`, `gad_dst3_adv_r.F`, `gad_dst3fl_adv_x.F`, `gad_dst3fl_adv_y.F`, `gad_dst3fl_adv_r.F`, `gTracer`, `gad_advection.F`, `multiDimAdvection`
  - L22 `doc/algorithm/algorithm.rst` — Discretization and Algorithm — `implicitViscosity`, `timestep.F`, `solve_for_pressure.F`, `correction_step.F`, `forward_step.F`, `dynamics.F`, `calc_div_ghat.F`, `cg2d.F`, `momentum_correction_step.F`, `calc_grad_phi_surf.F`, `freesurfFac`, `forcing_In_AB`, `tracForcingOutAB`, `momForcingOutAB`, `momDissip_In_AB`, `ALLOW_ADAMSBASHFORTH_3`, `CPP_OPTIONS.h`, `thermodynamics.F`, `temp_integrate.F`, `gad_calc_rhs.F`, `external_forcing.F`, `adams_bashforth2.F`, `timestep_tracer.F`, `impldiff.F`, `hFac`, `external_fields_load.F`, `do_atmospheric_phys.F`, `do_oceanic_phys.F`, `calc_gt.F`, `calc_phi_hyd.F` (+95)
  - L75 `doc/algorithm/c-grid.rst` — C grid staggering of variables — `ini_grid.F`, `ini_vertical_grid.F`, `ini_cartesian_grid.F`, `ini_spherical_polar_grid.F`, `ini_curvilinear_grid.F`, `ini_masks_etc.F`, `initialise_fixed.F`, `GRID.h`
  - L79 `doc/algorithm/crank-nicol.rst` — Crank-Nicolson barotropic time stepping — `implicSurfPress`, `implicDiv2DFlow`, `useRealFreshWaterFlux`, `implicitNHPress`, `correction_step.F`
  - L82 `doc/algorithm/finitevol-meth.rst` — The finite volume method: finite volumes versus finite difference
  - L85 `doc/algorithm/horiz-grid.rst` — Horizontal grid — `ini_cartesian_grid.F`, `ini_spherical_polar_grid.F`, `ini_curvilinear_grid.F`, `GRID.h`, `recip_rA`, `recip_rAz`, `recip_rAw`, `recip_rAs`, `recip_dxG`, `recip_dyG`, `recip_dxC`, `recip_dyC`, `recip_dxF`, `recip_dyF`, `recip_dxV`, `recip_dyU`, `recip_`, `ini_masks_etc.F`, `usingCartesianGrid`, `dXspacing`, `dYspacing`, `usingSphericalPolarGrid`, `usingCurvilinearGrid`
  - L92 `doc/algorithm/nonlinear-freesurf.rst` — Non-linear free-surface — `nonlinFreeSurf`, `uniformLin_PhiSurf`, `ini_linear_phisurf.F`, `EXACT_CONSERV`, `select_rStar`, `exactConserv`, `cg2dTargetResidual`, `etaH`, `etaHnm1`, `dEtaHdt`, `solve_for_pressure.F`, `integr_continuity.F`, `timestep.F`, `calc_gt.F`, `calc_gs.F`, `forward_step.F`, `DYNVARS.h`, `SURFACE.h`, `hFacInf`
  - L100 `doc/algorithm/vert-grid.rst` — Vertical grid — `delR`, `delZ`, `delP`, `delRc`, `recip_drF`, `recip_drC`, `ini_vertical_grid.F`, `GRID.h`, `hFacC`, `hFacW`, `hFacS`, `recip_hFacC`, `recip_hFacW`, `recip_hFacS`, `hFacMin`, `hFacMinDr`, `ini_masks_etc.F`
  - L104 `doc/autodiff/autodiff.rst` — Automatic Differentiation — `ALLOW_AUTODIFF_TAMC`, `ALLOW_ADJOINT_RUN`, `ALLOW_TANGENTLINEAR_RUN`, `ALLOW_GRDCHK`, `the_model_main.F`, `the_main_loop.F`, `ctrl_unpack.F`, `ctrl_pack.F`, `grdchk_main.F`, `AUTODIFF_OPTIONS.h`, `COST_OPTIONS.h`, `CTRL_OPTIONS.h`, `tamc.h`, `myThid`, `AD_FILES`, `AD_FLOW_FILES`, `ALLOW_COST`, `ALLOW_COST_TRACER`, `mult_tracer`, `objf_tracer`, `cost_readparms.F`, `cost_init_varia.F`, `cost_tile.F`, `cost_tracer.F`, `cost_final.F`, `nWetCtile`, `nWetWtile`, `nWetStile`, `ALLOW_NONDIMENSIONAL_CONTROL_IO`, `ctrl_readparms.F` (+34)
  - L155 `doc/contributing/contributing.rst` — Contributing to the MITgcm — `CPP_OPTIONS.h`, `GMREDI_OPTIONS.h`, `packages_boot.F`, `packages_readparms.F`, `initialise_fixed.F`, `packages_init_fixed.F`, `packages_check.F`, `packages_init_variables.F`, `initialise_varia.F`, `gmredi_calc_tensor.F`, `do_oceanic_phys.F`, `do_the_model_io.F`, `packages_write_pickup.F`, `the_model_main.F`, `gmredi_diagnostics_init.F`, `_mpi`, `mitgcmuv_ad`, `genmake_local`, `prepare_run`, `SIZE.h`, `cg2d_init_res`, `tr_checklist`, `global_sum.F`, `dynamics.F`, `subroutine_name.F`, `where_var1_defined.h`, `where_var2_defined.h`, `where_var3_defined.h`
  - L190 `doc/examples/advection_in_gyre/advection_in_gyre.rst` — Ocean Gyre Advection Schemes — `PTRACERS_ALLOW_DYN_STATE`, `PTRACERS_OPTIONS.h`
  - L197 `doc/examples/baroclinic_gyre/baroclinic_gyre.rst` — Baroclinic Ocean Gyre — `gBaro`, `SIZE.h`, `DIAGNOSTICS_SIZE.h`, `numDiags`, `useMNC`, `useDiagnostics`, `monitor_mnc`, `mnc_use_outdir`, `mnc_outdir_str`, `mnc_test_`, `endTime`, `monitorFreq`, `HFacC`, `HFacW`, `HFacS`, `nIter0`, `dumpFreq`, `deltaT`, `mnc_test_0001`, `mnc_test_0002`, `DIAG_STATIS_PARMS`, `PHrefC`, `PHrefF`, `RhoRef`, `phiHyd`, `build_mpi`, `mnc_test_00001`, `mnc_test_00004`, `globalFiles`, `build_openmp` (+6)
  - L221 `doc/examples/barotropic_gyre/barotropic_gyre.rst` — Barotropic Ocean Gyre — `SIZE.h`, `monitorFreq`, `debugLevel`, `etaN`, `dumpFreq`, `uVel`, `wVel`, `pChkptFreq`, `ChkptFreq`, `hFacC`, `hFacS`, `hFacW`, `PHrefC`, `PHrefF`
  - L240 `doc/examples/cfc_offline/cfc_offline.rst` — Offline Experiments — `PTRACERS_SIZE.h`, `GMREDI_OPTIONS.h`, `SIZE.h`, `nIter0`, `nTimesteps`, `deltaTtracer`, `deltatTracer`, `deltaToffline`, `deltaTClock`, `deltatTtracer`, `periodicExternalForcing`, `externForcingPeriod`, `externForcingCycle`, `UvelFile`, `VvelFile`, `WvelFile`, `ConvFile`, `offlineForcingPeriod`, `offlineForcingCycle`, `offlineIter0`, `uVeltave`, `vVeltave`, `wVeltave`, `deltatToffline`, `usePTRACERS`, `PTRACERS_numInUse`, `PTRACERS_Iter0`, `PTRACERS_advScheme`, `PTRACERS_diffKh`, `PTRACERS_diffKr` (+10)
  - L264 `doc/examples/deep_convection/deep_convection.rst` — Deep Convection — `SIZE.h`, `CPP_OPTIONS.h`
  - L278 `doc/examples/examples.rst` — MITgcm Tutorial Example Experiments — `SIZE.h`, `CPP_OPTIONS.h`, `useLANGMUIR`, `code_ad`, `input_ad`, `code_tap`, `input_tap`, `cg2d_nsa.F`, `cg2d.F`, `cg2d_mad.F`
  - L283 `doc/examples/global_oce_biogeo/global_oce_biogeo.rst` — Biogeochemistry Simulation — `nTimeSteps`, `taveFreq`, `PTRACERS_numInUse`, `PTRACERS_Iter0`, `SIZE.h`, `PTRACERS_SIZE.h`, `DIAGNOSTICS_SIZE.h`, `DIC_Biotave`, `DIC_Cartave`, `DIC_fluxCO2ave`, `DIC_pCO2tave`, `DIC_pHtave`, `DIC_SurOtave`, `DIC_Surtave`
  - L290 `doc/examples/global_oce_in_p/global_oce_in_p.rst` — Global Ocean Simulation in Pressure Coordinates — `CPP_OPTIONS.h`, `SIZE.h`, `startTime`, `ATMOSPHERIC_LOADING`, `EXACT_CONSERV`, `NONLIN_FRSURF`
  - L308 `doc/examples/global_oce_latlon/global_oce_latlon.rst` — Global Ocean Simulation — `SIZE.h`, `nIter0`, `startTime`, `endTime`
  - L322 `doc/examples/global_oce_optim/global_oce_optim.rst` — Global Ocean State Estimation — `ALLOW_HFLUXM_CONTROL`, `cHFLUXM_CONTROL`, `CTRL_OPTIONS.h`, `FFIELDS.h`, `ini_forcing.F`, `external_forcing_surf.F`, `ctrl_init.F`, `ctrl_pack.F`, `ctrl_unpack.F`, `ctrl_map_forcing.F`, `ALLOW_COST_TEMP`, `objf_temp`, `mult_temp`, `ALLOW_COST_HFLUXM`, `objf_hfluxm`, `mult_hflux`, `cost_temp.F`, `COST_OPTIONS.h`, `cost_tile.F`, `cost_hflux.F`, `cost_final.F`, `cost_weights.F`, `cost_readparms.F`, `ALLOW_ECCO_OPTIMIZATION`, `cg2dTargetResidual`, `useGrdchk`, `mitgcmuv_ad`, `mitcgmuv_ad`, `xx_hfluxm_file`, `mitgcm_ad` (+2)
  - L336 `doc/examples/held_suarez_cs/held_suarez_cs.rst` — Held-Suarez Atmosphere — `CPP_OPTIONS.h`, `SIZE.h`, `DIAGNOSTICS_SIZE.h`, `apply_forcing.F`, `nIter0`, `startTime`, `grid_cs32`, `nShapUV`, `nShapUVPhys`, `useEXCH2`, `EXTERNAL_FORCING_U`, `EXTERNAL_FORCING_V`, `EXTERNAL_FORCING_T`, `EXTERNAL_FORCING_S`, `GRID.h`
  - L352 `doc/examples/plume_on_slope/plume_on_slope.rst` — Gravity Plume On a Continental Slope — `ALLOW_NONHYDROSTATIC`, `nonHydrostatic`, `useOBCS`, `SIZE.h`, `CPP_OPTIONS.h`, `eosType`, `sBeta`
  - L359 `doc/examples/reentrant_channel/reentrant_channel.rst` — Southern Ocean Reentrant Channel Example — `tempAdvScheme`, `implicitViscosity`, `SIZE.h`, `LAYERS_SIZE.h`, `DIAGNOSTICS_SIZE.h`, `diffKh`, `numDiags`, `background_K`, `GM_isopycK`, `useRBCtemp`, `tauRelaxT`, `relaxMaskFile`, `relaxTFile`, `layers_maxNum`, `layers_name`, `layers_bounds`, `hFacMin`, `hFacMinDr`, `hFacC`, `nTimeSteps`, `monitorFreq`, `useGMRedi`, `diffKhT`, `stdiags_bylev`, `stdiags_2D`, `DeltaT`, `viscAh`, `viscC2Leith`, `useFullLeith`, `viscAhGridMax` (+3)
  - L387 `doc/examples/rotating_tank/rotating_tank.rst` — Rotating Tank — `CPP_OPTIONS.h`, `SIZE.h`
  - L400 `doc/examples/tracer_adjsens/tracer_adjsens.rst` — Adjoint Sensitivity Analysis for Tracer Injection — `COST_OPTIONS.h`, `CTRL_OPTIONS.h`, `CPP_OPTIONS.h`, `AUTODIFF_OPTIONS.h`, `CTRL_SIZE.h`, `GAD_OPTIONS.h`, `GMREDI_OPTIONS.h`, `SIZE.h`, `tamc.h`, `ctrl_map_ini_genarr.F`, `ptracers_forcing_surf.F`, `the_main_loop.F`, `ECCO_OPTIONS.h`, `ALLOW_AUTODIFF_TAMC`, `ALLOW_TAMC_CHECKPOINTING`, `ALLOW_AUTODIFF_MONITOR`, `ALLOW_DIVIDED_ADJOINT`, `ALLOW_COST`, `ALLOW_COST_TRACER`, `ALLOW_THETA0_CONTROL`, `ALLOW_SALT0_CONTROL`, `ALLOW_TR10_CONTROL`, `ALLOW_TAUU0_CONTROL`, `ALLOW_TAUV0_CONTROL`, `ALLOW_SFLUX0_CONTROL`, `ALLOW_HFLUX0_CONTROL`, `ALLOW_DIFFKR_CONTROL`, `ALLOW_KAPGM_CONTROL`, `checkpoint_lev3_directives.h`, `checkpoint_lev2_directives.h` (+18)
  - L421 `doc/getting_started/getting_started.rst` — Getting Started with MITgcm — `main.F`, `the_model_main.F`, `SIZE.h`, `genmake_local`, `EEPARAMS.h`, `prepare_run`, `OMP_NUM_THREADS`, `OMP_STACKSIZE`, `ckptA`, `ckptB`, `mnc_output_`, `SHORTWAVE_HEATING`, `ALLOW_GEOTHERMAL_FLUX`, `ALLOW_FRICTION_HEATING`, `ALLOW_ADDFLUID`, `ATMOSPHERIC_LOADING`, `ALLOW_BALANCE_FLUXES`, `ALLOW_BALANCE_RELAX`, `CHECK_SALINITY_FOR_NEGATIVE_VALUES`, `EXCLUDE_FFIELDS_LOAD`, `INCLUDE_PHIHYD_CALCULATION_CODE`, `INCLUDE_CONVECT_CALL`, `INCLUDE_CALC_DIFFUSIVITY_CALL`, `ALLOW_3D_DIFFKR`, `ALLOW_BL79_LAT_VARY`, `EXCLUDE_PCELL_MIX_CODE`, `ALLOW_SMAG_3D_DIFFUSIVITY`, `ALLOW_SOLVE4_PS_AND_DRAG`, `INCLUDE_IMPLVERTADV_CODE`, `ALLOW_ADAMSBASHFORTH_3` (+341)
  - L490 `doc/index.rst` — Welcome to MITgcm's user manual
  - L493 `doc/ocean_state_est/ocean_state_est.rst` — Packages III - Ocean State Estimation — `gencost_datafile`, `gencost_barfile`, `gencost_outputlevel`, `gencost_name`, `gencost_errfile`, `gencost_avgperiod`, `gencost_preproc`, `gencost_preproc_i`, `gencost_posproc`, `gencost_posproc_c`, `gencost_posproc_i`, `gencost_kLev_select`, `gencost_is3d`, `gencost_mask`, `mult_gencost`, `gencost_preproc_c`, `gencost_preproc_r`, `gencost_posproc_r`, `gencost_spmin`, `gencost_spmax`, `gencost_spzero`, `gencost_startdate1`, `gencost_startdate2`, `gencost_enddate1`, `gencost_enddate2`, `gencost_useDensityMask`, `gencost_sigmaLow`, `gencost_sigmaHigh`, `gencost_refPressure`, `gencost_tanhScale` (+147)
  - L538 `doc/outp_pkgs/flt.rst` — Introduction — `FLT_OPTIONS.h`, `ALLOW_3D_FLT`, `USE_FLT_ALT_NOISE`, `ALLOW_FLT_3D_NOISE`, `FLT_SECOND_ORDER_RUNGE_KUTTA`, `FLT_WITHOUT_X_PERIODICITY`, `FLT_WITHOUT_Y_PERIODICITY`, `DEVEL_FLT_EXCH2`, `FLT_SIZE.h`
  - L548 `doc/outp_pkgs/outp_pkgs.rst` — Packages II - Diagnostics and I/O — `diagnostics_addtolist.F`, `diagnostics_fill.F`, `diagnostics_scale_fill.F`, `diagnostics_fract_fill.F`, `diagnostics_frac_fill.F`, `diagnostics_is_on.F`, `diagnostics_count.F`, `DIAGNOSTICS.h`, `useDiagnostics`, `fileFlags`, `writeBinaryPrec`, `numDiags`, `numLists`, `numperList`, `diagSt_size`, `DIAGNOSTICS_SIZE.h`, `diagMate`, `diagNum`, `diagCode`, `diagnostics_main_init.F`, `precFloat64`, `precFloat32`, `globalFile`, `SAFE_IO`, `ALLOW_WHIO`, `_BYTESWAPIO`, `RW_SAFE_MFLDS`, `RW_DISABLE_SMALL_OVERLAP`, `MNC_COMMON.h`, `useMNC` (+28)
  - L609 `doc/overview/adjoint.rst` — Adjoint
  - L612 `doc/overview/atmosphere.rst` — Atmosphere
  - L615 `doc/overview/bound_forc_inter_waves.rst` — Boundary forced internal waves
  - L618 `doc/overview/coordinate_sys.rst` — Coordinate systems
  - L622 `doc/overview/cvct_mixing_topo.rst` — Convection and mixing over topography
  - L625 `doc/overview/eqn_motion_ocn.rst` — Equations of Motion for the Ocean
  - L632 `doc/overview/finding_pressure.rst` — Finding the pressure field
  - L639 `doc/overview/forcing_dissip.rst` — Forcing/dissipation
  - L646 `doc/overview/global_atmos_hs.rst` — Global atmosphere: ‘Held-Suarez’ benchmark
  - L649 `doc/overview/global_ocean_circ.rst` — Global ocean circulation
  - L652 `doc/overview/global_state_est.rst` — Global state estimation of the ocean
  - L655 `doc/overview/hydro_prim_eqn.rst` — Hydrostatic Primitive Equations for the Atmosphere in Pressure Coordinates
  - L660 `doc/overview/hydrostatic.rst` — Hydrostatic, Quasi-hydrostatic, Quasi-nonhydrostatic and Non-hydrostatic forms
  - L675 `doc/overview/kinematic_bound.rst` — Kinematic Boundary conditions
  - L680 `doc/overview/ocean.rst` — Ocean
  - L683 `doc/overview/ocean_biogeo_cyc.rst` — Ocean biogeochemical cycles
  - L686 `doc/overview/ocean_gyres.rst` — Ocean gyres
  - L689 `doc/overview/overview.rst` — Overview
  - L698 `doc/overview/parm_sens.rst` — Parameter sensitivity using the adjoint of MITgcm
  - L701 `doc/overview/sim_lab_exp.rst` — Simulations of laboratory experiments
  - L704 `doc/overview/soln_strategy.rst` — Solution strategy
  - L707 `doc/overview/vector_invar.rst` — Vector invariant form
  - L710 `doc/phys_pkgs/aim.rst` — Atmospheric Intermediate Physics: AIM — `aim_v23`
  - L716 `doc/phys_pkgs/bulk_force.rst` — BULK_FORCE: Bulk Formula Package
  - L729 `doc/phys_pkgs/cal.rst` — CAL: The calendar package
  - L739 `doc/phys_pkgs/dic.rst` — DIC Package
  - L747 `doc/phys_pkgs/exch2.rst` — exch2: Extended Cubed Sphere Topology — `W2_EXCH2_TOPOLOGY.h`, `w2_e2setup.F`, `exch2_nTiles`, `W2_maxNeighbours`, `exch2_tnx`, `exch2_tny`, `exch2_tbasex`, `exch2_tbasey`, `exch2_txglobalo`, `exch2_myFace`, `exch2_nNeighbours`, `exch2_tProc`, `exch2_isWedge`, `exch2_isEedge`, `exch2_isSedge`, `exch2_isNedge`, `exch2_pi`, `exch2_pj`, `exch2_oi`, `exch2_oj`, `exch2_oi_f`, `exch2_oj_f`, `exch2_itlo_c`, `exch2_ithi_c`, `exch2_jtlo_c`, `exch2_jthi_c`
  - L760 `doc/phys_pkgs/exf.rst` — EXF: The external forcing package — `_YYYY`, `exf_inscal_`, `exf_outscal_`, `EXF_USE_INTERPOLATION`, `_lon0`, `_lon_inc`, `_lat0`, `_lat_inc`, `_nlon`, `_nlat`, `global_oce_latlon`
  - L777 `doc/phys_pkgs/fizhi.rst` — Fizhi: High-end Atmospheric Physics
  - L890 `doc/phys_pkgs/gchem.rst` — GCHEM Package
  - L899 `doc/phys_pkgs/generic_advdiff.rst` — Generic Advection/Diffusion — `gad_advection.F`, `COSINEMETH_III`, `ISOTROPIC_COS_SCALING`, `DISABLE_MULTIDIM_ADVECTION`, `GAD_MULTIDIM_COMPRESSIBLE`, `GAD_ALLOW_TS_SOM_ADV`, `GAD_SMOLARKIEWICZ_HACK`, `gad_calc_rhs.F`
  - L906 `doc/phys_pkgs/ggl90.rst` — GGL90: a TKE vertical mixing scheme
  - L911 `doc/phys_pkgs/gmredi.rst` — GMREDI: Gent-McWilliams/Redi Eddy Parameterization — `GM_AdvForm`, `GM_BOLUS_ADVEC`, `GM_EXTRA_DIAGONAL`, `uVel`, `vVel`, `wVel`, `GM_AdvSeparate`, `GM_PsiX`, `GM_PsiY`, `gmredi_calc_tensor.F`, `GEOM_alpha`, `GEOM_vert_struc`, `GEOM_lmbda`, `GEOM_diffKh_EKE`, `pickup_gmredi`, `GM_taper_scheme`, `gmredi_slope_limit.F`, `GM_NON_UNITY_DIAGONAL`, `GM_BOLUS_BVP`, `ALLOW_GM_LEITH_QG`, `GM_VISBECK_VARIABLE_K`, `GM_GEOM_VARIABLE_K`, `GMREDI_OPTIONS.h`, `gmredi_readparms.F`, `useGMREDI`, `GM_background_K`, `GM_isopycK`, `GM_maxSlope`, `GM_Kmin_horiz`, `GM_Small_Number` (+27)
  - L934 `doc/phys_pkgs/gridalt.rst` — Gridalt - Alternate Grid Package — `gridalt_initialise.F`, `make_phys_grid.F`, `gridalt_mapping.F`, `gridalt_update.F`
  - L946 `doc/phys_pkgs/kl10.rst` — KL10: Vertical Mixing Due to Breaking Internal Waves — `dumpFreq`, `taveFreq`, `viscAz`, `diffKzT`
  - L962 `doc/phys_pkgs/kpp.rst` — KPP: Nonlocal K-Profile Parameterization for Vertical Mixing — `_KPP_RL`, `FRUGAL_KPP`, `KPP_SMOOTH_SHSQ`, `KPP_SMOOTH_DVSQ`, `KPP_SMOOTH_DENS`, `KPP_SMOOTH_VISC`, `KPP_SMOOTH_DIFF`, `KPP_ESTIMATE_UREF`, `INCLUDE_DIAGNOSTICS_INTERFACE_CODE`, `KPP_GHAT`, `EXCLUDE_KPP_SHEAR_MIX`, `deltaTClock`, `dumpFreq`, `taveFreq`
  - L986 `doc/phys_pkgs/land.rst` — Land package
  - L994 `doc/phys_pkgs/mom_packages.rst` — Momentum Packages — `COSINEMETH_III`, `ISOTROPIC_COS_SCALING`, `ALLOW_SMAG_3D`, `ALLOW_3D_VISCAH`, `ALLOW_3D_VISCA4`, `ALLOW_BOTTOMDRAG_ROUGHNESS`, `MOM_BOUNDARY_CONSERVE`, `MOM_COMMON_OPTIONS.h`, `MOM_FLUXFORM_OPTIONS.h`
  - L997 `doc/phys_pkgs/obcs.rst` — OBCS: Open boundary conditions for regional modeling — `obcs_calc.F`, `obcs_fields_load.F`, `obcs_prescribe_read.F`, `OBCS_OPTIONS.h`, `packages_readparms.F`, `obcs_readparms.F`, `exf_readparms.F`, `useOBCS`, `OB_Jnorth`, `OB_Jsouth`, `OB_Ieast`, `OB_Iwest`, `useOBCSprescribe`, `useOBCSsponge`, `useOBCSbalance`, `OBCS_balanceFacN`, `OBCS_balanceFacS`, `OBCS_balanceFacE`, `OBCS_balanceFacW`, `OBCSbalanceSurf`, `useOrlanskiNorth`, `useOrlanskiSouth`, `useOrlanskiEast`, `useOrlanskiWest`, `useStevensNorth`, `useStevensSouth`, `useStevensEast`, `useStevensWest`, `cvelTimeScale`, `useFixedCEast` (+36)
  - L1015 `doc/phys_pkgs/opps.rst` — OPPS: Ocean Penetrative Plume Scheme
  - L1020 `doc/phys_pkgs/packages_overview.rst` — Using MITgcm Packages — `gmRedi`, `default_pkg_list`, `usePackageName`, `PackB`, `PackA`, `PACKAGES_CONFIG.h`, `ALLOW_GMREDI`, `useGMRedi`, `PARM_PACKAGES`
  - L1031 `doc/phys_pkgs/phys_pkgs.rst` — Packages I - Physical Parameterizations
  - L1041 `doc/phys_pkgs/ptracers.rst` — PTRACERS Package — `PTRACERS_num`, `PTRACERS_Iter0`, `PTRACERS_numInUse`, `PTRACERS_dumpFreq`, `dumpFreq`, `PTRACERS_taveFreq`, `taveFreq`, `PTRACERS_monitorFreq`, `monitorFreq`, `PTRACERS_timeave_mnc`, `useMNC`, `timeave_mnc`, `PTRACERS_snapshot_mnc`, `snapshot_mnc`, `PTRACERS_monitor_mnc`, `monitor_mnc`, `PTRACERS_pickup_write_mnc`, `pickup_write_mnc`, `PTRACERS_pickup_read_mnc`, `pickup_read_mnc`, `PTRACERS_useRecords`, `PTRACERS_advScheme`, `saltAdvScheme`, `PTRACERS_ImplVertAdv`, `PTRACERS_diffKh`, `diffKhS`, `PTRACERS_diffK4`, `diffK4S`, `PTRACERS_diffKr`, `PTRACERS_diffKrNr` (+10)
  - L1050 `doc/phys_pkgs/rbcs.rst` — RBCS Package — `maskLEN`, `rbcsForcingPeriod`, `rbcsForcingCycle`, `rbcsForcingOffset`, `rbcsSingleTimeFiles`, `deltaTrbcs`, `deltaTclock`, `rbcsVanishingTime`, `myTime`, `rbcsIter0`, `useRBCtemp`, `useRBCsalt`, `useRBCuVel`, `useRBCvVel`, `tauRelaxT`, `tauRelaxS`, `tauRelaxU`, `tauRelaxV`, `relaxMaskFile`, `relaxMaskUFile`, `relaxMaskVFile`, `relaxTFile`, `relaxSFile`, `relaxUFile`, `relaxVFile`, `useRBCpTrNum`, `tauRelaxPTR`, `relaxPtracerFile`, `useRBCxxx`, `RBCS_SIZE.h`
  - L1065 `doc/phys_pkgs/remesh.rst` — SHELFICE Remeshing — `NONLIN_FRSURF`, `ALLOW_SHELFICE_REMESHING`, `SHI_ALLOW_GAMMAFRICT`, `SHI_withBL_uStarTopDz`, `CPP_OPTIONS.h`, `SHELFICE_OPTIONS.h`, `nonlinFreeSurf`, `select_rstar`, `SHI_withBL_realFWflux`, `SHELFICEboundaryLayer`, `SHELFICEmassFile`, `SHELFICEMassStepping`, `SHELFICEMassDynTendFile`, `SHELFICEDynMassOnly`, `shelficeMass`, `SHELFICERemeshFrequency`, `SHELFICESplitThreshold`, `hFacC`, `SHELFICEMergeThreshold`, `implicitFreeSurface`, `useRealFreshWaterFlux`, `SHELFICEuseGammaFrict`, `conserve_ssh`, `OBCS_BALANCE_FLOW`, `useOBCSbalance`
  - L1077 `doc/phys_pkgs/seaice.rst` — SEAICE Package — `SEAICE_BGRID_DYNAMICS`, `SEAICE_DEBUG`, `SEAICE_CGRID`, `SEAICE_ALLOW_EVP`, `SEAICE_ALLOW_JFNK`, `SEAICE_ALLOW_KRYLOV`, `SEAICE_ALLOW_TEM`, `SEAICE_ALLOW_MCS`, `SEAICE_ALLOW_MCE`, `SEAICE_ALLOW_TD`, `SEAICE_LSR_ZEBRA`, `SEAICE_ALLOW_FREEDRIFT`, `SEAICE_EXTERNAL_FLUXES`, `SEAICE_ZETA_SMOOTHREG`, `SEAICE_DELTA_SMOOTHREG`, `SEAICE_ALLOW_BOTTOMDRAG`, `SEAICE_ALLOW_SIDEDRAG`, `SEAICE_BICE_STRESS`, `EXPLICIT_SSH_SLOPE`, `SEAICE_LSRBNEW`, `SEAICE_ITD`, `SEAICE_VARIABLE_SALINITY`, `SEAICE_CAP_ICELOAD`, `siceLoad`, `ALLOW_SITRACER`, `SEAICE_USE_GROWTH_ADX`, `SEAICE_OPTIONS.h`, `seaice_growth_adx.F`, `seaice_growth.F`, `seaice_readparms.F` (+100)
  - L1113 `doc/phys_pkgs/shap_filt.rst` — Shapiro Filter
  - L1118 `doc/phys_pkgs/shelfice.rst` — SHELFICE Package — `ALLOW_SHELFICE_DEBUG`, `ALLOW_ISOMIP_TD`, `SHI_ALLOW_GAMMAFRICT`, `SHI_SALTBAL_FWFLX`, `SHELFICE_OPTIONS.h`, `useSHELFICE`, `SHELFICEtopoFile`, `SHELFICEloadAnomalyFile`, `rhoConst`, `useISOMIPTD`, `SHELFICEboundaryLayer`, `SHI_withBL_realFWflux`, `SHI_withBL_uStarTopDz`, `SHELFICEmassFile`, `SHELFICEMassDynTendFile`, `SHELFICETransCoeffTFile`, `SHELFICElatentHeat`, `SHELFICEHeatCapacity_Cp`, `rhoShelfIce`, `SHELFICEheatTransCoeff`, `SHELFICEsaltTransCoeff`, `SHELFICEsaltToHeatRatio`, `SHELFICEthetaSurface`, `no_slip_shelfice`, `no_slip_bottom`, `SHELFICEDragLinear`, `bottomDragLinear`, `SHELFICEDragQuadratic`, `bottomDragQuadratic`, `SHELFICEselectDragQuadr` (+21)
  - L1133 `doc/phys_pkgs/streamice.rst` — STREAMICE Package — `STREAMICE_CONSTRUCT_MATRIX`, `STREAMICE_HYBRID_STRESS`, `USE_ALT_RLOW`, `STREAMICE_GEOM_FILE_SETUP`, `STREAMICE_PARM03`, `ALLOW_PETSC`, `STREAMICE_COULOMB_SLIDING`, `STREAMICE_SMOOTH_FLOATATION`, `STREAMICE_OPTIONS.h`, `useSTREAMICE`, `STREAMICE_PARM01`, `streamice_density`, `streamice_density_ocean_avg`, `n_glen`, `eps_glen_min`, `eps_u_min`, `n_basal_friction`, `streamice_cg_tol`, `streamice_lower_cg_tol`, `streamice_max_cg_iter`, `streamice_maxcgiter_cpl`, `streamice_nonlin_tol`, `streamice_max_nl_iter`, `streamice_maxnliter_cpl`, `streamice_nonlin_tol_fp`, `streamice_err_norm`, `streamice_chkfixedptconvergence`, `streamice_chkresidconvergence`, `streamicethickInit`, `streamicethickFile` (+83)
  - L1156 `doc/phys_pkgs/thsice.rst` — THSICE: The Thermodynamic Sea Ice Package
  - L1173 `doc/phys_pkgs/zonal_filt.rst` — FFT Filtering Code
  - L1178 `doc/references.rst` —
  - L1180 `doc/related_projects/related_projects.rst` — Related Projects and Highlighted Papers
  - L1193 `doc/software_arch/software_arch.rst` — Software Architecture — `myThid`, `cg2d_r`, `SIZE.h`, `dynamics.F`, `main.F`, `the_model_main.F`, `myProcId`, `MPI_COMM_MODEL`, `myXGlobalLo`, `myYGlobalLo`, `pidW`, `pidE`, `pidS`, `pidN`, `eeboot_minimal.F`, `ini_procs.F`, `EESUPPORT.h`, `tileNo`, `tileNoN`, `tileNoS`, `tileNoE`, `tileNoW`, `tileCommModeN`, `tileCommModeS`, `tileCommModeE`, `tileCommModeW`, `nThreads`, `exchNeedsMemSync`, `cacheLineSize`, `lShare1` (+30)
  - L1230 `doc/utilities/utilities.rst` — Utilities

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/packages.md`
- L1 MITgcm packages — `aim_v23`, `atm_common`, `atm_compon_interf`, `atm_ocn_coupler`, `atm_phys`, `bulk_force`, `cd_code`, `compon_communic`, `down_slope`, `embed_files`, `generic_advdiff`, `mom_common`, `mom_fluxform`, `mom_vecinv`, `ocn_compon_interf`, `salt_plume`, `shap_filt`, `steep_icecavity`, `zonal_filt`
  - L114 Core (model/) changes in the user's branches

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/tools.md`
- L1 tools: genmake2 and testreport
  - L3 genmake2 -help
  - L181 testreport -help
  - L265 Optfiles shipped (tools/build_options) — `cygwin_ia32_g77`, `darwin_absoft_f77`, `darwin_amd64_gfortran`, `darwin_arm64_gfortran`, `darwin_ia32_g95`, `darwin_ia32_gfortran`, `darwin_ia32_ifort`, `darwin_ia32_pgf95_trane`, `darwin_ppc_f95`, `darwin_ppc_g77`, `darwin_ppc_xlf`, `darwin_ppc_xlf_panther`, `darwin_ppc_xlf_panther_baylor`, `darwin_ppc_xlf_tiger_baylor`, `linux_alpha_g77`, `linux_amd64_absoft`, `linux_amd64_g77`, `linux_amd64_g95`, `linux_amd64_gfortran`, `linux_amd64_gfortran_albedo`, `linux_amd64_gfortran_avx2`, `linux_amd64_gfortran_greenplanet`, `linux_amd64_ifort`, `linux_amd64_ifort11`, `linux_amd64_ifort_albedo`, `linux_amd64_ifort_beagle`, `linux_amd64_ifort_discover`, `linux_amd64_ifort_fimm_default`, `linux_amd64_ifort_fimm_emic`, `linux_amd64_ifort_uv100` (+37)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification.md`
- L1 Verification experiments — `advect_cs`, `advect_xy`, `advect_xz`, `aim.5l_Equatorial_Channel`, `aim.5l_LatLon`, `aim.5l_cs`, `atm_gray`, `bottom_ctrl_5x5`, `cfc_example`, `cheapAML_box`, `deep_anelastic`, `front_relax`, `global_oce_biogeo_bling`, `global_oce_latlon`, `global_ocean.90x40x15`, `global_ocean.cs32x15`, `halfpipe_streamice`, `ideal_2D_oce`, `internal_wave`, `inverted_barometer`, `lab_sea`, `matrix_example`, `obcs_ctrl`, `offline_exf_seaice`, `seaice_itd`, `seaice_obcs`, `shelfice_2d_remesh`, `short_surf_wave`, `so_box_biogeo`, `tutorial_advection_in_gyre` (+14)
  - L70 Experiments in the user's branches (new or modified vs upstream)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/admtlm.md`
- L1 pkg/admtlm — `useADMTLM`
  - L7 Headers — `ADMTLM_OPTIONS.h`, `arpack_debug.h`
  - L11 Routines (8) — `admtlm_bypassad.F`, `admtlm_driver.F`, `admtlm_dsvd.F`, `admtlm_dsvd2model.F`, `admtlm_init_fixed.F`, `admtlm_map.F`, `admtlm_metric.F`, `admtlm_model2dsvd.F`
  - L14 Called from outside the package — `ADMTLM_DSVD`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/aim_v23.md`
- L1 pkg/aim_v23 — `useAIM_V23`
  - L10 Namelist parameters
    - L11 AIM_PARAMS — `aim_useFMsurfBC`, `aim_useMMsurfFc`, `aim_surfForc_TimePeriod`, `aim_surfForc_NppCycle`, `aim_surfForc_TransRatio`, `aim_surfPotTemp`, `aim_energPrecip`, `aim_splitSIOsFx`, `aim_MMsufx`, `aim_MMsufxLength`, `aim_LandFile`, `aim_albFile`, `aim_vegFile`, `aim_sstFile`, `aim_lstFile`, `aim_oiceFile`, `aim_snowFile`, `aim_swcFile`, `aim_qfxFile`, `aim_dragStrato`, `aim_selectOceAlbedo`, `aim_select_pCO2`, `aim_abs_pCO2`, `aim_fixed_pCO2`, `atmpCO2init`, `aim_clrSkyDiag`, `aim_taveFreq`, `aim_diagFreq`, `aim_tendFreq`, `aim_timeave_mnc` (+3)
    - L45 AIM_PAR_FOR — `hSnowWetness`
    - L55 AIM_PAR_SFL — `dTstab`
    - L70 AIM_PAR_CNV
    - L78 AIM_PAR_LSC
    - L83 AIM_PAR_RAD
    - L104 AIM_PAR_VDI
  - L111 CPP options (defaults as shipped) — `ALLOW_DEW_ON_LAND`, `ALLOW_INSOLATION`, `ALLOW_CLOUD_3D`, `ALLOW_AIM_CO2`, `ALLOW_CLR_SKY_DIAG`, `_KD2KA`
  - L119 Headers — `AIM2DYN.h`, `AIM_CO2.h`, `AIM_FFIELDS.h`, `AIM_GRID.h`, `AIM_OPTIONS.h`, `AIM_PARAMS.h`, `AIM_SIZE.h`, `com_cnvcon.h`, `com_forcing.h`, `com_forcon.h`, `com_lsccon.h`, `com_physcon.h`, `com_physvar.h`, `com_radcon.h`, `com_radvar.h`, `com_sflcon.h`, `com_vdicon.h`, `phy_const.h`
  - L139 Routines (42) — `aim_aim2dyn.F`, `aim_aim2land.F`, `aim_aim2sioce.F`, `aim_diagnostics.F`, `aim_diagnostics_init.F`, `aim_do_co2.F`, `aim_do_physics.F`, `aim_dyn2aim.F`, `aim_fields_load.F`, `aim_initialise.F`, `aim_land2aim.F`, `aim_land_impl.F`, `aim_mnc_init.F`, `aim_readparms.F`, `aim_sice2aim.F`, `aim_sice_impl.F`, `aim_surf_bc.F`, `aim_tendency_apply.F`, `aim_write_local.F`, `aim_write_phys.F`, `aim_write_tave.F`, `phy_convmf.F`, `phy_driver.F`, `phy_inphys.F`, `phy_lscond.F`, `phy_radiat.F`, `phy_shtorh.F`, `phy_snow_precip.F`, `phy_suflux_land.F`, `phy_suflux_ocean.F` (+4)
  - L142 Called from outside the package — `AIM_TENDENCY_APPLY_S`, `AIM_TENDENCY_APPLY_T`, `AIM_TENDENCY_APPLY_U`, `AIM_TENDENCY_APPLY_V`, `AIM_DO_PHYSICS`, `AIM_FIELDS_LOAD`, `AIM_INITIALISE`
  - L155 Verification experiments compiling it (4) — `aim.5l_Equatorial_Channel`, `aim.5l_LatLon`, `aim.5l_cs`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/atm2d.md`
- L1 pkg/atm2d — `useATM2D`
  - L8 Namelist parameters
    - L9 COUPLE_PARM — `startYear`, `endYear`, `taveDump`
    - L16 PARM01_ATM2D — `atmosTauuFile`, `atmosTauvFile`, `atmosWindFile`
    - L20 PARM02_ATM2D — `tauuFile`, `tauvFile`, `windFile`, `qnetFile`, `evapFile`, `precipFile`
    - L27 PARM03_ATM2D — `thetaRelaxFile`, `saltRelaxFile`, `tauThetaRelax`, `tauSaltRelax`
    - L34 PARM04_ATM2D — `runoffFile`, `runoffMapFile`
    - L39 PARM05_ATM2D — `useObsEmP`, `useObsRunoff`, `useAltDeriv`
  - L44 CPP options (defaults as shipped) — `ATM2D_MPI_ON`, `ALLOW_DBUG_ATM2D`, `CPL_CHEM`, `CPL_NEM`, `CPL_OCEANCO2`, `CPL_TEM`, `IPCC_EMI`, `ML_2D`, `OCEAN_3D`, `ORBITAL_FOR`, `PREDICTED_AEROSOL`
  - L63 Headers — `AGRID.h`, `ATM2D_OPTIONS.h`, `ATM2D_VARS.h`, `ATMSIZE.h`, `CPLIDS.h`, `DRIVER.h`, `OCNIDS.h`, `OCNSIZE.h`, `OCNVARS.h`, `ctrparam.h`
  - L75 Routines (31) — `atm2d_finish.F`, `atm2d_init_fixed.F`, `atm2d_init_vars.F`, `atm2d_read_pickup.F`, `atm2d_readparms.F`, `atm2d_write_pickup.F`, `atm2ocn_main.F`, `calc_1dto2d.F`, `calc_fileload.F`, `calc_zonal_means.F`, `fixed_flux_add.F`, `forward_step_atm2d.F`, `get_ocnvars.F`, `init_atm2d.F`, `init_sumvars.F`, `month_end_diags.F`, `norm_ocn_fluxes.F`, `pass_thsice_fluxes.F`, `put_ocnvars.F`, `read_atmos.F`, `relax_add.F`, `subtract_means.F`, `sum_ocn_fluxes.F`, `sum_thsice_out.F`, `sum_yr_end_diags.F`, `tave_end_diags.F`, `yr_end_diags.F`
  - L78 Called from outside the package — `FORWARD_STEP_ATM2D`, `ATM2D_INIT_FIXED`, `ATM2D_INIT_VARS`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/atm_common.md`
- L1 pkg/atm_common — `useATM_COMMON`
  - L7 Headers — `ATM_COMMON_OPTIONS.h`
  - L10 Routines (1) — `atm_diagnostics_init.F`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/atm_compon_interf.md`
- L1 pkg/atm_compon_interf — `useATM_COMPON_INTERF`
  - L8 Namelist parameters
    - L9 CPL_ATM_PARAM — `cpl_earlyExpImpCall`, `cpl_oldPickup`, `useImportMxlD`, `useImportSST`, `useImportSSS`, `useImportVsq`, `useImportThSIce`, `useImportFlxCO2`, `cpl_atmSendFrq`, `maxNumberPrint`
  - L21 Headers — `ATMCPL.h`, `ATM_CPL_OPTIONS.h`, `CPL_PARAMS.h`
  - L26 Routines (17) — `atm_apply_import.F`, `atm_check_cplconfig.F`, `atm_cpl_read_pickup.F`, `atm_export_atmconfig.F`, `atm_export_fields.F`, `atm_export_fld.F`, `atm_get_atmconfig.F`, `atm_import_fields.F`, `atm_import_ocnconfig.F`, `atm_store_aim_fields.F`, `atm_store_aim_wndstr.F`, `atm_store_dynvars.F`, `atm_store_land.F`, `atm_store_my_data.F`, `atm_store_surfflux.F`, `atm_store_thsice.F`, `cpl_output.F`
  - L29 Called from outside the package — `ATM_STORE_MY_DATA`, `ATM_APPLY_IMPORT`
  - L35 Verification experiments compiling it (1)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/atm_ocn_coupler.md`
- L1 pkg/atm_ocn_coupler — `useATM_OCN_COUPLER`
  - L7 README
  - L41 Namelist parameters
    - L42 COUPLER_PARAMS — `cpl_sequential`, `cpl_exchange_RunOff`, `cpl_exchange1W_sIce`, `cpl_exchange2W_sIce`, `cpl_exchange_SaltPl`, `cpl_exchange_DIC`, `runOffMapSize`, `runOffMapFile`
  - L52 Headers — `ATMIDS.h`, `ATMVARS.h`, `CPLIDS.h`, `CPL_MAP2GRIDS.h`, `CPL_PARAMS.h`, `CPP_EEOPTIONS.h`, `CPP_OPTIONS.h`, `OCNIDS.h`, `OCNVARS.h`
  - L63 Routines (25) — `accept_component_registrations.F`, `atm_to_ocn_maprunoff.F`, `atm_to_ocn_mapxyr8.F`, `cpl_check_cplconfig.F`, `cpl_init_atm_vars.F`, `cpl_init_ocn_vars.F`, `cpl_read_params.F`, `cpl_recv_atm_atmconfig.F`, `cpl_recv_atm_fields.F`, `cpl_recv_ocn_fields.F`, `cpl_recv_ocn_ocnconfig.F`, `cpl_register_atm.F`, `cpl_register_ocn.F`, `cpl_send_atm_cplparms.F`, `cpl_send_atm_fields.F`, `cpl_send_atm_ocnconfig.F`, `cpl_send_ocn_atmconfig.F`, `cpl_send_ocn_cplparms.F`, `cpl_send_ocn_fields.F`, `exch_component_configs.F`, `initialise.F`, `mds_byteswap.F`, `ocn_to_atm_mapxyr8.F`, `set_runoffmap.F`
  - L66 Called from outside the package — `CPL_REGISTER_OCN`, `ACCEPT_COMPONENT_REGISTRATIONS`, `EXCH_COMPONENT_CONFIGS`, `CPL_RECV_OCN_OCNCONFIG`, `CPL_SEND_OCN_ATMCONFIG`, `CPL_RECV_OCN_FIELDS`, `CPL_SEND_OCN_FIELDS`, `MDS_BYTESWAPR8`, `MDS_BYTESWAPR4`
  - L116 Verification experiments compiling it (1)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/atm_phys.md`
- L1 pkg/atm_phys — `useATM_PHYS`, `data.atm_gray`, `data.atm_phys`
  - L9 Namelist parameters
    - L10 ATMOSPHERE_NML — `ldry_convection`, `lwet_convection`, `do_virtual`, `two_stream`, `mixed_layer_bc`, `roughness_heat`, `roughness_moist`, `roughness_mom`
    - L20 ATM_PHYS_PARM01 — `atmPhys_addTendT`, `atmPhys_addTendS`, `atmPhys_addTendU`, `atmPhys_addTendV`, `atmPhys_tauDampUV`, `atmPhys_dampUVfac`, `atmPhys_stepSST`, `atmPhys_sstFile`, `atmPhys_qFlxFile`, `atmPhys_mxldFile`, `atmPhys_albedoFile`, `atmPhys_ozoneFile`
    - L33 MY25_TURB_NML
    - L34 SHALLOW_CONV_NML
    - L35 VERT_TURB_DRIVER_NML — `do_shallow_conv`, `do_mellor_yamada`
  - L39 Headers — `ATM_PHYS_OPTIONS.h`, `ATM_PHYS_PARAMS.h`, `ATM_PHYS_VARS.h`
  - L44 Routines (36) — `atm_phys_check.F`, `atm_phys_diagnostics_init.F`, `atm_phys_driver.F`, `atm_phys_dyn2phys.F`, `atm_phys_init_fixed.F`, `atm_phys_init_varia.F`, `atm_phys_read_pickup.F`, `atm_phys_readparms.F`, `atm_phys_tendency_apply.F`, `atm_phys_write_pickup.F`, `dargan_bettsmiller_mod.F90`, `lscale_cond_mod.F90`, `my25_turb_mod.F90`, `shallow_conv_mod.F90`, `simple_sat_vapor_pres_mod.F90`
  - L47 Called from outside the package — `ATM_PHYS_TENDENCY_APPLY_S`, `ATM_PHYS_TENDENCY_APPLY_T`, `ATM_PHYS_TENDENCY_APPLY_U`, `ATM_PHYS_TENDENCY_APPLY_V`, `ATM_PHYS_DRIVER`, `ATM_PHYS_CHECK`, `ATM_PHYS_INIT_FIXED`, `ATM_PHYS_INIT_VARIA`, `ATM_PHYS_READPARMS`, `ATM_PHYS_WRITE_PICKUP`
  - L63 Verification experiments compiling it (1) — `atm_gray`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/autodiff.md`
- L1 pkg/autodiff — `useAUTODIFF`
  - L12 README
  - L46 Namelist parameters
    - L47 AUTODIFF_PARM01 — `dumpAdVarExch`, `mon_AdVarExch`, `dumpAdByRec`, `useApproxAdvectionInAdMode`, `useKPPinAdMode`, `useGGL90inAdMode`, `useGMRediInAdMode`, `useSALT_PLUMEinAdMode`, `useSEAICEinAdMode`, `useSmoothCorrel2DinAdMode`, `inAdExact`, `SEAICEapproxLevInAd`, `viscFacInAd`, `viscFacInFw`, `SIregFacInAd`, `SIregFacInFw`, `SEAICEuseFREEDRIFTswitchInAd`, `SEAICEuseDYNAMICSswitchInAd`, `cg2dFullAdjoint`
  - L68 CPP options (defaults as shipped) — `ALLOW_AUTODIFF_TAMC`, `AUTODIFF_TAMC_COMPATIBILITY`, `ALLOW_TAMC_CHECKPOINTING`, `ALLOW_AUTODIFF_MONITOR`, `ALLOW_AUTODIFF_MONITOR_DIAG`, `AUTODIFF_2_LEVEL_CHECKPOINT`, `AUTODIFF_4_LEVEL_CHECKPOINT`, `ALLOW_DIVIDED_ADJOINT`, `ALLOW_AUTODIFF_WHTAPEIO`, `ALLOW_INIT_WHTAPEIO`, `AUTODIFF_USE_STORE_RESTORE`, `AUTODIFF_USE_STORE_RESTORE_OBCS`, `AUTODIFF_ALLOW_VISCFACADJ`, `AUTODIFF_DISABLE_LEITH`, `AUTODIFF_DISABLE_REYNOLDS_SCALE`, `AUTODIFF_EXCLUDE_ADEXCH_RS`
  - L86 Headers — `AUTODIFF.h`, `AUTODIFF_OPTIONS.h`, `AUTODIFF_PARAMS.h`, `AUTODIFF_STORE.h`, `ECCO_CPPOPTIONS.h`, `adcommon.h`, `check_lev1_dir_forcing.h`, `check_lev1_dir_ocephys.h`, `checkpoint_lev1_directives.h`, `checkpoint_lev1_template.h`, `checkpoint_lev2_directives.h`, `checkpoint_lev3_directives.h`, `checkpoint_lev4_directives.h`, `g_common.h`, `tamc.h`
  - L103 Routines (139) — `active_file.F`, `active_file_ad.F`, `active_file_control.F`, `active_file_control_slice.F`, `active_file_g.F`, `active_file_gen.F`, `active_file_gen_ad.F`, `active_file_gen_g.F`, `active_file_loc.F`, `active_file_loc_ad.F`, `active_file_loc_g.F`, `adautodiff_whtapeio_sync.F`, `add_prefix.F`, `addamp_adj.F`, `addummy_for_etan.F`, `addummy_in_dynamics.F`, `addummy_in_stepping.F`, `adopen_adclose.F`, `adread_adwrite.F`, `adread_adwrite_i.F`, `adzero_adj.F`, `autodiff_check.F`, `autodiff_findunit.F`, `autodiff_inadmode_set.F`, `autodiff_inadmode_set_ad.F`, `autodiff_inadmode_set_g.F`, `autodiff_inadmode_unset.F`, `autodiff_inadmode_unset_ad.F`, `autodiff_inadmode_unset_g.F`, `autodiff_ini_model_io.F` (+26)
  - L106 Called from outside the package — `ZERO_ADJ_LOC`, `DUMMY_IN_DYNAMICS`, `AUTODIFF_INADMODE_SET`, `AUTODIFF_INADMODE_UNSET`, `DUMMY_IN_STEPPING`, `AUTODIFF_INI_MODEL_IO`, `AUTODIFF_INIT_VARIA`, `DUMMY_FOR_ETAN`, `AUTODIFF_CHECK`, `AUTODIFF_READPARMS`, `CG2D_STORE`, `AUTODIFF_RESTORE`, `AUTODIFF_STORE`, `AUTODIFF_WHTAPEIO_SYNC`, `ADACTIVE_READ_XYZ`, `ADACTIVE_READ_XYZ_LOC`, `ADACTIVE_READ_XY_LOC`, `ACTIVE_WRITE_XY`, `ACTIVE_READ_XY`, `ACTIVE_READ_XYZ`, `ACTIVE_READ_YZ`, `ACTIVE_READ_XZ`, `ACTIVE_WRITE_XYZ`, `ACTIVE_WRITE_GEN_RS`, `ACTIVE_WRITE_XZ`, `ACTIVE_WRITE_YZ`, `DUMP_ADJ_XY`, `DUMP_ADJ_XY_UV`, `COPY_ADVAR_OUTP`, `COPY_AD_UV_OUTP` (+9)
  - L211 Verification experiments compiling it (16) — `bottom_ctrl_5x5`, `global_oce_biogeo_bling`, `global_oce_latlon`, `global_ocean.90x40x15`, `global_ocean.cs32x15`, `halfpipe_streamice`, `lab_sea`, `obcs_ctrl`, `offline_exf_seaice`, `tutorial_dic_adjoffline`, `tutorial_global_oce_biogeo`, `tutorial_global_oce_optim`, `tutorial_tracer_adjsens`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/autodiff@bbl.md`
- L1 pkg/autodiff  (bbl: ~/Documents/research/ECCO/BBL/MITgcm) — `useAUTODIFF`
  - L14 README
  - L48 Namelist parameters
    - L49 AUTODIFF_PARM01 — `dumpAdVarExch`, `mon_AdVarExch`, `dumpAdByRec`, `useApproxAdvectionInAdMode`, `useKPPinAdMode`, `useGGL90inAdMode`, `useGMRediInAdMode`, `useSALT_PLUMEinAdMode`, `useSEAICEinAdMode`, `useSmoothCorrel2DinAdMode`, `inAdExact`, `SEAICEapproxLevInAd`, `viscFacInAd`, `viscFacInFw`, `SIregFacInAd`, `SIregFacInFw`, `SEAICEuseFREEDRIFTswitchInAd`, `SEAICEuseDYNAMICSswitchInAd`, `cg2dFullAdjoint`
  - L70 CPP options (defaults as shipped) — `ALLOW_AUTODIFF_TAMC`, `AUTODIFF_TAMC_COMPATIBILITY`, `ALLOW_TAMC_CHECKPOINTING`, `ALLOW_AUTODIFF_MONITOR`, `ALLOW_AUTODIFF_MONITOR_DIAG`, `AUTODIFF_2_LEVEL_CHECKPOINT`, `AUTODIFF_4_LEVEL_CHECKPOINT`, `ALLOW_DIVIDED_ADJOINT`, `ALLOW_AUTODIFF_WHTAPEIO`, `ALLOW_INIT_WHTAPEIO`, `AUTODIFF_USE_STORE_RESTORE`, `AUTODIFF_USE_STORE_RESTORE_OBCS`, `AUTODIFF_ALLOW_VISCFACADJ`, `AUTODIFF_DISABLE_LEITH`, `AUTODIFF_DISABLE_REYNOLDS_SCALE`, `AUTODIFF_EXCLUDE_ADEXCH_RS`
  - L88 Headers — `AUTODIFF.h`, `AUTODIFF_OPTIONS.h`, `AUTODIFF_PARAMS.h`, `AUTODIFF_STORE.h`, `ECCO_CPPOPTIONS.h`, `adcommon.h`, `check_lev1_dir_forcing.h`, `check_lev1_dir_ocephys.h`, `checkpoint_lev1_directives.h`, `checkpoint_lev1_template.h`, `checkpoint_lev2_directives.h`, `checkpoint_lev3_directives.h`, `checkpoint_lev4_directives.h`, `g_common.h`, `tamc.h`
  - L105 Routines (139) — `active_file.F`, `active_file_ad.F`, `active_file_control.F`, `active_file_control_slice.F`, `active_file_g.F`, `active_file_gen.F`, `active_file_gen_ad.F`, `active_file_gen_g.F`, `active_file_loc.F`, `active_file_loc_ad.F`, `active_file_loc_g.F`, `adautodiff_whtapeio_sync.F`, `add_prefix.F`, `addamp_adj.F`, `addummy_for_etan.F`, `addummy_in_dynamics.F`, `addummy_in_stepping.F`, `adopen_adclose.F`, `adread_adwrite.F`, `adread_adwrite_i.F`, `adzero_adj.F`, `autodiff_check.F`, `autodiff_findunit.F`, `autodiff_inadmode_set.F`, `autodiff_inadmode_set_ad.F`, `autodiff_inadmode_set_g.F`, `autodiff_inadmode_unset.F`, `autodiff_inadmode_unset_ad.F`, `autodiff_inadmode_unset_g.F`, `autodiff_ini_model_io.F` (+26)
  - L108 Called from outside the package — `ZERO_ADJ_LOC`, `DUMMY_IN_DYNAMICS`, `AUTODIFF_INADMODE_SET`, `AUTODIFF_INADMODE_UNSET`, `DUMMY_IN_STEPPING`, `AUTODIFF_INI_MODEL_IO`, `AUTODIFF_INIT_VARIA`, `DUMMY_FOR_ETAN`, `AUTODIFF_CHECK`, `AUTODIFF_READPARMS`, `CG2D_STORE`, `AUTODIFF_RESTORE`, `AUTODIFF_STORE`, `AUTODIFF_WHTAPEIO_SYNC`, `ADACTIVE_READ_XYZ`, `ADACTIVE_READ_XYZ_LOC`, `ADACTIVE_READ_XY_LOC`, `ACTIVE_WRITE_XY`, `ACTIVE_READ_XY`, `ACTIVE_READ_XYZ`, `ACTIVE_READ_YZ`, `ACTIVE_READ_XZ`, `ACTIVE_WRITE_XYZ`, `ACTIVE_WRITE_GEN_RS`, `ACTIVE_WRITE_XZ`, `ACTIVE_WRITE_YZ`, `DUMP_ADJ_XY`, `DUMP_ADJ_XY_UV`, `COPY_ADVAR_OUTP`, `COPY_AD_UV_OUTP` (+9)
  - L213 Verification experiments compiling it (1)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/autodiff@bbl_c68g.md`
- L1 pkg/autodiff  (bbl_c68g: ~/Documents/research/ECCO/BBL/MITgcm_c68g) — `useAUTODIFF`
  - L14 README
  - L48 Namelist parameters
    - L49 AUTODIFF_PARM01 — `dumpAdVarExch`, `mon_AdVarExch`, `dumpAdByRec`, `useApproxAdvectionInAdMode`, `useKPPinAdMode`, `useGGL90inAdMode`, `useGMRediInAdMode`, `useSALT_PLUMEinAdMode`, `useSEAICEinAdMode`, `useSmoothCorrel2DinAdMode`, `inAdExact`, `SEAICEapproxLevInAd`, `viscFacInAd`, `viscFacInFw`, `SIregFacInAd`, `SIregFacInFw`, `SEAICEuseFREEDRIFTswitchInAd`, `SEAICEuseDYNAMICSswitchInAd`, `cg2dFullAdjoint`
  - L70 CPP options (defaults as shipped) — `ALLOW_AUTODIFF_TAMC`, `ALLOW_TAMC_CHECKPOINTING`, `ALLOW_AUTODIFF_MONITOR`, `ALLOW_AUTODIFF_MONITOR_DIAG`, `ALLOW_DIVIDED_ADJOINT`, `ALLOW_PACKUNPACK_METHOD2`, `ALLOW_AUTODIFF_WHTAPEIO`, `AUTODIFF_USE_STORE_RESTORE`, `AUTODIFF_USE_STORE_RESTORE_OBCS`, `AUTODIFF_PTRACERS_SPLIT_FILES`, `AUTODIFF_ALLOW_VISCFACADJ`, `AUTODIFF_EXCLUDE_ADEXCH_RS`
  - L84 Headers — `AUTODIFF.h`, `AUTODIFF_MYFIELDS.h`, `AUTODIFF_OPTIONS.h`, `AUTODIFF_PARAMS.h`, `AUTODIFF_STORE.h`, `ECCO_CPPOPTIONS.h`, `adcommon.h`, `check_lev1_dir_forcing.h`, `check_lev1_dir_ocephys.h`, `checkpoint_lev1_directives.h`, `checkpoint_lev1_template.h`, `checkpoint_lev2_directives.h`, `checkpoint_lev3_directives.h`, `checkpoint_lev4_directives.h`, `g_common.h`, `tamc.h`
  - L102 Routines (135) — `active_file.F`, `active_file_ad.F`, `active_file_control.F`, `active_file_control_slice.F`, `active_file_g.F`, `active_file_gen.F`, `active_file_gen_ad.F`, `active_file_gen_g.F`, `active_file_loc.F`, `active_file_loc_ad.F`, `active_file_loc_g.F`, `adautodiff_whtapeio_sync.F`, `addamp_adj.F`, `addummy_for_etan.F`, `addummy_in_dynamics.F`, `addummy_in_stepping.F`, `adopen_adclose.F`, `adread_adwrite.F`, `adread_adwrite_i.F`, `adzero_adj.F`, `autodiff_check.F`, `autodiff_findunit.F`, `autodiff_inadmode_set.F`, `autodiff_inadmode_set_ad.F`, `autodiff_inadmode_set_g.F`, `autodiff_inadmode_unset.F`, `autodiff_inadmode_unset_ad.F`, `autodiff_inadmode_unset_g.F`, `autodiff_ini_model_io.F`, `autodiff_init_varia.F` (+24)
  - L105 Called from outside the package — `ZERO_ADJ_LOC`, `DUMMY_IN_DYNAMICS`, `AUTODIFF_INADMODE_SET`, `AUTODIFF_INADMODE_UNSET`, `DUMMY_IN_STEPPING`, `AUTODIFF_INI_MODEL_IO`, `AUTODIFF_INIT_VARIA`, `DUMMY_FOR_ETAN`, `AUTODIFF_CHECK`, `AUTODIFF_READPARMS`, `CG2D_STORE`, `AUTODIFF_RESTORE`, `AUTODIFF_STORE`, `AUTODIFF_WHTAPEIO_SYNC`, `ADACTIVE_READ_XYZ`, `ADACTIVE_READ_XYZ_LOC`, `ADACTIVE_READ_XY_LOC`, `ACTIVE_WRITE_XY`, `ACTIVE_READ_XY`, `ACTIVE_READ_XYZ`, `ACTIVE_READ_YZ`, `ACTIVE_READ_XZ`, `ACTIVE_WRITE_XYZ`, `ACTIVE_WRITE_GEN_RS`, `ACTIVE_WRITE_XZ`, `ACTIVE_WRITE_YZ`, `ACTIVE_WRITE_XY_LOC`, `DUMP_ADJ_XY`, `DUMP_ADJ_XY_UV`, `COPY_ADVAR_OUTP` (+4)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/bbl.md`
- L1 pkg/bbl — `useBBL`
  - L9 README — `bbl_tracer`, `bbl_rho`, `bbl_eta`, `resTracer`
  - L43 Namelist parameters
    - L44 BBL_PARM01 — `bbl_wvel`, `bbl_hvel`, `bbl_initEta`, `bbl_thetaFile`, `bbl_saltFile`, `bbl_etaFile`
  - L52 Headers — `BBL.h`, `BBL_OPTIONS.h`
  - L56 Routines (17) — `bbl_calc_rho.F`, `bbl_calc_rhs.F`, `bbl_check.F`, `bbl_diagnostics_init.F`, `bbl_diagnostics_state.F`, `bbl_init_fixed.F`, `bbl_init_varia.F`, `bbl_read_pickup.F`, `bbl_readparms.F`, `bbl_tendency_apply.F`, `bbl_write_pickup.F`
  - L59 Called from outside the package — `BBL_TENDENCY_APPLY_S`, `BBL_TENDENCY_APPLY_T`, `BBL_CALC_RHO`, `BBL_CALC_RHS`, `BBL_DIAGNOSTICS_STATE`, `BBL_CHECK`, `BBL_INIT_FIXED`, `BBL_INIT_VARIA`, `BBL_READPARMS`, `BBL_WRITE_PICKUP`
  - L73 Verification experiments compiling it (1) — `global_oce_latlon`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/bbl@bbl.md`
- L1 pkg/bbl  (bbl: ~/Documents/research/ECCO/BBL/MITgcm) — `useBBL`
  - L11 README — `bbl_tracer`, `bbl_rho`, `bbl_eta`, `resTracer`
  - L45 Namelist parameters
    - L46 BBL_PARM01 — `bbl_wvel`, `bbl_hvel`, `bbl_initEta`, `bbl_useNofSpeed`, `bbl_nofCoeff`, `bbl_fMin`, `bbl_selectEntrain`, `bbl_selectDetrain`, `bbl_uDet`, `bbl_gpTaper`, `bbl_tauMix`, `bbl_coldStart`, `bbl_doPtracers`, `bbl_thetaFile`, `bbl_saltFile`, `bbl_etaFile`
  - L64 CPP options (defaults as shipped) — `BBL_CHECK_BUDGET`
  - L67 Headers — `BBL.h`, `BBL_OPTIONS.h`, `BBL_PTR.h`, `bbl_ad_check_lev1_dir.h`, `bbl_ad_check_lev2_dir.h`, `bbl_ad_check_lev3_dir.h`, `bbl_ad_check_lev4_dir.h`
  - L76 Routines (17) — `bbl_calc_rho.F`, `bbl_calc_rhs.F`, `bbl_check.F`, `bbl_diagnostics_init.F`, `bbl_diagnostics_state.F`, `bbl_init_fixed.F`, `bbl_init_varia.F`, `bbl_read_pickup.F`, `bbl_readparms.F`, `bbl_tendency_apply.F`, `bbl_write_pickup.F`
  - L79 Called from outside the package — `BBL_TENDENCY_APPLY_S`, `BBL_TENDENCY_APPLY_T`, `BBL_CALC_RHO`, `BBL_CALC_RHS`, `BBL_DIAGNOSTICS_STATE`, `BBL_CHECK`, `BBL_INIT_FIXED`, `BBL_INIT_VARIA`, `BBL_READPARMS`, `BBL_WRITE_PICKUP`, `BBL_TENDENCY_APPLY_PTR`
  - L94 Verification experiments compiling it (1)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/bbl@bbl_c68g.md`
- L1 pkg/bbl  (bbl_c68g: ~/Documents/research/ECCO/BBL/MITgcm_c68g) — `useBBL`
  - L11 README — `bbl_tracer`, `bbl_rho`, `bbl_eta`, `resTracer`
  - L45 Namelist parameters
    - L46 BBL_PARM01 — `bbl_wvel`, `bbl_hvel`, `bbl_initEta`, `bbl_useNofSpeed`, `bbl_nofCoeff`, `bbl_fMin`, `bbl_selectEntrain`, `bbl_selectDetrain`, `bbl_uDet`, `bbl_gpTaper`, `bbl_tauMix`, `bbl_coldStart`, `bbl_doPtracers`, `bbl_thetaFile`, `bbl_saltFile`, `bbl_etaFile`
  - L64 CPP options (defaults as shipped) — `BBL_CHECK_BUDGET`
  - L67 Headers — `BBL.h`, `BBL_OPTIONS.h`, `BBL_PTR.h`, `bbl_ad_check_lev1_dir.h`, `bbl_ad_check_lev2_dir.h`, `bbl_ad_check_lev3_dir.h`, `bbl_ad_check_lev4_dir.h`
  - L76 Routines (17) — `bbl_calc_rho.F`, `bbl_calc_rhs.F`, `bbl_check.F`, `bbl_diagnostics_init.F`, `bbl_diagnostics_state.F`, `bbl_init_fixed.F`, `bbl_init_varia.F`, `bbl_read_pickup.F`, `bbl_readparms.F`, `bbl_tendency_apply.F`, `bbl_write_pickup.F`
  - L79 Called from outside the package — `BBL_TENDENCY_APPLY_S`, `BBL_TENDENCY_APPLY_T`, `BBL_CALC_RHO`, `BBL_CALC_RHS`, `BBL_DIAGNOSTICS_STATE`, `BBL_CHECK`, `BBL_INIT_FIXED`, `BBL_INIT_VARIA`, `BBL_READPARMS`, `BBL_WRITE_PICKUP`, `BBL_TENDENCY_APPLY_PTR`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/bling.md`
- L1 pkg/bling — `useBLING`
  - L11 Namelist parameters
    - L12 ABIOTIC_PARMS — `selectBTconst`, `selectFTconst`, `selectHFconst`, `selectK1K2const`, `selectPHsolver`
    - L21 BIOTIC_PARMS — `CtoN`, `CtoP`, `NtoP`, `HtoC`, `NO3toN`, `O2toN`, `O2toP`, `CatoN`, `CatoP`, `masstoN`, `Pc_0_diaz`, `alpha_photo`, `gamma_DON`, `k_Fe_diaz`, `k_NO3`, `k_PtoN`, `k_FetoN`, `k_NO3_sm`, `k_NO3_lg`, `k_PO4_sm`, `k_PO4_lg`, `k_Fe_sm`, `k_Fe_lg`, `PtoN_min`, `PtoN_max`, `FetoN_min`, `FetoN_max`, `kappa_eppley_diaz`, `phi_dvm`, `sigma_dvm` (+46)
    - L106 BLING_FORCING — `bling_windFile`, `bling_atmospFile`, `bling_iceFile`, `bling_ironFile`, `bling_silicaFile`, `bling_psmFile`, `bling_plgFile`, `bling_PdiazFile`, `bling_forcingPeriod`, `bling_forcingCycle`, `bling_Pc_2dFile`, `bling_Pc_2d_diazFile`, `bling_k0_2dFile`, `bling_alpha_photo2dFile`, `bling_phi_DOM2dFile`, `bling_k_Fe2dFile`, `bling_k_Fe_diaz2dFile`, `bling_gamma_POM2dFile`, `bling_wsink0_2dFile`, `bling_phi_sm2dFile`, `bling_phi_lg2dFile`, `bling_pCO2`, `river_conc_po4`, `river_dom_to_nut`, `apco2RepCycle`, `apco2StartTime`, `exf_inscal_apco2`, `exf_outscal_apco2`, `apco2_exfremo_intercept`, `apco2_exfremo_slope` (+7)
  - L150 CPP options (defaults as shipped) — `USE_BLING_V1`, `USE_SIBLING`, `USE_BLING_DVM`, `ADVECT_PHYTO`, `BLING_NO_NEG`, `MIN_NUT_LIM`, `SIZE_DEP_LIM`, `ML_MEAN_LIGHT`, `ML_MEAN_PHYTO`, `BLING_USE_THRESHOLD_MLD`, `USE_QSW`, `PHYTO_SELF_SHADING`, `BLING_ADJOINT_SAFE`, `CARBONCHEM_SOLVESAPHE`, `CARBONCHEM_TOTALPHSCALE`, `NEW_FRAC_EXP`
  - L169 Headers — `BLING_LOAD.h`, `BLING_OPTIONS.h`, `BLING_VARS.h`, `bling_ad_check_lev1_dir.h`, `bling_ad_check_lev2_dir.h`, `bling_ad_check_lev3_dir.h`, `bling_ad_check_lev4_dir.h`
  - L178 Routines (19) — `bling_airseaflux.F`, `bling_bio.F`, `bling_bio_nitrogen.F`, `bling_carbonate_init.F`, `bling_carbonate_sys.F`, `bling_diagnostics_init.F`, `bling_fields_load.F`, `bling_ini_forcing.F`, `bling_init_fixed.F`, `bling_init_varia.F`, `bling_light.F`, `bling_main.F`, `bling_min_val.F`, `bling_mixedlayer.F`, `bling_read_pickup.F`, `bling_readparms.F`, `bling_sgs.F`, `bling_tr_register.F`, `bling_write_pickup.F`
  - L181 Called from outside the package — `BLING_FIELDS_LOAD`, `BLING_MAIN`, `BLING_INIT_FIXED`, `BLING_CARBONATE_INIT`, `BLING_INIT_VARIA`, `BLING_INI_FORCING`, `BLING_READPARMS`, `BLING_TR_REGISTER`, `BLING_WRITE_PICKUP`
  - L192 Verification experiments compiling it (1) — `global_oce_biogeo_bling`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/bulk_force.md`
- L1 pkg/bulk_force — `useBULK_FORCE`
  - L8 Namelist parameters
    - L9 BULKF_CONST — `rhoA`, `rhoFW`, `cpAir`, `cDrag_1`, `cDrag_2`, `cDrag_3`, `cStantonS`, `cStantonU`, `cDalton`, `humid_fac`, `saltQsFac`, `gamma_blk`, `atm_emissivity`, `ocean_emissivity`, `snow_emissivity`, `ice_emissivity`, `dTstab`, `ocean_albedo`
    - L43 BULKF_PARM01 — `useFluxFormula_AIM`, `blk_nIter`, `calcWindStress`, `blk_taveFreq`, `AirTempFile`, `AirHumidityFile`, `RainFile`, `SolarFile`, `LongwaveFile`, `UWindFile`, `VWindFile`, `RunoffFile`, `WSpeedFile`, `QnetFile`, `EmPFile`, `CloudFile`, `airPotTempFile`
    - L61 BULKF_PARM02 — `qnet_off`, `empmr_off`
  - L66 CPP options (defaults as shipped) — `CONSERV_BULKF`, `ALLOW_FORMULA_AIM`
  - L70 Headers — `BULKF.h`, `BULKF_CONSERV.h`, `BULKF_INT.h`, `BULKF_PARAMS.h`, `BULK_FORCE_OPTIONS.h`
  - L77 Routines (11) — `bulkf_ave.F`, `bulkf_fields_load.F`, `bulkf_flux_adjust.F`, `bulkf_forcing.F`, `bulkf_formula_aim.F`, `bulkf_formula_lanl.F`, `bulkf_formula_lay.F`, `bulkf_init_varia.F`, `bulkf_output.F`, `bulkf_readparms.F`, `bulkf_sh2rh_aim.F`
  - L80 Called from outside the package — `BULKF_FORCING`, `BULKF_FIELDS_LOAD`, `BULKF_INIT_VARIA`, `BULKF_READPARMS`, `BULKF_FORMULA_AIM`, `BULKF_FORMULA_LANL`, `BULKF_FORMULA_LAY`, `BULKF_FLUX_ADJUST`
  - L90 Verification experiments compiling it (1) — `global_ocean.cs32x15`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/cal.md`
- L1 pkg/cal — `useCAL`
  - L10 Namelist parameters
    - L11 CAL_NML — `TheCalendar`, `startDate_1`, `startDate_2`, `calendarDumps`
  - L17 Headers — `CAL_OPTIONS.h`, `cal.h`
  - L21 Routines (32) — `cal_addtime.F`, `cal_checkdate.F`, `cal_compdates.F`, `cal_convdate.F`, `cal_copydate.F`, `cal_daysformonth.F`, `cal_dayspermonth.F`, `cal_fulldate.F`, `cal_getdate.F`, `cal_getmonthsrec.F`, `cal_init_fixed.F`, `cal_intdays.F`, `cal_intmonths.F`, `cal_intyears.F`, `cal_isleap.F`, `cal_monthsforyear.F`, `cal_monthsperyear.F`, `cal_numints.F`, `cal_printdate.F`, `cal_printerror.F`, `cal_readparms.F`, `cal_set.F`, `cal_stepsforday.F`, `cal_stepsperday.F`, `cal_subdates.F`, `cal_summary.F`, `cal_time2dump.F`, `cal_timeinterval.F`, `cal_timepassed.F`, `cal_timestamp.F` (+2)
  - L24 Called from outside the package — `CAL_TIME2DUMP`, `CAL_GETDATE`, `CAL_INIT_FIXED`, `CAL_READPARMS`, `CAL_ADDTIME`, `CAL_COPYDATE`, `CAL_GETMONTHSREC`, `CAL_TIMEINTERVAL`, `CAL_TIMEPASSED`, `CAL_TOSECONDS`, `CAL_FULLDATE`, `CAL_CONVDATE`, `CAL_PRINTDATE`
  - L101 Verification experiments compiling it (2) — `global_oce_biogeo_bling`, `obcs_ctrl`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/cd_code.md`
- L1 pkg/cd_code — `useCD_CODE`
  - L10 CPP options (defaults as shipped) — `CD_CODE_NO_AB_MOMENTUM`, `CD_CODE_NO_AB_CORIOLIS`
  - L14 Headers — `CD_CODE_OPTIONS.h`, `CD_CODE_VARS.h`, `cd_code_ad_check_lev1_dir.h`, `cd_code_ad_check_lev2_dir.h`, `cd_code_ad_check_lev3_dir.h`, `cd_code_ad_check_lev4_dir.h`
  - L22 Routines (5) — `cd_code_ini_vars.F`, `cd_code_init_fixed.F`, `cd_code_read_pickup.F`, `cd_code_scheme.F`, `cd_code_write_pickup.F`
  - L25 Called from outside the package — `CD_CODE_INIT_FIXED`, `CD_CODE_INI_VARS`, `CD_CODE_WRITE_PICKUP`, `CD_CODE_SCHEME`
  - L31 Verification experiments compiling it (14) — `bottom_ctrl_5x5`, `cfc_example`, `global_oce_biogeo_bling`, `global_oce_latlon`, `global_ocean.90x40x15`, `ideal_2D_oce`, `lab_sea`, `so_box_biogeo`, `tutorial_global_oce_biogeo`, `tutorial_global_oce_latlon`, `tutorial_global_oce_optim`, `tutorial_tracer_adjsens`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/cfc.md`
- L1 pkg/cfc — `useCFC`
  - L11 Namelist parameters
    - L12 CFC_FORCING — `atmCFC_inpFile`, `atmCFC_recSepTime`, `atmCFC_timeOffset`, `atmCFC_yNorthBnd`, `atmCFC_ySouthBnd`, `CFC_windFile`, `CFC_atmospFile`, `CFC_iceFile`, `CFC_forcingPeriod`, `CFC_forcingCycle`
  - L24 Headers — `CFC.h`, `CFC_ATMOS.h`, `CFC_SIZE.h`, `cfc_ad_check_lev1_dir.h`, `cfc_ad_check_lev2_dir.h`, `cfc_ad_check_lev3_dir.h`, `cfc_ad_check_lev4_dir.h`
  - L33 Routines (10) — `cfc11_forcing.F`, `cfc11_surfforcing.F`, `cfc12_forcing.F`, `cfc12_surfforcing.F`, `cfc_atmos.F`, `cfc_check.F`, `cfc_fields_load.F`, `cfc_param.F`, `cfc_readparms.F`, `cfc_tr_register.F`
  - L36 Called from outside the package — `CFC11_FORCING`, `CFC12_FORCING`, `CFC_CHECK`, `CFC_FIELDS_LOAD`, `CFC_ATMOS`, `CFC_PARAM`, `CFC_READPARMS`, `CFC_TR_REGISTER`
  - L46 Verification experiments compiling it (2) — `cfc_example`, `tutorial_cfc_offline`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/cheapaml.md`
- L1 pkg/cheapaml — `useCHEAPAML`
  - L9 Namelist parameters
    - L10 CHEAPAML_CONST — `cheapaml_ntim`, `cheapaml_mask_width`, `cheapaml_h`, `cheapaml_kdiff`, `cheap_tauRelax`, `cheap_tauRelaxOce`, `cdrag_1`, `cdrag_2`, `cdrag_3`, `gasR`, `gamma_blk`, `humid_fac`, `cheap_pr1`, `cheap_pr2`, `cheapaml_taurelax`, `cheapaml_taurelaxocean`
    - L38 CHEAPAML_PARM01 — `periodicExternalForcing_cheap`, `externForcingPeriod_cheap`, `externForcingCycle_cheap`, `AirTempFile`, `SolarFile`, `UWindFile`, `VWindFile`, `TrFile`, `QrFile`, `AirQFile`, `UStressFile`, `VStressFile`, `WaveHFile`, `WavePFile`, `TracerFile`, `TracerRfile`, `cheapMaskFile`, `cheap_hFile`, `cheap_clFile`, `cheap_dlwFile`, `cheap_prFile`
    - L60 CHEAPAML_PARM02 — `cheapamlXperiodic`, `cheapamlYperiodic`, `useFreshWaterFlux`, `useFluxLimit`, `FluxFormula`, `WaveModel`, `useStressOption`, `useCheapTracer`, `useTimeVarBLH`, `useClouds`, `useDLongWave`, `usePrecip`, `useRelativeWind`
  - L75 CPP options (defaults as shipped) — `INCONSISTENT_WIND_LOCATION`
  - L78 Headers — `CHEAPAML.h`, `CHEAPAML_OPTIONS.h`
  - L82 Routines (16) — `cheapaml.F`, `cheapaml_calc_rhs.F`, `cheapaml_coare3_flux.F`, `cheapaml_copy_edges.F`, `cheapaml_diagnostics_init.F`, `cheapaml_fields_load.F`, `cheapaml_init_fixed.F`, `cheapaml_init_varia.F`, `cheapaml_lanl_flux.F`, `cheapaml_read_pickup.F`, `cheapaml_readparms.F`, `cheapaml_seaice.F`, `cheapaml_timestep.F`, `cheapaml_write_pickup.F`
  - L85 Called from outside the package — `CHEAPAML_FIELDS_LOAD`, `CHEAPAML_INIT_FIXED`, `CHEAPAML_INIT_VARIA`, `CHEAPAML_READPARMS`, `CHEAPAML_WRITE_PICKUP`
  - L93 Verification experiments compiling it (1) — `cheapAML_box`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/chronos.md`
- L1 pkg/chronos — `useCHRONOS`
  - L7 Headers — `chronos.h`
  - L10 Routines (16) — `chronos.F`
  - L13 Called from outside the package — `SET_ALARM`, `GET_TIME`, `GET_ALARM`, `INTERP_TIME`, `TIME_BOUND`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/compon_communic.md`
- L1 pkg/compon_communic — `useCOMPON_COMMUNIC`
  - L7 Headers — `CPLR_SIG.h`
  - L10 Routines (40) — `comprecv_i4vec.F`, `comprecv_r4.F`, `comprecv_r4tiles.F`, `comprecv_r8.F`, `comprecv_r8tiles.F`, `compsend_i4vec.F`, `compsend_r4.F`, `compsend_r4tiles.F`, `compsend_r8.F`, `compsend_r8tiles.F`, `couprecv_i4vec.F`, `couprecv_r4.F`, `couprecv_r4tiles.F`, `couprecv_r8.F`, `couprecv_r8tiles.F`, `coupsend_i4vec.F`, `coupsend_r4.F`, `coupsend_r4tiles.F`, `coupsend_r8.F`, `coupsend_r8tiles.F`, `generate_tag.F`, `mitcomponent_init.F`, `mitcomponent_register.F`, `mitcomponent_tile_register.F`, `mitcoupler_init.F`, `mitcoupler_register.F`, `mitcoupler_tile_register.F`, `mitcplr_all_check.F`, `mitcplr_char2dbl.F`, `mitcplr_char2int.F` (+10)
  - L13 Called from outside the package — `COUPRECV_R8TILES`, `COUPRECV_I4VEC`, `MITCOUPLER_TILE_REGISTER`, `COUPSEND_R8TILES`, `MITCPLR_ALL_CHECK`, `MITCOUPLER_INIT`, `COMPSEND_I4VEC`, `COMPSEND_R8TILES`, `COMPRECV_R8TILES`, `COMPRECV_I4VEC`, `MITCOMPONENT_INIT`, `MITCOMPONENT_TILE_REGISTER`, `COUPSEND_I4VEC`
  - L57 Verification experiments compiling it (1)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/cost.md`
- L1 pkg/cost — `useCOST`
  - L11 Namelist parameters
    - L12 COST_NML — `mult_atl`, `mult_test`, `mult_tracer`, `multTheta`, `multSalt`, `multUvel`, `multVvel`, `multEtan`, `mult_depth`, `mult_temp_tut`, `mult_hflux_tut`, `cost_mask_file`
  - L27 CPP options (defaults as shipped) — `ALLOW_COST_STATE_FINAL`, `ALLOW_COST_VECTOR`, `ALLOW_COST_ATLANTIC_HEAT`, `ALLOW_COST_ATLANTIC_HEAT_DOMASS`, `ALLOW_COST_TEST`, `ALLOW_COST_TSQUARED`, `ALLOW_COST_TRACER`, `ALLOW_COST_TEMP`, `ALLOW_COST_HFLUXM`, `ALLOW_DIC_COST`, `ALLOW_THSICE_COST_TEST`, `ALLOW_COST_SHELFICE`, `ALLOW_COST_STREAMICE`
  - L42 Headers — `COST_OPTIONS.h`, `COST_TAP_ADJ.h`, `adcost.h`, `cost.h`, `g_cost.h`
  - L49 Routines (18) — `cost_accumulate_mean.F`, `cost_atlantic_heat.F`, `cost_check.F`, `cost_copy_file.F`, `cost_dependent_init.F`, `cost_depth.F`, `cost_driver.F`, `cost_final.F`, `cost_final_restore.F`, `cost_final_store.F`, `cost_init_fixed.F`, `cost_init_varia.F`, `cost_readparms.F`, `cost_state_final.F`, `cost_test.F`, `cost_tile.F`, `cost_tracer.F`, `cost_vector.F`
  - L52 Called from outside the package — `COST_TILE`, `COST_INIT_VARIA`, `COST_CHECK`, `COST_INIT_FIXED`, `COST_READPARMS`, `COST_DRIVER`, `COST_FINAL`, `COST_DEPENDENT_INIT`, `COST_FINAL_RESTORE`, `COST_FINAL_STORE`
  - L65 Verification experiments compiling it (16) — `bottom_ctrl_5x5`, `global_oce_biogeo_bling`, `global_oce_latlon`, `global_ocean.90x40x15`, `global_ocean.cs32x15`, `halfpipe_streamice`, `lab_sea`, `obcs_ctrl`, `offline_exf_seaice`, `tutorial_dic_adjoffline`, `tutorial_global_oce_biogeo`, `tutorial_global_oce_optim`, `tutorial_tracer_adjsens`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/ctrl.md`
- L1 pkg/ctrl — `useCTRL`
  - L12 Namelist parameters
    - L13 CTRL_NML — `xx_theta_file`, `xx_salt_file`, `xx_hflux_file`, `xx_hflux_remo_intercept`, `xx_hflux_remo_slope`, `xx_hfluxstartdate1`, `xx_hfluxstartdate2`, `xx_hfluxperiod`, `xx_sflux_file`, `xx_sflux_remo_intercept`, `xx_sflux_remo_slope`, `xx_sfluxstartdate1`, `xx_sfluxstartdate2`, `xx_sfluxperiod`, `xx_tauu_file`, `xx_tauu_remo_intercept`, `xx_tauu_remo_slope`, `xx_tauustartdate1`, `xx_tauustartdate2`, `xx_tauuperiod`, `xx_tauv_file`, `xx_tauv_remo_intercept`, `xx_tauv_remo_slope`, `xx_tauvstartdate1`, `xx_tauvstartdate2`, `xx_tauvperiod`, `xx_atemp_file`, `xx_atemp_remo_intercept`, `xx_atemp_remo_slope`, `xx_atempstartdate1` (+154)
    - L199 CTRL_PACKNAMES
    - L209 CTRL_NML_GENARR — `xx_genarr2d_file`, `xx_genarr3d_file`, `xx_gentim2d_file`, `xx_genarr2d_weight`, `xx_genarr3d_weight`, `xx_gentim2d_weight`, `xx_genarr2d_bounds`, `xx_genarr3d_bounds`, `xx_gentim2d_bounds`, `genarr2dPrecond`, `genarr3dPrecond`, `gentim2dPrecond`, `mult_gentim2d`, `mult_genarr2d`, `mult_genarr3d`, `xx_genarr2d_preproc`, `xx_genarr2d_preproc_i`, `xx_genarr2d_preproc_r`, `xx_genarr2d_preproc_c`, `xx_genarr3d_preproc`, `xx_genarr3d_preproc_i`, `xx_genarr3d_preproc_r`, `xx_genarr3d_preproc_c`, `xx_gentim2d_preproc`, `xx_gentim2d_preproc_i`, `xx_gentim2d_preproc_r`, `xx_gentim2d_preproc_c`, `xx_gentim2d_period`, `xx_gentim2d_startdate1`, `xx_gentim2d_startdate2` (+2)
    - L242 OPTIM — `dfminFrac`
  - L257 CPP options (defaults as shipped) — `CTRL_SET_PREC_32`, `ALLOW_NONDIMENSIONAL_CONTROL_IO`, `EXCLUDE_CTRL_PACK`, `ALLOW_PACKUNPACK_METHOD2`, `CTRL_DO_PACK_UNPACK_ONLY`, `CTRL_PACK_PRECISE`, `CTRL_UNPACK_PRECISE`, `CTRL_DELZNORM`, `READ_OLD_CTRL_PACK_FILE`, `ALLOW_KAPGM_CONTROL`, `ALLOW_KAPREDI_CONTROL`, `ALLOW_DIFFKR_CONTROL`, `ALLOW_BOTTOMDRAG_CONTROL`, `ALLOW_DIC_CONTROL`, `ALLOW_DEPTH_CONTROL`, `USE_SMOOTH_MIN`, `ALLOW_GENARR2D_CONTROL`, `ALLOW_GENARR3D_CONTROL`, `ALLOW_GENTIM2D_CONTROL`, `ALLOW_UVEL0_CONTROL`, `ALLOW_VVEL0_CONTROL`, `ALLOW_OBCSN_CONTROL`, `ALLOW_OBCSS_CONTROL`, `ALLOW_OBCSE_CONTROL`, `ALLOW_OBCSW_CONTROL`, `ALLOW_OBCS_CONTROL`, `ALLOW_OBCS_CONTROL_MODES`, `ALLOW_OBCS_WEIGHTS2D`, `ALLOW_ADCTRLBOUND`, `ALLOW_ROTATE_UV_CONTROLS` (+3)
  - L292 Headers — `CTRL.h`, `CTRL_DUMMY.h`, `CTRL_FIELDS.h`, `CTRL_GENARR.h`, `CTRL_OBCS.h`, `CTRL_OPTIONS.h`, `CTRL_SIZE.h`, `OPTIMCYCLE.h`, `ctrl_ad_check_lev1_dir.h`, `ctrl_ad_check_lev2_dir.h`, `ctrl_ad_check_lev3_dir.h`, `ctrl_ad_check_lev4_dir.h`, `ctrl_local_params.h`
  - L307 Routines (71) — `adctrl_bound.F`, `ctrl_bound.F`, `ctrl_bound_tl.F`, `ctrl_check.F`, `ctrl_convert_header.F`, `ctrl_cost_driver.F`, `ctrl_cost_final.F`, `ctrl_cost_gen.F`, `ctrl_depth_ini.F`, `ctrl_get_gen.F`, `ctrl_get_gen_rec.F`, `ctrl_get_mask.F`, `ctrl_getobcse.F`, `ctrl_getobcsn.F`, `ctrl_getobcss.F`, `ctrl_getobcsw.F`, `ctrl_getrec.F`, `ctrl_hfacc_ini.F`, `ctrl_init.F`, `ctrl_init_ctrlvar.F`, `ctrl_init_fixed.F`, `ctrl_init_obcs_variables.F`, `ctrl_init_rec.F`, `ctrl_init_variables.F`, `ctrl_init_wet.F`, `ctrl_map_forcing.F`, `ctrl_map_genarr.F`, `ctrl_map_gentim2d.F`, `ctrl_map_ini.F`, `ctrl_map_ini_ecco.F` (+26)
  - L310 Called from outside the package — `CTRL_MAP_FORCING`, `CTRL_DEPTH_INI`, `CTRL_MAP_GENTIM2D`, `CTRL_CHECK`, `CTRL_INIT_FIXED`, `CTRL_INIT_VARIABLES`, `CTRL_READPARMS`, `OPTIM_READPARMS`, `CTRL_PACK`, `CTRL_UNPACK`, `CTRL_INIT`, `CTRL_SET_FNAME`, `CTRL_SET_UNPACK_XY`, `CTRL_SET_UNPACK_XYZ`, `CTRL_SET_UNPACK_XZ`, `CTRL_SET_UNPACK_YZ`, `CTRL_SET_PACK_XY`, `CTRL_SET_PACK_XYZ`, `CTRL_SET_PACK_XZ`, `CTRL_SET_PACK_YZ`, `CTRL_COST_DRIVER`, `CTRL_COST_FINAL`, `CTRL_GETOBCSE`, `CTRL_GETOBCSN`, `CTRL_GETOBCSS`, `CTRL_GETOBCSW`
  - L341 Verification experiments compiling it (16) — `bottom_ctrl_5x5`, `global_oce_biogeo_bling`, `global_oce_latlon`, `global_ocean.90x40x15`, `global_ocean.cs32x15`, `halfpipe_streamice`, `lab_sea`, `obcs_ctrl`, `offline_exf_seaice`, `tutorial_dic_adjoffline`, `tutorial_global_oce_biogeo`, `tutorial_global_oce_optim`, `tutorial_tracer_adjsens`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/darwin.md`
- L1 pkg/darwin  (darwin3: ~/Documents/GitHub/darwin3) — `useDARWIN`
  - L8 Namelist parameters
    - L9 DARWIN_INDICES — `i_DIC`, `i_NO3`, `i_NO2`, `i_NH4`, `i_PO4`, `i_FeT`, `i_SiO2`, `i_DOC`, `i_DON`, `i_DOP`, `i_DOFe`, `i_POC`, `i_PON`, `i_POP`, `i_POFe`, `i_POSi`, `i_PIC`, `i_ALK`, `i_O2`, `i_CDOM`, `e_CDOM`, `i_fe`, `e_fe`, `i_si`, `e_si`, `i_Chl`, `e_Chl`, `i_ch`, `e_ch`, `n_DARWIN`
    - L46 DARWIN_FORCING_PARAMS — `darwin_useEXFwind`, `darwin_useQsw`, `darwin_useSEAICE`, `R_DOFe_DOP_runoff`, `R_NO3_DIN_runoff`, `R_NO2_DIN_runoff`, `R_NH4_DIN_runoff`, `R_DFe_DIP_runoff`, `R_DIP_IP_runoff`, `R_POFe_POP_runoff`, `R_ALK_DIC_runoff`, `PARRepCycle`, `PARStartTime`, `PAR_exfremo_intercept`, `PAR_exfremo_slope`, `darwin_inscal_PAR`, `ironRepCycle`, `ironStartTime`, `iron_exfremo_intercept`, `iron_exfremo_slope`, `darwin_inscal_iron`, `iceRepCycle`, `iceStartTime`, `ice_exfremo_intercept`, `ice_exfremo_slope`, `darwin_inscal_ice`, `windRepCycle`, `windStartTime`, `wind_exfremo_intercept`, `wind_exfremo_slope` (+73)
    - L234 DARWIN_INTERP_PARAMS — `PAR_lon0`, `PAR_lat0`, `PAR_nlon`, `PAR_nlat`, `PAR_lon_inc`, `PAR_interpMethod`, `PAR_lat_inc`, `iron_lon0`, `iron_lat0`, `iron_nlon`, `iron_nlat`, `iron_lon_inc`, `iron_interpMethod`, `iron_lat_inc`, `ice_lon0`, `ice_lat0`, `ice_nlon`, `ice_nlat`, `ice_lon_inc`, `ice_interpMethod`, `ice_lat_inc`, `wind_lon0`, `wind_lat0`, `wind_nlon`, `wind_nlat`, `wind_lon_inc`, `wind_interpMethod`, `wind_lat_inc`, `pCO2_lon0`, `pCO2_lat0` (+82)
    - L347 DARWIN_CONSTANTS
    - L349 CARBON_CONSTANTS
    - L371 DARWIN_PARAMS — `darwin_seed`, `darwin_strict_check`, `iDEBUG`, `jDEBUG`, `kDEBUG`, `darwin_pickupSuff`, `darwin_linFSConserve`, `darwin_read_phos`, `darwin_chlInitBalanced`, `darwin_chlIter0`, `katten_w`, `katten_chl`, `TempAeArr`, `TemprefArr`, `TempCoeffArr`, `reminTempAe`, `mortTempAe`, `mort2TempAe`, `uptakeTempAe`, `ligand_tot`, `ligand_stab`, `scav_rat`, `scav_inter`, `scav_exp`, `scav_R_POPPOC`, `scav_tau`, `scav_POC_wgt`, `scav_PSi_wgt`, `scav_PIC_wgt`, `scav_degrPOM` (+60)
    - L494 DARWIN_CDOM_PARAMS — `fracCDOM`, `R_NC_CDOM`, `R_FeC_CDOM`, `R_PC_CDOM`, `R_NP_CDOM`, `R_FeP_CDOM`, `R_CP_CDOM`
    - L506 DARWIN_DEPENDENT_PARAMS — `laCDOM`, `darwin_dependent_i_dummy`
    - L509 DARWIN_RADTRANS_PARAMS — `darwin_waterAbsorbFile`, `darwin_phytoAbsorbFile`, `darwin_particleAbsorbFile`, `darwin_part_size_P`, `darwin_bbmin`, `darwin_bbw`, `darwin_lambda_aCDOM`, `darwin_Sdom`, `darwin_aCDOM_fac`, `darwin_rCDOM`, `darwin_RPOC`, `darwin_allomSpectra`, `darwin_aCarCell`, `darwin_bCarCell`, `darwin_absorpSlope`, `darwin_bbbSlope`, `darwin_scatSwitchSizeLog`, `darwin_scatSlopeSmall`, `darwin_scatSlopeLarge`
    - L529 DARWIN_RT_DEPPARAMS — `aphy_chl_type`, `aphy_chl_ps_type`, `aphy_mgC_type`, `bphy_mgC_type`, `bbphy_mgC_type`, `asize_mgC`, `apart_P`, `bpart_P`, `bbpart_P`, `exCDOM`
    - L549 DARWIN_RANDOM_PARAMS — `diaz_growfac`, `cocco_growfac`, `diatom_growfac`, `tempcoeff2_small`, `tempcoeff2_big`, `val_R_NC`, `val_R_NC_diaz`, `val_R_PC`, `val_R_SiC_diatom`, `val_R_FeC`, `val_R_FeC_diaz`, `val_R_PICPOC`, `val_R_ChlC`, `val_R_NC_zoo`, `val_R_PC_zoo`, `val_R_SiC_zoo`, `val_R_FeC_zoo`, `val_R_PICPOC_zoo`, `val_R_ChlC_zoo`, `SmallSink`, `BigSink`, `SmallPsat`, `BigPsat`, `ProcPsat`, `UniDzPsat`, `CoccoPsat`, `SmallPsatrange`, `BigPsatrange`, `ProcPsatrange`, `UniDzPsatrange` (+31)
    - L649 DARWIN_TRAIT_PARAMS — `darwin_sort_biovol`, `darwin_effective_ksat`, `darwin_select_kn_allom`, `grp_logvolind`, `grp_biovol`, `grp_names`, `grp_nplank`, `grp_photo`, `grp_bacttype`, `grp_aerobic`, `grp_denit`, `grp_pred`, `grp_prey`, `grp_hasSi`, `grp_hasPIC`, `grp_diazo`, `grp_useNH4`, `grp_useNO2`, `grp_useNO3`, `grp_combNO`, `grp_aptype`, `grp_tempMort`, `grp_tempMort2`, `grp_tempGraz`, `a_Xmin`, `a_R_NC`, `a_R_PC`, `a_R_SiC`, `a_R_FeC`, `a_R_ChlC` (+117)
    - L802 DARWIN_TRAITS — `isPhoto`, `bactType`, `isAerobic`, `isDenit`, `hasSi`, `hasPIC`, `useNH4`, `useNO2`, `useNO3`, `combNO`, `isPrey`, `isPred`, `tempMort`, `tempMort2`, `tempGraz`, `ExportFracMort`, `ExportFracMort2`, `ExportFracExude`, `FracExudeC`, `phytoTempCoeff`, `pkg_darwin_temperature_params`, `phytoTempExp1`, `phytoTempAe`, `phytoTempExp2`, `phytoTempOptimum`, `phytoDecayPower`, `hetTempAe`, `hetTempExp2`, `hetTempOptimum`, `hetDecayPower` (+39)
    - L899 DARWIN_DEPENDENT_TRAITS — `normI`, `biovol_bygroup`, `alpha_mean`
    - L909 DARWIN_RADTRANS_TRAITS — `aphy_chl`, `aphy_chl_ps`, `aphy_mgC`, `bphy_mgC`, `bbphy_mgC`
    - L915 DARWIN_RT_DEPTRAITS
  - L918 CPP options (defaults as shipped) — `DARWIN_ALLOW_NQUOTA`, `DARWIN_ALLOW_PQUOTA`, `DARWIN_ALLOW_FEQUOTA`, `DARWIN_ALLOW_SIQUOTA`, `DARWIN_ALLOW_CHLQUOTA`, `DARWIN_ALLOW_CSTORE`, `DARWIN_ALLOW_CSTORE_DIAGS`, `DARWIN_ALLOW_CDOM`, `DARWIN_ALLOW_CARBON`, `DARWIN_TOTALPHSCALE`, `DARWIN_SOLVESAPHE`, `DARWIN_USE_PLOAD`, `DARWIN_ALLOW_RADI`, `DARWIN_ALLOW_DENIT`, `DARWIN_ALLOW_EXUDE`, `ALLOW_OLD_VIRTUALFLUX`, `DARWIN_NITRATE_FELIMIT`, `DARWIN_BOTTOM_SINK`, `DARWIN_CDOM_UNITS_CARBON`, `DARWIN_NUTRIENT_RUNOFF`, `DARWIN_AVPAR`, `DARWIN_ALLOW_GEIDER`, `DARWIN_GEIDER_RHO_SYNTH`, `DARWIN_CHL_INIT_LEGACY`, `DARWIN_SCATTER_CHL`, `DARWIN_DIAG_IOP`, `DARWIN_GRAZING_SWITCH`, `DARWIN_ALLOMETRIC_PALAT`, `DARWIN_NOZOOTEMP`, `DARWIN_NOTEMP` (+17)
  - L967 Headers — `DARWIN_DIAGS.h`, `DARWIN_EXF_FIELDS.h`, `DARWIN_EXF_PARAMS.h`, `DARWIN_FIELDS.h`, `DARWIN_FLUX.h`, `DARWIN_INDICES.h`, `DARWIN_OPTIONS.h`, `DARWIN_PARAMS.h`, `DARWIN_RADTRANS.h`, `DARWIN_SIZE.h`, `DARWIN_TRAITPARAMS.h`, `DARWIN_TRAITS.h`
  - L981 Routines (47) — `darwin_add_botforc.F`, `darwin_add_surfforc.F`, `darwin_atmos.F`, `darwin_carbon_chem.F`, `darwin_check.F`, `darwin_cons.F`, `darwin_cons_reset.F`, `darwin_conserve_surface.F`, `darwin_diagnostics_init.F`, `darwin_diags.F`, `darwin_exf_init_fixed.F`, `darwin_exf_init_varia.F`, `darwin_exf_load.F`, `darwin_exf_readparms.F`, `darwin_fe_chem.F`, `darwin_fields_load.F`, `darwin_forcing.F`, `darwin_generate_allometric.F`, `darwin_generate_random.F`, `darwin_init_chl.F`, `darwin_init_fixed.F`, `darwin_init_varia.F`, `darwin_insol.F`, `darwin_light.F`, `darwin_light_radtrans.F`, `darwin_monitor.F`, `darwin_nut_supply.F`, `darwin_plankton.F`, `darwin_random.F`, `darwin_read_params.F` (+11)
  - L984 Called from outside the package — `DARWIN_CHECK`, `DARWIN_CONS`, `DARWIN_FIELDS_LOAD`, `DARWIN_ATMOS`, `DARWIN_CONSERVE_SURFACE`, `DARWIN_CONS_RESET`, `DARWIN_FORCING`, `DARWIN_NUT_SUPPLY`, `DARWIN_INIT_FIXED`, `DARWIN_INIT_VARIA`, `DARWIN_DIAGS`, `DARWIN_READPARMS`, `DARWIN_TR_REGISTER`, `DARWIN_WRITE_PICKUP`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/darwin@d3backport.md`
- L1 pkg/darwin  (d3backport: ~/Documents/research/debug/darwin3) — `useDARWIN`
  - L10 Namelist parameters
    - L11 DARWIN_INDICES — `i_DIC`, `i_NO3`, `i_NO2`, `i_NH4`, `i_PO4`, `i_FeT`, `i_SiO2`, `i_DOC`, `i_DON`, `i_DOP`, `i_DOFe`, `i_POC`, `i_PON`, `i_POP`, `i_POFe`, `i_POSi`, `i_PIC`, `i_ALK`, `i_O2`, `i_CDOM`, `e_CDOM`, `i_fe`, `e_fe`, `i_si`, `e_si`, `i_Chl`, `e_Chl`, `i_ch`, `e_ch`, `n_DARWIN`
    - L48 DARWIN_FORCING_PARAMS — `darwin_useEXFwind`, `darwin_useQsw`, `darwin_useSEAICE`, `R_DOFe_DOP_runoff`, `R_NO3_DIN_runoff`, `R_NO2_DIN_runoff`, `R_NH4_DIN_runoff`, `R_DFe_DIP_runoff`, `R_DIP_IP_runoff`, `R_POFe_POP_runoff`, `R_ALK_DIC_runoff`, `PARRepCycle`, `PARStartTime`, `PAR_exfremo_intercept`, `PAR_exfremo_slope`, `darwin_inscal_PAR`, `ironRepCycle`, `ironStartTime`, `iron_exfremo_intercept`, `iron_exfremo_slope`, `darwin_inscal_iron`, `iceRepCycle`, `iceStartTime`, `ice_exfremo_intercept`, `ice_exfremo_slope`, `darwin_inscal_ice`, `windRepCycle`, `windStartTime`, `wind_exfremo_intercept`, `wind_exfremo_slope` (+73)
    - L236 DARWIN_INTERP_PARAMS — `PAR_lon0`, `PAR_lat0`, `PAR_nlon`, `PAR_nlat`, `PAR_lon_inc`, `PAR_interpMethod`, `PAR_lat_inc`, `iron_lon0`, `iron_lat0`, `iron_nlon`, `iron_nlat`, `iron_lon_inc`, `iron_interpMethod`, `iron_lat_inc`, `ice_lon0`, `ice_lat0`, `ice_nlon`, `ice_nlat`, `ice_lon_inc`, `ice_interpMethod`, `ice_lat_inc`, `wind_lon0`, `wind_lat0`, `wind_nlon`, `wind_nlat`, `wind_lon_inc`, `wind_interpMethod`, `wind_lat_inc`, `pCO2_lon0`, `pCO2_lat0` (+82)
    - L349 DARWIN_CONSTANTS
    - L351 CARBON_CONSTANTS
    - L373 DARWIN_PARAMS — `darwin_seed`, `darwin_strict_check`, `iDEBUG`, `jDEBUG`, `kDEBUG`, `darwin_pickupSuff`, `darwin_linFSConserve`, `darwin_read_phos`, `darwin_chlInitBalanced`, `darwin_chlIter0`, `katten_w`, `katten_chl`, `TempAeArr`, `TemprefArr`, `TempCoeffArr`, `reminTempAe`, `mortTempAe`, `mort2TempAe`, `uptakeTempAe`, `ligand_tot`, `ligand_stab`, `scav_rat`, `scav_inter`, `scav_exp`, `scav_R_POPPOC`, `scav_tau`, `scav_POC_wgt`, `scav_PSi_wgt`, `scav_PIC_wgt`, `scav_degrPOM` (+92)
    - L528 DARWIN_CDOM_PARAMS — `fracCDOM`, `R_NC_CDOM`, `R_FeC_CDOM`, `R_PC_CDOM`, `R_NP_CDOM`, `R_FeP_CDOM`, `R_CP_CDOM`
    - L540 DARWIN_DEPENDENT_PARAMS — `laCDOM`, `darwin_dependent_i_dummy`
    - L543 DARWIN_RADTRANS_PARAMS — `darwin_waterAbsorbFile`, `darwin_phytoAbsorbFile`, `darwin_particleAbsorbFile`, `darwin_part_size_P`, `darwin_bbmin`, `darwin_bbw`, `darwin_lambda_aCDOM`, `darwin_Sdom`, `darwin_aCDOM_fac`, `darwin_rCDOM`, `darwin_RPOC`, `darwin_allomSpectra`, `darwin_aCarCell`, `darwin_bCarCell`, `darwin_absorpSlope`, `darwin_bbbSlope`, `darwin_scatSwitchSizeLog`, `darwin_scatSlopeSmall`, `darwin_scatSlopeLarge`
    - L563 DARWIN_RT_DEPPARAMS — `aphy_chl_type`, `aphy_chl_ps_type`, `aphy_mgC_type`, `bphy_mgC_type`, `bbphy_mgC_type`, `asize_mgC`, `apart_P`, `bpart_P`, `bbpart_P`, `exCDOM`
    - L583 DARWIN_RANDOM_PARAMS — `diaz_growfac`, `cocco_growfac`, `diatom_growfac`, `tempcoeff2_small`, `tempcoeff2_big`, `val_R_NC`, `val_R_NC_diaz`, `val_R_PC`, `val_R_SiC_diatom`, `val_R_FeC`, `val_R_FeC_diaz`, `val_R_PICPOC`, `val_R_ChlC`, `val_R_NC_zoo`, `val_R_PC_zoo`, `val_R_SiC_zoo`, `val_R_FeC_zoo`, `val_R_PICPOC_zoo`, `val_R_ChlC_zoo`, `SmallSink`, `BigSink`, `SmallPsat`, `BigPsat`, `ProcPsat`, `UniDzPsat`, `CoccoPsat`, `SmallPsatrange`, `BigPsatrange`, `ProcPsatrange`, `UniDzPsatrange` (+31)
    - L683 DARWIN_TRAIT_PARAMS — `darwin_sort_biovol`, `darwin_effective_ksat`, `darwin_select_kn_allom`, `grp_logvolind`, `grp_biovol`, `grp_names`, `grp_nplank`, `grp_photo`, `grp_bacttype`, `grp_aerobic`, `grp_denit`, `grp_pred`, `grp_prey`, `grp_hasSi`, `grp_hasPIC`, `grp_diazo`, `grp_useNH4`, `grp_useNO2`, `grp_useNO3`, `grp_combNO`, `grp_aptype`, `grp_tempMort`, `grp_tempMort2`, `grp_tempGraz`, `a_Xmin`, `a_R_NC`, `a_R_PC`, `a_R_SiC`, `a_R_FeC`, `a_R_ChlC` (+117)
    - L836 DARWIN_TRAITS — `isPhoto`, `bactType`, `isAerobic`, `isDenit`, `hasSi`, `hasPIC`, `useNH4`, `useNO2`, `useNO3`, `combNO`, `isPrey`, `isPred`, `tempMort`, `tempMort2`, `tempGraz`, `ExportFracMort`, `ExportFracMort2`, `ExportFracExude`, `FracExudeC`, `phytoTempCoeff`, `pkg_darwin_temperature_params`, `phytoTempExp1`, `phytoTempAe`, `phytoTempExp2`, `phytoTempOptimum`, `phytoDecayPower`, `hetTempAe`, `hetTempExp2`, `hetTempOptimum`, `hetDecayPower` (+39)
    - L933 DARWIN_DEPENDENT_TRAITS — `normI`, `biovol_bygroup`, `alpha_mean`
    - L943 DARWIN_RADTRANS_TRAITS — `aphy_chl`, `aphy_chl_ps`, `aphy_mgC`, `bphy_mgC`, `bbphy_mgC`
    - L949 DARWIN_RT_DEPTRAITS
  - L952 CPP options (defaults as shipped) — `DARWIN_ALLOW_NQUOTA`, `DARWIN_ALLOW_PQUOTA`, `DARWIN_ALLOW_FEQUOTA`, `DARWIN_ALLOW_SIQUOTA`, `DARWIN_ALLOW_CHLQUOTA`, `DARWIN_ALLOW_CSTORE`, `DARWIN_ALLOW_CSTORE_DIAGS`, `DARWIN_ALLOW_CDOM`, `DARWIN_ALLOW_CARBON`, `DARWIN_TOTALPHSCALE`, `DARWIN_SOLVESAPHE`, `DARWIN_USE_PLOAD`, `DARWIN_ALLOW_RADIv1`, `DARWIN_ALLOW_RADIv2`, `DARWIN_ALLOW_DENIT`, `DARWIN_ALLOW_EXUDE`, `ALLOW_OLD_VIRTUALFLUX`, `DARWIN_NITRATE_FELIMIT`, `DARWIN_BOTTOM_SINK`, `DARWIN_CDOM_UNITS_CARBON`, `DARWIN_NUTRIENT_RUNOFF`, `DARWIN_AVPAR`, `DARWIN_ALLOW_GEIDER`, `DARWIN_GEIDER_RHO_SYNTH`, `DARWIN_CHL_INIT_LEGACY`, `DARWIN_SCATTER_CHL`, `DARWIN_DIAG_IOP`, `DARWIN_GRAZING_SWITCH`, `DARWIN_ALLOMETRIC_PALAT`, `DARWIN_NOZOOTEMP` (+18)
  - L1002 Headers — `DARWIN_DIAGS.h`, `DARWIN_EXF_FIELDS.h`, `DARWIN_EXF_PARAMS.h`, `DARWIN_FIELDS.h`, `DARWIN_FLUX.h`, `DARWIN_INDICES.h`, `DARWIN_OPTIONS.h`, `DARWIN_PARAMS.h`, `DARWIN_RADTRANS.h`, `DARWIN_SIZE.h`, `DARWIN_TRAITPARAMS.h`, `DARWIN_TRAITS.h`
  - L1016 Routines (47) — `darwin_add_botforc.F`, `darwin_add_surfforc.F`, `darwin_atmos.F`, `darwin_carbon_chem.F`, `darwin_check.F`, `darwin_cons.F`, `darwin_cons_reset.F`, `darwin_conserve_surface.F`, `darwin_diagnostics_init.F`, `darwin_diags.F`, `darwin_exf_init_fixed.F`, `darwin_exf_init_varia.F`, `darwin_exf_load.F`, `darwin_exf_readparms.F`, `darwin_fe_chem.F`, `darwin_fields_load.F`, `darwin_forcing.F`, `darwin_generate_allometric.F`, `darwin_generate_random.F`, `darwin_init_chl.F`, `darwin_init_fixed.F`, `darwin_init_varia.F`, `darwin_insol.F`, `darwin_light.F`, `darwin_light_radtrans.F`, `darwin_monitor.F`, `darwin_nut_supply.F`, `darwin_plankton.F`, `darwin_random.F`, `darwin_read_params.F` (+11)
  - L1019 Called from outside the package — `DARWIN_CHECK`, `DARWIN_CONS`, `DARWIN_FIELDS_LOAD`, `DARWIN_ATMOS`, `DARWIN_CONSERVE_SURFACE`, `DARWIN_CONS_RESET`, `DARWIN_FORCING`, `DARWIN_NUT_SUPPLY`, `DARWIN_INIT_FIXED`, `DARWIN_INIT_VARIA`, `DARWIN_DIAGS`, `DARWIN_READPARMS`, `DARWIN_TR_REGISTER`, `DARWIN_WRITE_PICKUP`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/debug.md`
- L1 pkg/debug
  - L9 Headers — `DEBUG_OPTIONS.h`
  - L12 Routines (13) — `chksum_tiled.F`, `debug_call.F`, `debug_cs_corner_uv.F`, `debug_enter.F`, `debug_fld_stats_rl.F`, `debug_fld_stats_rs.F`, `debug_leave.F`, `debug_msg.F`, `debug_stats_rl.F`, `debug_stats_rs.F`, `fill_in_corners_rl.F`, `write_fullarray_rl.F`, `write_fullarray_rs.F`
  - L15 Called from outside the package — `DEBUG_ENTER`, `DEBUG_LEAVE`, `DEBUG_CALL`, `DEBUG_MSG`, `DEBUG_STATS_RL`, `DEBUG_CS_CORNER_UV`, `WRITE_FULLARRAY_RS`
  - L180 Verification experiments compiling it (58) — `advect_cs`, `advect_xz`, `aim.5l_Equatorial_Channel`, `aim.5l_LatLon`, `aim.5l_cs`, `atm_gray`, `bottom_ctrl_5x5`, `cfc_example`, `cheapAML_box`, `deep_anelastic`, `front_relax`, `global_oce_biogeo_bling`, `global_oce_latlon`, `global_ocean.90x40x15`, `global_ocean.cs32x15`, `halfpipe_streamice`, `ideal_2D_oce`, `internal_wave`, `inverted_barometer`, `lab_sea`, `matrix_example`, `obcs_ctrl`, `offline_exf_seaice`, `seaice_itd`, `seaice_obcs`, `shelfice_2d_remesh`, `short_surf_wave`, `so_box_biogeo`, `tutorial_advection_in_gyre`, `tutorial_baroclinic_gyre` (+12)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/diagnostics.md`
- L1 pkg/diagnostics — `useDIAGNOSTICS`
  - L10 Namelist parameters
    - L11 DIAGNOSTICS_LIST — `timePhase`, `averagingFreq`, `averagingPhase`, `repeatCycle`, `missing_value`, `missing_value_int`, `fileName`, `fileFlags`, `dumpAtLast`, `diag_dBugLevel`, `diag_mnc`, `useMissingValue`, `diagCG_maxIters`, `diagCG_resTarget`, `diagCG_pcOffDFac`, `diagCG_prtResFrq`, `xPsi0`, `yPsi0`, `diag_pickup_read`, `diag_pickup_write`, `diag_pickup_read_mnc`, `diag_pickup_write_mnc`, `diagMdsDir`, `diagMdsDirCreate`
    - L39 DIAG_STATIS_PARMS — `stat_freq`, `stat_phase`, `stat_region`, `stat_fields`, `stat_fName`, `diagSt_mnc`, `set_regMask`, `val_regMask`, `diagSt_regMaskFile`, `nSetRegMskFile`
  - L51 CPP options (defaults as shipped) — `DIAGSTATS_REGION_MASK`, `DIAGNOSTICS_HAS_PICKUP`
  - L55 Headers — `DIAGNOSTICS.h`, `DIAGNOSTICS_CALC.h`, `DIAGNOSTICS_P2SHARE.h`, `DIAGNOSTICS_SIZE.h`, `DIAGSTATS_REGIONS.h`, `DIAG_OPTIONS.h`
  - L63 Routines (65) — `diag_calc_psivel.F`, `diag_cg2d.F`, `diag_vegtile_fill.F`, `diagnostics_addtolist.F`, `diagnostics_calc_phivel.F`, `diagnostics_check.F`, `diagnostics_clear.F`, `diagnostics_fill.F`, `diagnostics_fill_field.F`, `diagnostics_fill_rs.F`, `diagnostics_fill_state.F`, `diagnostics_fract_fill.F`, `diagnostics_ini_io.F`, `diagnostics_init_early.F`, `diagnostics_init_fixed.F`, `diagnostics_init_varia.F`, `diagnostics_interp_p2p.F`, `diagnostics_interp_vert.F`, `diagnostics_is_on.F`, `diagnostics_list_check.F`, `diagnostics_main_init.F`, `diagnostics_mnc_out.F`, `diagnostics_out.F`, `diagnostics_read_pickup.F`, `diagnostics_readparms.F`, `diagnostics_scale_fill.F`, `diagnostics_scale_fill_rs.F`, `diagnostics_set_calc.F`, `diagnostics_set_levels.F`, `diagnostics_set_pointers.F` (+23)
  - L66 Called from outside the package — `DIAGNOSTICS_SCALE_FILL`, `DIAGNOSTICS_FILL`, `DIAGNOSTICS_COUNT`, `DIAGNOSTICS_FILL_RS`, `DIAGNOSTICS_SCALE_FILL_RS`, `DIAGNOSTICS_FILL_STATE`, `DIAGNOSTICS_WRITE`, `DIAGNOSTICS_SWITCH_ONOFF`, `DIAGNOSTICS_INI_IO`, `DIAGNOSTICS_CHECK`, `DIAGNOSTICS_INIT_EARLY`, `DIAGNOSTICS_INIT_FIXED`, `DIAGNOSTICS_MAIN_INIT`, `DIAGNOSTICS_INIT_VARIA`, `DIAGNOSTICS_READPARMS`, `DIAGNOSTICS_WRITE_PICKUP`, `DIAGNOSTICS_ADDTOLIST`, `DIAGNOSTICS_FRACT_FILL`, `DIAGNOSTICS_SETKLEV`, `DIAGNOSTICS_WRITE_ADJ`, `DIAG_VEGTILE_FILL`
  - L296 Verification experiments compiling it (41) — `advect_cs`, `advect_xz`, `aim.5l_cs`, `atm_gray`, `cfc_example`, `cheapAML_box`, `deep_anelastic`, `front_relax`, `global_oce_biogeo_bling`, `global_oce_latlon`, `global_ocean.90x40x15`, `global_ocean.cs32x15`, `halfpipe_streamice`, `ideal_2D_oce`, `inverted_barometer`, `lab_sea`, `obcs_ctrl`, `offline_exf_seaice`, `seaice_itd`, `shelfice_2d_remesh`, `short_surf_wave`, `so_box_biogeo`, `tutorial_advection_in_gyre`, `tutorial_baroclinic_gyre`, `tutorial_deep_convection`, `tutorial_global_oce_biogeo`, `tutorial_held_suarez_cs`, `tutorial_reentrant_channel`, `tutorial_rotating_tank`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/dic.md`
- L1 pkg/dic — `useDIC`
  - L11 Namelist parameters
    - L12 ABIOTIC_PARMS — `selectBTconst`, `selectFTconst`, `selectHFconst`, `selectK1K2const`, `selectPHsolver`, `useCalciteSaturation`, `calcOmegaCalciteFreq`, `nIterCO3`, `selectCalciteDissolution`, `WsinkPIC`, `selectCalciteBottomRemin`, `calciteDissolRate`, `calciteDissolExp`
    - L29 BIOTIC_PARMS — `R_OP`, `R_CP`, `R_NP`, `R_FeP`, `fesedflux_pcm`, `FeIntSec`, `ligand_stab`, `ligand_tot`, `alphaUniform`, `rainRatioUniform`
    - L54 DIC_FORCING — `DIC_windFile`, `DIC_atmospFile`, `DIC_silicaFile`, `DIC_deepSilicaFile`, `DIC_iceFile`, `DIC_parFile`, `DIC_chlaFile`, `DIC_ironFile`, `DIC_forcingPeriod`, `DIC_forcingCycle`, `dic_int1`, `dic_int2`, `dic_int3`, `dic_int4`, `dic_pCO2`
  - L71 CPP options (defaults as shipped) — `CARBONCHEM_SOLVESAPHE`, `CARBONCHEM_TOTALPHSCALE`, `DIC_BIOTIC`, `ALLOW_O2`, `ALLOW_FE`, `READ_PAR`, `DIC_NO_NEG`, `DIC_BOUNDS`, `USE_QSW`, `USE_QSW_UNDERICE`, `USE_PLOAD`, `ALLOW_OLD_VIRTUALFLUX`, `WATERVAP_BUG`, `DIC_CALCITE_SAT`, `LIGHT_CHL`, `DIC_AD_SAFE`
  - L91 Headers — `DIC_ATMOS.h`, `DIC_COST.h`, `DIC_CTRL.h`, `DIC_LOAD.h`, `DIC_OPTIONS.h`, `DIC_VARS.h`, `dic_ad_check_lev1_dir.h`, `dic_ad_check_lev2_dir.h`, `dic_ad_check_lev3_dir.h`, `dic_ad_check_lev4_dir.h`
  - L103 Routines (44) — `alk_surfforcing.F`, `bio_export.F`, `calcite_saturation.F`, `car_flux.F`, `car_flux_omega_top.F`, `carbon_chem.F`, `dic_atmos.F`, `dic_biotic_diags.F`, `dic_biotic_forcing.F`, `dic_biotic_init.F`, `dic_cost.F`, `dic_diagnostics_init.F`, `dic_fields_load.F`, `dic_fields_update.F`, `dic_ini_atmos.F`, `dic_ini_forcing.F`, `dic_init_fixed.F`, `dic_init_varia.F`, `dic_mnc_init.F`, `dic_read_co2_pickup.F`, `dic_read_pickup.F`, `dic_readparms.F`, `dic_set_control.F`, `dic_solvesaphe.F`, `dic_store_fluxco2.F`, `dic_surfforcing.F`, `dic_surfforcing_init.F`, `dic_tr_register.F`, `dic_write_pickup.F`, `fe_chem.F` (+3)
  - L106 Called from outside the package — `CALC_PCO2_APPROX`, `CALC_PCO2_SOLVESAPHE`, `CARBON_COEFFS`, `DIC_COEFFS_SURF`, `AHINI_FOR_AT`, `CARBON_COEFFS_PRESSURE_DEP`, `DIC_COEFFS_DEEP`, `ANW_INFSUP`, `EQUATION_AT`, `SOLVE_AT_FAST`, `SOLVE_AT_GENERAL`, `SOLVE_AT_GENERAL_SEC`, `DIC_FIELDS_LOAD`, `DIC_ATMOS`, `DIC_BIOTIC_FORCING`, `DIC_COST`, `DIC_STORE_FLUXCO2`, `DIC_INIT_FIXED`, `DIC_INIT_VARIA`, `DIC_SURFFORCING_INIT`, `DIC_READPARMS`, `DIC_TR_REGISTER`, `DIC_WRITE_PICKUP`
  - L139 Verification experiments compiling it (3) — `so_box_biogeo`, `tutorial_dic_adjoffline`, `tutorial_global_oce_biogeo`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/down_slope.md`
- L1 pkg/down_slope — `useDOWN_SLOPE`, `data.down_slope`
  - L10 Namelist parameters
    - L11 DWNSLP_PARM01 — `DWNSLP_slope`, `DWNSLP_rec_mu`, `DWNSLP_drFlow`, `temp_useDWNSLP`, `salt_useDWNSLP`
  - L18 Headers — `DWNSLP_OPTIONS.h`, `DWNSLP_PARAMS.h`, `DWNSLP_SIZE.h`, `DWNSLP_VARS.h`
  - L24 Routines (7) — `dwnslp_apply.F`, `dwnslp_calc_flow.F`, `dwnslp_calc_rho.F`, `dwnslp_diagnostics_init.F`, `dwnslp_init_fixed.F`, `dwnslp_init_varia.F`, `dwnslp_readparms.F`
  - L27 Called from outside the package — `DWNSLP_CALC_FLOW`, `DWNSLP_CALC_RHO`, `DWNSLP_INIT_FIXED`, `DWNSLP_INIT_VARIA`, `DWNSLP_READPARMS`, `DWNSLP_APPLY`
  - L37 Verification experiments compiling it (2) — `global_ocean.90x40x15`, `lab_sea`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/ebm.md`
- L1 pkg/ebm — `useEBM`
  - L10 Namelist parameters
    - L11 EBM_PARM01 — `tauThetaZonRelax`, `scale_runoff`, `RunoffFile`
  - L16 CPP options (defaults as shipped) — `EBM_WIND_PERT`, `EBM_CLIMATE_CHANGE`, `EBM_VERSION_1BASIN`
  - L21 Headers — `EBM.h`, `EBM_OPTIONS.h`, `ebm_ad_check_lev1_dir.h`, `ebm_ad_check_lev2_dir.h`, `ebm_ad_check_lev3_dir.h`
  - L28 Routines (8) — `ebm_area_t.F`, `ebm_atmosphere.F`, `ebm_driver.F`, `ebm_ini_vars.F`, `ebm_load_climatology.F`, `ebm_readparms.F`, `ebm_wind_perturb.F`, `ebm_zonalmean.F`
  - L31 Called from outside the package — `EBM_DRIVER`, `EBM_INI_VARS`, `EBM_READPARMS`
  - L36 Verification experiments compiling it (1) — `global_oce_latlon`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/ecco.md`
- L1 pkg/ecco — `useECCO`
  - L11 Namelist parameters
    - L12 ECCO_COST_NML — `data_errfile`, `costTranspDataFile`, `using_cost_altim`, `using_cost_bp`, `using_cost_scat`, `using_cost_sst`, `using_cost_seaice`, `using_cost_transp`, `using_mdt`, `using_tpj`, `using_ers`, `using_gfo`, `using_topex`, `tpTimeMaskFile`, `ersTimeMaskFile`, `gfoTimeMaskFile`, `scatx_errfile`, `scaty_errfile`, `ssh_errfile`, `sshv4cost_scalefile`, `sshv4cost_errfile`, `sshv4cost_errfactor`, `tp_errfile`, `ers_errfile`, `gfo_errfile`, `hflux_errfile`, `sflux_errfile`, `tauu_errfile`, `tauv_errfile`, `atemp_errfile` (+103)
    - L254 ECCO_GENCOST_NML — `using_gencost`, `gencost_barfile`, `gencost_datafile`, `gencost_name`, `gencost_scalefile`, `gencost_errfile`, `gencost_itracer`, `gencost_kLev_select`, `gencost_preproc`, `gencost_preproc_c`, `gencost_preproc_i`, `gencost_preproc_r`, `gencost_posproc`, `gencost_posproc_c`, `gencost_posproc_i`, `gencost_posproc_r`, `gencost_outputlevel`, `gencost_mask`, `gencost_spmin`, `gencost_spmax`, `gencost_spzero`, `gencost_wei1d`, `gencost_avgperiod`, `gencost_nrecperiod`, `gencost_startdate1`, `gencost_startdate2`, `gencost_enddate1`, `gencost_enddate2`, `gencost_smooth2Ddiffnbt`, `gencost_is1d` (+9)
  - L295 CPP options (defaults as shipped) — `ALLOW_GENCOST_CONTRIBUTION`, `ALLOW_GENCOST3D`, `ALLOW_GENCOST_1D`, `ALLOW_GENCOST_SSTV4_OUTPUT`, `ECCO_VARIABLE_AREAVOLGLOB`, `ALLOW_PSBAR_STERIC`, `ALLOW_IB_CORR`, `ALLOW_SHALLOW_ALTIMETRY`, `ALLOW_HIGHLAT_ALTIMETRY`, `ALLOW_ECCO_OLD_FC_PRINT`, `ECCO_VERBOSE`, `ALLOW_ECCO_DEBUG`, `ALLOW_SSH_COST_CONTRIBUTION`, `ALLOW_SST_COST_CONTRIBUTION`, `ALLOW_SEAICE_COST_CONTRIBUTION`
  - L312 Headers — `ECCO.h`, `ECCO_OPTIONS.h`, `ECCO_SIZE.h`, `ecco_cost.h`, `ecco_local_params.h`
  - L319 Routines (69) — `cost_averagesfields.F`, `cost_averagesflags.F`, `cost_averagesgeneric.F`, `cost_averagesinit.F`, `cost_bp_read.F`, `cost_forcing.F`, `cost_gencal.F`, `cost_gencost_all.F`, `cost_gencost_assignperiod.F`, `cost_gencost_boxmean.F`, `cost_gencost_bpv4.F`, `cost_gencost_customize.F`, `cost_gencost_glbmean.F`, `cost_gencost_moc.F`, `cost_gencost_seaicev4.F`, `cost_gencost_sshv4.F`, `cost_gencost_sstv4.F`, `cost_gencost_transp.F`, `cost_generic.F`, `cost_genread.F`, `cost_hyd.F`, `cost_internal_params.F`, `cost_obcs.F`, `cost_obcs_ageos.F`, `cost_obcse.F`, `cost_obcsn.F`, `cost_obcss.F`, `cost_obcsvol.F`, `cost_obcsw.F`, `cost_sla_read.F` (+20)
  - L322 Called from outside the package — `STERGLOH_OUTPUT`, `ECCO_PHYS`, `ECCO_CHECK`, `ECCO_INIT_FIXED`, `ECCO_INIT_VARIA`, `ECCO_READPARMS`, `ECCO_WRITE_PICKUP`, `COST_AVERAGESFIELDS`, `COST_GENCOST_ALL`, `ECCO_COST_FINAL`
  - L334 Verification experiments compiling it (4) — `global_oce_biogeo_bling`, `lab_sea`, `obcs_ctrl`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/embed_files.md`
- L1 pkg/embed_files — `useEMBED_FILES`
  - L7 README
  - L20 Headers — `EMBED_FILES_OPTIONS.h`
  - L23 Routines (1) — `embed_files_init.F`
  - L26 Called from outside the package — `EMBED_FILES_INIT`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/exch2.md`
- L1 pkg/exch2 — `useEXCH2`
  - L10 Namelist parameters
    - L11 W2_EXCH2_PARM01 — `preDefTopol`, `dimsFacets`, `facetEdgeLink`, `blankList`, `W2_mapIO`, `W2_printMsg`, `W2_useE2ioLayOut`
  - L20 CPP options (defaults as shipped) — `W2_USE_E2_SAFEMODE`, `W2_E2_DEBUG_ON`, `W2_USE_R1_ONLY`, `W2_FILL_NULL_REGIONS`, `W2_CUMSUM_USE_MATRIX`
  - L27 Headers — `W2_EXCH2_BUFFER.h`, `W2_EXCH2_PARAMS.h`, `W2_EXCH2_SIZE.h`, `W2_EXCH2_TOPOLOGY.h`, `W2_OPTIONS.h`
  - L34 Routines (127) — `exch2_3d_r4.F`, `exch2_3d_r8.F`, `exch2_3d_rl.F`, `exch2_3d_rs.F`, `exch2_ad_get_r41.F`, `exch2_ad_get_r42.F`, `exch2_ad_get_r81.F`, `exch2_ad_get_r82.F`, `exch2_ad_get_rl1.F`, `exch2_ad_get_rl2.F`, `exch2_ad_get_rs1.F`, `exch2_ad_get_rs2.F`, `exch2_ad_put_r41.F`, `exch2_ad_put_r42.F`, `exch2_ad_put_r81.F`, `exch2_ad_put_r82.F`, `exch2_ad_put_rl1.F`, `exch2_ad_put_rl2.F`, `exch2_ad_put_rs1.F`, `exch2_ad_put_rs2.F`, `exch2_check_depths.F`, `exch2_get_r41.F`, `exch2_get_r42.F`, `exch2_get_r81.F`, `exch2_get_r82.F`, `exch2_get_rl1.F`, `exch2_get_rl2.F`, `exch2_get_rs1.F`, `exch2_get_rs2.F`, `exch2_get_scal_bounds.F` (+93)
  - L37 Called from outside the package — `W2_CUMULSUM_Z_TILE_RL`, `W2_EEBOOT`, `EXCH2_CHECK_DEPTHS`
  - L42 Verification experiments compiling it (15) — `advect_cs`, `aim.5l_cs`, `atm_gray`, `global_oce_biogeo_bling`, `global_ocean.90x40x15`, `global_ocean.cs32x15`, `tutorial_held_suarez_cs`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/exch2@bbl.md`
- L1 pkg/exch2  (bbl: ~/Documents/research/ECCO/BBL/MITgcm) — `useEXCH2`
  - L12 Namelist parameters
    - L13 W2_EXCH2_PARM01 — `preDefTopol`, `dimsFacets`, `facetEdgeLink`, `blankList`, `W2_mapIO`, `W2_printMsg`, `W2_useE2ioLayOut`
  - L22 CPP options (defaults as shipped) — `W2_USE_E2_SAFEMODE`, `W2_E2_DEBUG_ON`, `W2_USE_R1_ONLY`, `W2_FILL_NULL_REGIONS`, `W2_CUMSUM_USE_MATRIX`
  - L29 Headers — `W2_EXCH2_BUFFER.h`, `W2_EXCH2_PARAMS.h`, `W2_EXCH2_SIZE.h`, `W2_EXCH2_TOPOLOGY.h`, `W2_OPTIONS.h`
  - L36 Routines (127) — `exch2_3d_r4.F`, `exch2_3d_r8.F`, `exch2_3d_rl.F`, `exch2_3d_rs.F`, `exch2_ad_get_r41.F`, `exch2_ad_get_r42.F`, `exch2_ad_get_r81.F`, `exch2_ad_get_r82.F`, `exch2_ad_get_rl1.F`, `exch2_ad_get_rl2.F`, `exch2_ad_get_rs1.F`, `exch2_ad_get_rs2.F`, `exch2_ad_put_r41.F`, `exch2_ad_put_r42.F`, `exch2_ad_put_r81.F`, `exch2_ad_put_r82.F`, `exch2_ad_put_rl1.F`, `exch2_ad_put_rl2.F`, `exch2_ad_put_rs1.F`, `exch2_ad_put_rs2.F`, `exch2_check_depths.F`, `exch2_get_r41.F`, `exch2_get_r42.F`, `exch2_get_r81.F`, `exch2_get_r82.F`, `exch2_get_rl1.F`, `exch2_get_rl2.F`, `exch2_get_rs1.F`, `exch2_get_rs2.F`, `exch2_get_scal_bounds.F` (+93)
  - L39 Called from outside the package — `W2_CUMULSUM_Z_TILE_RL`, `W2_EEBOOT`, `EXCH2_3D_R4`, `EXCH2_3D_R8`, `EXCH2_3D_RL`, `EXCH2_3D_RS`, `EXCH2_S3D_R4`, `EXCH2_S3D_R8`, `EXCH2_S3D_RL`, `EXCH2_S3D_RS`, `EXCH2_SM_3D_R4`, `EXCH2_SM_3D_R8`, `EXCH2_SM_3D_RL`, `EXCH2_SM_3D_RS`, `EXCH2_UV_3D_R4`, `EXCH2_UV_CGRID_3D_R4`, `EXCH2_UV_3D_R8`, `EXCH2_UV_CGRID_3D_R8`, `EXCH2_UV_3D_RL`, `EXCH2_UV_CGRID_3D_RL`, `EXCH2_UV_3D_RS`, `EXCH2_UV_CGRID_3D_RS`, `EXCH2_UV_AGRID_3D_R4`, `EXCH2_UV_AGRID_3D_R8`, `EXCH2_UV_AGRID_3D_RL`, `EXCH2_UV_AGRID_3D_RS`, `EXCH2_UV_BGRID_3D_R4`, `EXCH2_UV_BGRID_3D_R8`, `EXCH2_UV_BGRID_3D_RL`, `EXCH2_UV_BGRID_3D_RS` (+9)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/exf.md`
- L1 pkg/exf — `useEXF`
  - L11 Namelist parameters
    - L12 EXF_NML_01 — `repeatPeriod`, `exf_albedo`, `ocean_emissivity`, `ice_emissivity`, `snow_emissivity`, `exf_iceCd`, `exf_iceCe`, `exf_iceCh`, `exf_scal_BulkCdn`, `useExfCheckRange`, `exf_iprec`, `exf_iprec_obcs`, `exf_yftype`, `exf_verbose`, `exf_debugLev`, `exf_monFreq`, `exf_adjMonFreq`, `exf_adjMonSelect`, `useExfYearlyFields`, `twoDigitYear`, `readStressOnAgrid`, `rotateStressOnAgrid`, `readStressOnCgrid`, `useAtmWind`, `useRelativeWind`, `noNegativeEvap`, `useStabilityFct_overIce`, `diags_opOceWeighted`, `select_ZenAlbedo`, `useExfZenIncoming` (+16)
    - L74 EXF_NML_02 — `hfluxRepCycle`, `hfluxStartTime`, `atempRepCycle`, `atempStartTime`, `aqhRepCycle`, `aqhStartTime`, `hs_file`, `hs_startdate1`, `hs_startdate2`, `hs_RepCycle`, `hs_period`, `hs_StartTime`, `hl_file`, `hl_startdate1`, `hl_startdate2`, `hl_RepCycle`, `hl_period`, `hl_StartTime`, `sfluxRepCycle`, `sfluxStartTime`, `evapRepCycle`, `evapStartTime`, `precipRepCycle`, `precipStartTime`, `snowprecipRepCycle`, `snowprecipStartTime`, `runoffRepCycle`, `runoffStartTime`, `saltflxRepCycle`, `saltflxStartTime` (+41)
    - L243 EXF_NML_03 — `exf_inscal_hflux`, `exf_inscal_sflux`, `exf_inscal_evap`, `exf_inscal_ustress`, `exf_inscal_vstress`, `exf_inscal_uwind`, `exf_inscal_vwind`, `exf_inscal_wspeed`, `exf_inscal_atemp`, `exf_offset_atemp`, `exf_inscal_aqh`, `exf_inscal_hs`, `exf_inscal_hl`, `exf_inscal_sst`, `exf_inscal_sss`, `exf_inscal_swflux`, `exf_inscal_lwflux`, `exf_inscal_precip`, `exf_inscal_runoff`, `exf_inscal_apressure`, `exf_inscal_snowprecip`, `exf_inscal_runoftemp`, `exf_inscal_saltflx`, `exf_inscal_swdown`, `exf_inscal_lwdown`, `exf_inscal_climsst`, `exf_inscal_climsss`, `exf_inscal_climustr`, `exf_inscal_climvstr`, `exf_outscal_hflux` (+70)
    - L369 EXF_NML_04 — `ustress_lon0`, `ustress_lon_inc`, `ustress_lat0`, `ustress_lat_inc`, `ustress_nlon`, `ustress_nlat`, `ustress_interpMethod`, `vstress_lon0`, `vstress_lon_inc`, `vstress_lat0`, `vstress_lat_inc`, `vstress_nlon`, `vstress_nlat`, `vstress_interpMethod`, `hflux_lon0`, `hflux_lon_inc`, `hflux_lat0`, `hflux_lat_inc`, `hflux_nlon`, `hflux_nlat`, `hflux_interpMethod`, `sflux_lon0`, `sflux_lon_inc`, `sflux_lat0`, `sflux_lat_inc`, `sflux_nlon`, `sflux_nlat`, `sflux_interpMethod`, `swflux_lon0`, `swflux_lon_inc` (+160)
    - L560 EXF_NML_OBCS — `useOBCSYearlyFields`, `obcsNstartdate1`, `obcsNstartdate2`, `obcsNstartTime`, `obcsNperiod`, `obcsNrepCycle`, `obcsSstartdate1`, `obcsSstartdate2`, `obcsSstartTime`, `obcsSperiod`, `obcsSrepCycle`, `obcsEstartdate1`, `obcsEstartdate2`, `obcsEstartTime`, `obcsEperiod`, `obcsErepCycle`, `obcsWstartdate1`, `obcsWstartdate2`, `obcsWstartTime`, `obcsWperiod`, `obcsWrepCycle`, `siobNstartdate1`, `siobNstartdate2`, `siobNstartTime`, `siobNperiod`, `siobNrepCycle`, `siobSstartdate1`, `siobSstartdate2`, `siobSstartTime`, `siobSperiod` (+11)
  - L603 CPP options (defaults as shipped) — `ALLOW_ATM_TEMP`, `ALLOW_ATM_WIND`, `ALLOW_DOWNWARD_RADIATION`, `ALLOW_BULKFORMULAE`, `ALLOW_BULK_LARGEYEAGER04`, `ALLOW_DRAG_LARGEYEAGER09`, `EXF_READ_EVAP`, `ALLOW_READ_TURBFLUXES`, `ALLOW_RUNOFF`, `ALLOW_RUNOFTEMP`, `ALLOW_SALTFLX`, `EXF_CALC_ATMRHO`, `ALLOW_ZENITHANGLE`, `EXF_LWDOWN_WITH_EMISSIVITY`, `ALLOW_CLIMSST_RELAXATION`, `ALLOW_CLIMSSS_RELAXATION`, `EXF_ALLOW_TIDES`, `EXF_SEAICE_FRACTION`, `USE_EXF_INTERPOLATION`, `EXF_USE_OLD_VEC_ROTATION`, `EXF_USE_OLD_INTERP_POLE`, `EXF_INTERP_USE_DYNALLOC`, `EXF_IREAD_USE_GLOBAL_POINTER`, `ALLOW_BULK_OFFLINE`, `ALLOW_CLIMSTRESS_RELAXATION`
  - L630 Headers — `EXF_CONSTANTS.h`, `EXF_FIELDS.h`, `EXF_INTERP_PARAM.h`, `EXF_INTERP_SIZE.h`, `EXF_OPTIONS.h`, `EXF_PARAM.h`, `exf_ad_check_lev1_dir.h`, `exf_ad_check_lev2_dir.h`, `exf_ad_check_lev3_dir.h`, `exf_ad_check_lev4_dir.h`
  - L642 Routines (52) — `exf_adjoint_snapshots.F`, `exf_adjoint_snapshots__g.F`, `exf_adjoint_snapshots_ad.F`, `exf_bulkformulae.F`, `exf_check.F`, `exf_check_range.F`, `exf_diagnostics_fill.F`, `exf_diagnostics_init.F`, `exf_filter_rl.F`, `exf_getclim.F`, `exf_getffield_start.F`, `exf_getffieldrec.F`, `exf_getffields.F`, `exf_getfield_start.F`, `exf_getforcing.F`, `exf_getmonthsrec.F`, `exf_getsurfacefluxes.F`, `exf_getyearlyfieldname.F`, `exf_init_fixed.F`, `exf_init_fld.F`, `exf_init_interp.F`, `exf_init_varia.F`, `exf_interp.F`, `exf_interp_read.F`, `exf_interp_uv.F`, `exf_interpolate.F`, `exf_mapfields.F`, `exf_monitor.F`, `exf_monitor_ad.F`, `exf_radiation.F` (+11)
  - L645 Called from outside the package — `EXF_GETFORCING`, `EXF_CHECK`, `EXF_INIT_FIXED`, `EXF_INIT_VARIA`, `EXF_READPARMS`, `EXF_SET_FLD`, `EXF_INIT_FLD`, `EXF_GETFFIELD_START`, `EXF_GETFFIELDREC`, `EXF_GETMONTHSREC`, `EXF_SET_OBCS_XZ`, `EXF_SET_OBCS_YZ`, `EXF_ADJOINT_SNAPSHOTS`, `EXF_WEIGHT_SFX_DIAGS`
  - L662 Verification experiments compiling it (8) — `global_oce_latlon`, `global_ocean.cs32x15`, `lab_sea`, `obcs_ctrl`, `offline_exf_seaice`, `seaice_itd`, `seaice_obcs`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/fizhi.md`
- L1 pkg/fizhi — `useFIZHI`
  - L10 Namelist parameters
    - L11 FIZHI_LIST — `fizhi_mnc_write_pickup`, `fizhi_mnc_read_pickup`
  - L20 CPP options (defaults as shipped) — `FIZHI_USE_FIXED_DAY`, `TRY_NEW_GETPWHERE`, `FIZHI_F77_COMPIL`, `FIZHI_CRAY`, `FIZHI_SGI`, `FIZHI_TRBFLX_OLD_BUG`
  - L28 Headers — `FIZHI_OPTIONS.h`, `chronos.h`, `fizhi_SHP.h`, `fizhi_SIZE.h`, `fizhi_chemistry_coms.h`, `fizhi_coms.h`, `fizhi_earth_coms.h`, `fizhi_io_comms.h`, `fizhi_land_SIZE.h`, `fizhi_land_coms.h`, `fizhi_ocean_coms.h`, `sibber.h`, `snwmid.h`
  - L48 Routines (156) — `AtoC.F`, `CtoA.F`, `do_fizhi.F`, `fizhi_alarms.F`, `fizhi_clockstuff.F`, `fizhi_diagalarms.F`, `fizhi_diagnostics_init.F`, `fizhi_driver.F`, `fizhi_fillnegs.F`, `fizhi_gwdrag.F`, `fizhi_init_chem.F`, `fizhi_init_fixed.F`, `fizhi_init_vars.F`, `fizhi_init_veg.F`, `fizhi_init_vegsurftiles.F`, `fizhi_lsm.F`, `fizhi_lwrad.F`, `fizhi_mnc_init.F`, `fizhi_moist.F`, `fizhi_mpistuff.F`, `fizhi_rayleigh.F`, `fizhi_read_pickup.F`, `fizhi_readparms.F`, `fizhi_readwrite_vegtiles.F`, `fizhi_step_diag.F`, `fizhi_swrad.F`, `fizhi_tendency_apply.F`, `fizhi_turb.F`, `fizhi_update_time.F`, `fizhi_utils.F` (+13)
  - L51 Called from outside the package — `FIZHI_TENDENCY_APPLY_S`, `FIZHI_TENDENCY_APPLY_T`, `FIZHI_TENDENCY_APPLY_U`, `FIZHI_TENDENCY_APPLY_V`, `FIZHI_UPDATE_TIME`, `FIZHI_WRAPPER`, `STEP_FIZHI_FG`, `UPDATE_CHEMISTRY_EXPORTS`, `UPDATE_EARTH_EXPORTS`, `UPDATE_OCEAN_EXPORTS`, `FIZHI_WRITE_STATE`, `STEP_FIZHI_CORR`, `FIZHI_INIT_FIXED`, `FIZHI_INIT_VARS`, `FIZHI_READPARMS`, `FIZHI_WRITE_DATETIME`, `FIZHI_WRITE_PICKUP`, `FIZHI_WRITE_VEGTILES`, `MY_EXIT`, `MY_FINALIZE`, `FIZHI_DIAGALARMS`
  - L79 Verification experiments compiling it (3)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/flt.md`
- L1 pkg/flt — `useFLT`
  - L10 README
  - L44 Namelist parameters
    - L45 FLT_NML — `flt_int_traj`, `flt_int_prof`, `flt_selectTrajOutp`, `flt_selectProfOutp`, `flt_noise`, `flt_deltaT`, `FLT_Iter0`, `flt_file`, `mapIniPos2Index`
  - L56 CPP options (defaults as shipped) — `ALLOW_3D_FLT`, `USE_FLT_ALT_NOISE`, `ALLOW_FLT_3D_NOISE`, `FLT_SECOND_ORDER_RUNGE_KUTTA`, `FLT_WITHOUT_X_PERIODICITY`, `FLT_WITHOUT_Y_PERIODICITY`, `DEVEL_FLT_EXCH2`
  - L65 Headers — `FLT.h`, `FLT_BUFF.h`, `FLT_OPTIONS.h`, `FLT_SIZE.h`
  - L71 Routines (25) — `exch2_recv_get_vec.F`, `exch2_send_put_vec.F`, `exch_recv_get_vec.F`, `exch_send_put_vec.F`, `flt_down.F`, `flt_exch2.F`, `flt_exchg.F`, `flt_init_fixed.F`, `flt_init_varia.F`, `flt_interp_linear.F`, `flt_main.F`, `flt_mapping.F`, `flt_readparms.F`, `flt_runga2.F`, `flt_runga4.F`, `flt_traj.F`, `flt_up.F`, `flt_write_pickup.F`
  - L74 Called from outside the package — `FLT_MAIN`, `FLT_INIT_FIXED`, `FLT_INIT_VARIA`, `FLT_READPARMS`, `FLT_WRITE_PICKUP`
  - L81 Verification experiments compiling it (2)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/frazil.md`
- L1 pkg/frazil — `useFRAZIL`
  - L8 Headers — `FRAZIL.h`, `FRAZIL_OPTIONS.h`
  - L12 Routines (5) — `frazil_calc_rhs.F`, `frazil_diagnostics_init.F`, `frazil_init_fixed.F`, `frazil_init_varia.F`, `frazil_tendency_apply.F`
  - L15 Called from outside the package — `FRAZIL_TENDENCY_APPLY_T`, `FRAZIL_CALC_RHS`, `FRAZIL_INIT_FIXED`, `FRAZIL_INIT_VARIA`
  - L22 Verification experiments compiling it (1) — `global_oce_latlon`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/gchem.md`
- L1 pkg/gchem — `useGCHEM`
  - L11 Namelist parameters
    - L12 GCHEM_PARM01 — `useCFC`, `useDIC`, `useBLING`, `useSPOIL`, `useDARWIN`, `fileName1`, `fileName2`, `fileName3`, `fileName4`, `fileName5`, `gchem_int1`, `gchem_int2`, `gchem_int3`, `gchem_int4`, `gchem_int5`, `gchem_rl1`, `gchem_rl2`, `gchem_rl3`, `gchem_rl4`, `gchem_rl5`, `gchem_ForcingPeriod`, `gchem_ForcingCycle`, `gchem_secondsPerYear`, `tIter0`
  - L39 CPP options (defaults as shipped) — `GCHEM_SEPARATE_FORCING`, `GCHEM_ADD2TR_TENDENCY`
  - L45 Headers — `GCHEM.h`, `GCHEM_FIELDS.h`, `GCHEM_OPTIONS.h`, `GCHEM_SIZE.h`
  - L51 Routines (15) — `gchem_add_tendency.F`, `gchem_calc_tendency.F`, `gchem_check.F`, `gchem_cons.F`, `gchem_diagnostics_init.F`, `gchem_fields_load.F`, `gchem_forcing_sep.F`, `gchem_init_fixed.F`, `gchem_init_vari.F`, `gchem_insolation.F`, `gchem_output.F`, `gchem_readparms.F`, `gchem_surfmean.F`, `gchem_tr_register.F`, `gchem_write_pickup.F`
  - L54 Called from outside the package — `GCHEM_OUTPUT`, `GCHEM_CALC_TENDENCY`, `GCHEM_FORCING_SEP`, `GCHEM_FIELDS_LOAD`, `GCHEM_CHECK`, `GCHEM_INIT_FIXED`, `GCHEM_INIT_VARI`, `GCHEM_READPARMS`, `GCHEM_WRITE_PICKUP`, `GCHEM_INSOLATION`, `GCHEM_ADD_TENDENCY`
  - L68 Verification experiments compiling it (6) — `cfc_example`, `global_oce_biogeo_bling`, `so_box_biogeo`, `tutorial_cfc_offline`, `tutorial_dic_adjoffline`, `tutorial_global_oce_biogeo`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/generic_advdiff.md`
- L1 pkg/generic_advdiff
  - L10 CPP options (defaults as shipped) — `COSINEMETH_III`, `ISOTROPIC_COS_SCALING`, `DISABLE_MULTIDIM_ADVECTION`, `GAD_MULTIDIM_COMPRESSIBLE`, `GAD_ALLOW_TS_SOM_ADV`, `GAD_SMOLARKIEWICZ_HACK`
  - L19 Headers — `GAD.h`, `GAD_FLUX_LIMITER.h`, `GAD_OPTIONS.h`, `GAD_SOM_VARS.h`, `gad_ad_check_lev1_dir.h`, `gad_ad_check_lev2_dir.h`, `gad_ad_check_lev3_dir.h`, `gad_ad_check_lev4_dir.h`
  - L29 Routines (101) — `gad_advection.F`, `gad_advscheme.F`, `gad_biharm_r.F`, `gad_biharm_x.F`, `gad_biharm_y.F`, `gad_c2_adv_r.F`, `gad_c2_adv_x.F`, `gad_c2_adv_y.F`, `gad_c2_impl_r.F`, `gad_c4_adv_r.F`, `gad_c4_adv_x.F`, `gad_c4_adv_y.F`, `gad_calc_rhs.F`, `gad_check.F`, `gad_del2.F`, `gad_diagnostics_init.F`, `gad_diagnostics_state.F`, `gad_diff_r.F`, `gad_diff_x.F`, `gad_diff_y.F`, `gad_dst2u1_adv_r.F`, `gad_dst2u1_adv_x.F`, `gad_dst2u1_adv_y.F`, `gad_dst2u1_impl_r.F`, `gad_dst3_adv_r.F`, `gad_dst3_adv_x.F`, `gad_dst3_adv_y.F`, `gad_dst3fl_adv_r.F`, `gad_dst3fl_adv_x.F`, `gad_dst3fl_adv_y.F` (+62)
  - L32 Called from outside the package — `GAD_SOM_EXCHANGES`, `GAD_DIAGNOSTICS_STATE`, `GAD_CHECK`, `GAD_INIT_FIXED`, `GAD_INIT_VARIA`, `GAD_WRITE_PICKUP`, `GAD_ADVECTION`, `GAD_CALC_RHS`, `GAD_IMPLICIT_R`, `GAD_SOM_ADVECT`, `SALT_FILL`, `GAD_C2_ADV_X`, `GAD_C2_ADV_Y`, `GAD_DIFF_X`, `GAD_DIFF_Y`, `GAD_DST3FL_ADV_X`, `GAD_DST3FL_ADV_Y`, `GAD_EXCH_SOM`, `GAD_DST2U1_ADV_X`, `GAD_DST2U1_ADV_Y`, `GAD_DST3_ADV_X`, `GAD_DST3_ADV_Y`, `GAD_FLUXLIMIT_ADV_X`, `GAD_FLUXLIMIT_ADV_Y`, `GAD_OS7MP_ADV_X`, `GAD_OS7MP_ADV_Y`, `GAD_PPM_ADV_X`, `GAD_PPM_ADV_Y`, `GAD_PQM_ADV_X`, `GAD_PQM_ADV_Y`
  - L85 Verification experiments compiling it (56) — `advect_cs`, `advect_xz`, `aim.5l_Equatorial_Channel`, `aim.5l_LatLon`, `aim.5l_cs`, `atm_gray`, `bottom_ctrl_5x5`, `cfc_example`, `cheapAML_box`, `deep_anelastic`, `front_relax`, `global_oce_biogeo_bling`, `global_oce_latlon`, `global_ocean.90x40x15`, `global_ocean.cs32x15`, `ideal_2D_oce`, `internal_wave`, `inverted_barometer`, `lab_sea`, `matrix_example`, `obcs_ctrl`, `offline_exf_seaice`, `seaice_itd`, `seaice_obcs`, `shelfice_2d_remesh`, `short_surf_wave`, `so_box_biogeo`, `tutorial_advection_in_gyre`, `tutorial_baroclinic_gyre`, `tutorial_cfc_offline` (+11)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/ggl90.md`
- L1 pkg/ggl90 — `useGGL90`
  - L10 Namelist parameters
    - L11 GGL90_PARM01 — `GGL90dumpFreq`, `GGL90taveFreq`, `GGL90diffTKEh`, `GGL90mixingMaps`, `GGL90writeState`, `GGL90TKEsurfMin`, `GGL90mixingLengthMin`, `mxlMaxFlag`, `adMxlMaxFlag`, `mxlSurfFlag`, `GGL90viscMax`, `GGL90diffMax`, `GGL90_dirichlet`, `calcMeanVertShear`, `useIDEMIX`, `useLANGMUIR`
    - L35 GGL90_PARM02 — `IDEMIX_tau_v`, `IDEMIX_tau_h`, `IDEMIX_gamma`, `IDEMIX_jstar`, `IDEMIX_mu0`, `IDEMIX_tidal_file`, `IDEMIX_wind_file`, `IDEMIX_mixing_efficiency`, `IDEMIX_diff_max`, `IDEMIX_diff_min`, `IDEMIX_frac_F_b`, `IDEMIX_frac_F_s`, `IDEMIX_include_GM`, `IDEMIX_include_GM_bottom`
    - L50 GGL90_PARM03 — `LC_Gamma`, `LC_num`, `LC_lambda`
  - L55 CPP options (defaults as shipped) — `ALLOW_GGL90_HORIZDIFF`, `ALLOW_GGL90_SMOOTH`, `ALLOW_GGL90_IDEMIX`, `GGL90_IDEMIX_CVMIX_VERSION`, `ALLOW_GGL90_LANGMUIR`, `GGL90_REGULARIZE_MIXINGLENGTH`, `GGL90_MISSING_HFAC_BUG`
  - L64 Headers — `GGL90.h`, `GGL90_OPTIONS.h`, `ggl90_ad_check_lev1_dir.h`, `ggl90_ad_check_lev2_dir.h`, `ggl90_ad_check_lev3_dir.h`, `ggl90_ad_check_lev4_dir.h`
  - L72 Routines (17) — `ggl90_add_stokesdrift.F`, `ggl90_calc.F`, `ggl90_calc_diff.F`, `ggl90_calc_visc.F`, `ggl90_check.F`, `ggl90_diagnostics_init.F`, `ggl90_exchanges.F`, `ggl90_idemix.F`, `ggl90_init_fixed.F`, `ggl90_init_varia.F`, `ggl90_mixinglength.F`, `ggl90_output.F`, `ggl90_read_pickup.F`, `ggl90_readparms.F`, `ggl90_write_pickup.F`
  - L75 Called from outside the package — `GGL90_CALC_DIFF`, `GGL90_CALC_VISC`, `GGL90_CALC`, `GGL90_EXCHANGES`, `GGL90_OUTPUT`, `GGL90_CHECK`, `GGL90_INIT_FIXED`, `GGL90_INIT_VARIA`, `GGL90_READPARMS`, `GGL90_WRITE_PICKUP`, `GGL90_ADD_STOKESDRIFT`
  - L89 Verification experiments compiling it (5) — `global_oce_latlon`, `global_ocean.90x40x15`, `global_ocean.cs32x15`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/gmredi.md`
- L1 pkg/gmredi — `useGMREDI`
  - L11 README
  - L31 Namelist parameters
    - L32 GM_PARM01 — `GM_AdvForm`, `GM_AdvSeparate`, `GM_InMomAsStress`, `GM_isopycK`, `GM_background_K`, `GM_iso2dFile`, `GM_iso1dFile`, `GM_bol2dFile`, `GM_bol1dFile`, `GM_K3dRediFile`, `GM_K3dGMFile`, `GM_background_K3dFile`, `GM_isopycK3dFile`, `GM_taper_scheme`, `GM_maxSlope`, `GM_Kmin_horiz`, `GM_Small_Number`, `GM_slopeSqCutoff`, `GM_Scrit`, `GM_Sd`, `GM_isoFac_calcK`, `GM_facTrL2dz`, `GM_facTrL2ML`, `GM_maxTransLay`, `GM_UseBVP`, `GM_BVP_cMin`, `GM_BVP_ModeNumber`, `GM_useSubMeso`, `subMeso_Ceff`, `subMeso_invTau` (+45)
  - L109 CPP options (defaults as shipped) — `GM_READ_K3D_REDI`, `GM_READ_K3D_GM`, `GM_VISBECK_VARIABLE_K`, `GM_GEOM_VARIABLE_K`, `GM_BATES_K3D`, `GM_BATES_PASSIVE`, `GM_NON_UNITY_DIAGONAL`, `GM_EXTRA_DIAGONAL`, `GM_BOLUS_ADVEC`, `GM_BOLUS_BVP`, `ALLOW_GM_LEITH_QG`, `GM_AUTODIFF_EXCESSIVE_STORE`, `GMREDI_MASK_SLOPES`
  - L124 Headers — `GMREDI.h`, `GMREDI_OPTIONS.h`
  - L128 Routines (30) — `gmredi_calc_bates_k.F`, `gmredi_calc_diff.F`, `gmredi_calc_eigs.F`, `gmredi_calc_geom.F`, `gmredi_calc_psi_b.F`, `gmredi_calc_psi_bolus.F`, `gmredi_calc_psi_bvp.F`, `gmredi_calc_qgleith.F`, `gmredi_calc_tensor.F`, `gmredi_calc_urms.F`, `gmredi_check.F`, `gmredi_diagnostics_fill.F`, `gmredi_diagnostics_impl.F`, `gmredi_diagnostics_init.F`, `gmredi_do_exch.F`, `gmredi_init_fixed.F`, `gmredi_init_varia.F`, `gmredi_mnc_init.F`, `gmredi_output.F`, `gmredi_read_pickup.F`, `gmredi_readparms.F`, `gmredi_residual_flow.F`, `gmredi_rtransport.F`, `gmredi_slope_limit.F`, `gmredi_slope_psi.F`, `gmredi_write_pickup.F`, `gmredi_xtransport.F`, `gmredi_ytransport.F`, `submeso_calc_psi.F`
  - L131 Called from outside the package — `GMREDI_CALC_DIFF`, `GMREDI_CALC_TENSOR`, `GMREDI_CALC_TENSOR_DUMMY`, `GMREDI_DO_EXCH`, `GMREDI_DIAGNOSTICS_IMPL`, `GMREDI_OUTPUT`, `GMREDI_CHECK`, `GMREDI_INIT_FIXED`, `GMREDI_INIT_VARIA`, `GMREDI_READPARMS`, `GMREDI_WRITE_PICKUP`, `GMREDI_RESIDUAL_FLOW`, `GMREDI_RTRANSPORT`, `GMREDI_XTRANSPORT`, `GMREDI_YTRANSPORT`
  - L148 Verification experiments compiling it (20) — `cfc_example`, `front_relax`, `global_oce_biogeo_bling`, `global_oce_latlon`, `global_ocean.90x40x15`, `global_ocean.cs32x15`, `ideal_2D_oce`, `lab_sea`, `seaice_obcs`, `so_box_biogeo`, `tutorial_cfc_offline`, `tutorial_dic_adjoffline`, `tutorial_global_oce_biogeo`, `tutorial_global_oce_latlon`, `tutorial_global_oce_optim`, `tutorial_reentrant_channel`, `tutorial_tracer_adjsens`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/gmredi@wad.md`
- L1 pkg/gmredi  (wad: ~/Documents/research/ECCO/wetting_drying/MITgcm) — `useGMREDI`
  - L13 README
  - L33 Namelist parameters
    - L34 GM_PARM01 — `GM_AdvForm`, `GM_AdvSeparate`, `GM_InMomAsStress`, `GM_isopycK`, `GM_background_K`, `GM_iso2dFile`, `GM_iso1dFile`, `GM_bol2dFile`, `GM_bol1dFile`, `GM_K3dRediFile`, `GM_K3dGMFile`, `GM_background_K3dFile`, `GM_isopycK3dFile`, `GM_taper_scheme`, `GM_maxSlope`, `GM_Kmin_horiz`, `GM_Small_Number`, `GM_slopeSqCutoff`, `GM_Scrit`, `GM_Sd`, `GM_isoFac_calcK`, `GM_facTrL2dz`, `GM_facTrL2ML`, `GM_maxTransLay`, `GM_UseBVP`, `GM_BVP_cMin`, `GM_BVP_ModeNumber`, `GM_useSubMeso`, `subMeso_Ceff`, `subMeso_invTau` (+45)
  - L111 CPP options (defaults as shipped) — `GM_READ_K3D_REDI`, `GM_READ_K3D_GM`, `GM_VISBECK_VARIABLE_K`, `GM_GEOM_VARIABLE_K`, `GM_BATES_K3D`, `GM_BATES_PASSIVE`, `GM_NON_UNITY_DIAGONAL`, `GM_EXTRA_DIAGONAL`, `GM_BOLUS_ADVEC`, `GM_BOLUS_BVP`, `ALLOW_GM_LEITH_QG`, `GM_AUTODIFF_EXCESSIVE_STORE`, `GMREDI_MASK_SLOPES`
  - L126 Headers — `GMREDI.h`, `GMREDI_OPTIONS.h`
  - L130 Routines (29) — `gmredi_calc_bates_k.F`, `gmredi_calc_diff.F`, `gmredi_calc_eigs.F`, `gmredi_calc_geom.F`, `gmredi_calc_psi_bolus.F`, `gmredi_calc_psi_bvp.F`, `gmredi_calc_qgleith.F`, `gmredi_calc_tensor.F`, `gmredi_calc_urms.F`, `gmredi_check.F`, `gmredi_diagnostics_fill.F`, `gmredi_diagnostics_impl.F`, `gmredi_diagnostics_init.F`, `gmredi_do_exch.F`, `gmredi_init_fixed.F`, `gmredi_init_varia.F`, `gmredi_mnc_init.F`, `gmredi_output.F`, `gmredi_read_pickup.F`, `gmredi_readparms.F`, `gmredi_residual_flow.F`, `gmredi_rtransport.F`, `gmredi_slope_limit.F`, `gmredi_slope_psi.F`, `gmredi_write_pickup.F`, `gmredi_xtransport.F`, `gmredi_ytransport.F`, `submeso_calc_psi.F`
  - L133 Called from outside the package — `GMREDI_CALC_DIFF`, `GMREDI_CALC_TENSOR`, `GMREDI_CALC_TENSOR_DUMMY`, `GMREDI_DO_EXCH`, `GMREDI_DIAGNOSTICS_IMPL`, `GMREDI_OUTPUT`, `GMREDI_CHECK`, `GMREDI_INIT_FIXED`, `GMREDI_INIT_VARIA`, `GMREDI_READPARMS`, `GMREDI_WRITE_PICKUP`, `GMREDI_RESIDUAL_FLOW`, `GMREDI_RTRANSPORT`, `GMREDI_XTRANSPORT`, `GMREDI_YTRANSPORT`
  - L150 Verification experiments compiling it (1)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/gmredi@wadcheckin.md`
- L1 pkg/gmredi  (wadcheckin: ~/Documents/research/ECCO/MITgcm_wad_checkin) — `useGMREDI`
  - L13 README
  - L33 Namelist parameters
    - L34 GM_PARM01 — `GM_AdvForm`, `GM_AdvSeparate`, `GM_InMomAsStress`, `GM_isopycK`, `GM_background_K`, `GM_iso2dFile`, `GM_iso1dFile`, `GM_bol2dFile`, `GM_bol1dFile`, `GM_K3dRediFile`, `GM_K3dGMFile`, `GM_background_K3dFile`, `GM_isopycK3dFile`, `GM_taper_scheme`, `GM_maxSlope`, `GM_Kmin_horiz`, `GM_Small_Number`, `GM_slopeSqCutoff`, `GM_Scrit`, `GM_Sd`, `GM_isoFac_calcK`, `GM_facTrL2dz`, `GM_facTrL2ML`, `GM_maxTransLay`, `GM_UseBVP`, `GM_BVP_cMin`, `GM_BVP_ModeNumber`, `GM_useSubMeso`, `subMeso_Ceff`, `subMeso_invTau` (+45)
  - L111 CPP options (defaults as shipped) — `GM_READ_K3D_REDI`, `GM_READ_K3D_GM`, `GM_VISBECK_VARIABLE_K`, `GM_GEOM_VARIABLE_K`, `GM_BATES_K3D`, `GM_BATES_PASSIVE`, `GM_NON_UNITY_DIAGONAL`, `GM_EXTRA_DIAGONAL`, `GM_BOLUS_ADVEC`, `GM_BOLUS_BVP`, `ALLOW_GM_LEITH_QG`, `GM_AUTODIFF_EXCESSIVE_STORE`, `GMREDI_MASK_SLOPES`
  - L126 Headers — `GMREDI.h`, `GMREDI_OPTIONS.h`
  - L130 Routines (29) — `gmredi_calc_bates_k.F`, `gmredi_calc_diff.F`, `gmredi_calc_eigs.F`, `gmredi_calc_geom.F`, `gmredi_calc_psi_bolus.F`, `gmredi_calc_psi_bvp.F`, `gmredi_calc_qgleith.F`, `gmredi_calc_tensor.F`, `gmredi_calc_urms.F`, `gmredi_check.F`, `gmredi_diagnostics_fill.F`, `gmredi_diagnostics_impl.F`, `gmredi_diagnostics_init.F`, `gmredi_do_exch.F`, `gmredi_init_fixed.F`, `gmredi_init_varia.F`, `gmredi_mnc_init.F`, `gmredi_output.F`, `gmredi_read_pickup.F`, `gmredi_readparms.F`, `gmredi_residual_flow.F`, `gmredi_rtransport.F`, `gmredi_slope_limit.F`, `gmredi_slope_psi.F`, `gmredi_write_pickup.F`, `gmredi_xtransport.F`, `gmredi_ytransport.F`, `submeso_calc_psi.F`
  - L133 Called from outside the package — `GMREDI_CALC_DIFF`, `GMREDI_CALC_TENSOR`, `GMREDI_CALC_TENSOR_DUMMY`, `GMREDI_DO_EXCH`, `GMREDI_DIAGNOSTICS_IMPL`, `GMREDI_OUTPUT`, `GMREDI_CHECK`, `GMREDI_INIT_FIXED`, `GMREDI_INIT_VARIA`, `GMREDI_READPARMS`, `GMREDI_WRITE_PICKUP`, `GMREDI_RESIDUAL_FLOW`, `GMREDI_RTRANSPORT`, `GMREDI_XTRANSPORT`, `GMREDI_YTRANSPORT`
  - L150 Verification experiments compiling it (1)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/grdchk.md`
- L1 pkg/grdchk — `useGRDCHK`
  - L11 Namelist parameters
    - L12 GRDCHK_NML — `grdchk_eps`, `useCentralDiff`, `iGloPos`, `jGloPos`, `kGloPos`, `iGloTile`, `jGloTile`
  - L31 Headers — `GRDCHK.h`, `GRDCHK_OPTIONS.h`
  - L35 Routines (14) — `grdchk_check.F`, `grdchk_ctrl_fname.F`, `grdchk_get_mask.F`, `grdchk_get_obcs_mask.F`, `grdchk_get_position.F`, `grdchk_getadxx.F`, `grdchk_getxx.F`, `grdchk_init.F`, `grdchk_loc.F`, `grdchk_main.F`, `grdchk_print.F`, `grdchk_readparms.F`, `grdchk_setxx.F`, `grdchk_summary.F`
  - L38 Called from outside the package — `GRDCHK_CHECK`, `GRDCHK_READPARMS`, `GRDCHK_MAIN`
  - L44 Verification experiments compiling it (16) — `bottom_ctrl_5x5`, `global_oce_biogeo_bling`, `global_oce_latlon`, `global_ocean.90x40x15`, `global_ocean.cs32x15`, `halfpipe_streamice`, `lab_sea`, `obcs_ctrl`, `offline_exf_seaice`, `tutorial_dic_adjoffline`, `tutorial_global_oce_biogeo`, `tutorial_global_oce_optim`, `tutorial_tracer_adjsens`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/gridalt.md`
- L1 pkg/gridalt — `useGRIDALT`
  - L8 Headers — `GRIDALT_OPTIONS.h`, `gridalt_mapping.h`
  - L12 Routines (6) — `dyn2phys.F`, `gridalt_diagnostics_init.F`, `gridalt_initialise.F`, `gridalt_update.F`, `make_phys_grid.F`, `phys2dyn.F`
  - L15 Called from outside the package — `GRIDALT_UPDATE`, `GRIDALT_INITIALISE`
  - L24 Verification experiments compiling it (3)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/icefront.md`
- L1 pkg/icefront — `useICEFRONT`
  - L9 Namelist parameters
    - L10 ICEFRONT_PARM01 — `rhoIcefront`, `ICEFRONTlatentHeat`, `ICEFRONTHeatCapacity_Cp`, `ICEFRONTthetaSurface`, `applyIcefrontTendT`, `applyIcefrontTendS`, `ICEFRONTdepthFile`, `ICEFRONTlengthFile`
    - L20 ICEFRONT_EXF_PARM02 — `SGRunOffFile`, `SGRunOffperiod`, `SGRunOffStartTime`, `SGRunOffstartdate1`, `SGRunOffstartdate2`, `SGRunOffconst`, `SGRunOff_inscal`, `SGRunOff_remov_intercept`, `SGRunOff_remov_slope`
  - L31 Headers — `ICEFRONT.h`, `ICEFRONT_OPTIONS.h`
  - L35 Routines (8) — `icefront_check.F`, `icefront_diagnostics_init.F`, `icefront_init_fixed.F`, `icefront_init_varia.F`, `icefront_readparms.F`, `icefront_tendency_apply.F`, `icefront_thermodynamics.F`
  - L38 Called from outside the package — `ICEFRONT_TENDENCY_APPLY_S`, `ICEFRONT_TENDENCY_APPLY_T`, `ICEFRONT_THERMODYNAMICS`, `ICEFRONT_CHECK`, `ICEFRONT_INIT_FIXED`, `ICEFRONT_INIT_VARIA`, `ICEFRONT_READPARMS`
  - L49 Verification experiments compiling it (1)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/kl10.md`
- L1 pkg/kl10 — `useKL10`
  - L9 Namelist parameters
    - L10 KL10_PARM01 — `KLviscMax`, `KLdumpFreq`, `KLwriteState`, `KLtaveFreq`
  - L16 Headers — `KL10.h`, `KL10_DIAGS.h`, `KL10_OPTIONS.h`
  - L21 Routines (9) — `kl10_calc.F`, `kl10_calc_diff.F`, `kl10_calc_visc.F`, `kl10_check.F`, `kl10_diagnostics_init.F`, `kl10_init_fixed.F`, `kl10_init_varia.F`, `kl10_output.F`, `kl10_readparms.F`
  - L24 Called from outside the package — `KL10_CALC_DIFF`, `KL10_CALC_VISC`, `KL10_CALC`, `KL10_OUTPUT`, `KL10_CHECK`, `KL10_INIT_FIXED`, `KL10_INIT_VARIA`, `KL10_READPARMS`
  - L34 Verification experiments compiling it (1) — `internal_wave`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/kpp.md`
- L1 pkg/kpp — `useKPP`
  - L11 Namelist parameters
    - L12 KPP_PARM01 — `kpp_freq`, `kpp_dumpFreq`, `kpp_taveFreq`, `KPPwriteState`, `KPP_ghatUseTotalDiffus`, `KPPuseDoubleDiff`, `LimitHblStable`, `KPPuseSWfrac3D`, `minKPPhbl`, `dB_dz`, `num_v_smooth_Ri`, `KPPmixingMaps`, `num_v_smooth_BV`, `num_z_smooth_sh`, `num_m_smooth_sh`
  - L62 CPP options (defaults as shipped) — `KPP_SMOOTH_SHSQ`, `KPP_SMOOTH_DVSQ`, `KPP_SMOOTH_DBLOC`, `KPP_SMOOTH_DENS`, `KPP_SMOOTH_VISC`, `KPP_SMOOTH_DIFF`, `KPP_ESTIMATE_UREF`, `KPP_DO_NOT_MATCH_DIFFUSIVITIES`, `KPP_DO_NOT_MATCH_DERIVATIVES`, `KPP_SMOOTH_REGULARISATION`, `KPP_SCALE_SHEARMIXING`, `KPP_GHAT`, `EXCLUDE_KPP_SHEAR_MIX`, `EXCLUDE_KPP_DOUBLEDIFF`, `KPP_AUTODIFF_EXCESSIVE_STORE`, `ALLOW_KPP_VERTICALLY_SMOOTH`
  - L81 Headers — `KPP.h`, `KPP_OPTIONS.h`, `KPP_PARAMS.h`
  - L86 Routines (27) — `kpp_calc.F`, `kpp_calc_diff_ptr.F`, `kpp_calc_diff_s.F`, `kpp_calc_diff_t.F`, `kpp_calc_visc.F`, `kpp_check.F`, `kpp_diagnostics_init.F`, `kpp_do_exch.F`, `kpp_forcing_surf.F`, `kpp_init_fixed.F`, `kpp_init_varia.F`, `kpp_output.F`, `kpp_readparms.F`, `kpp_routines.F`, `kpp_transport_ptr.F`, `kpp_transport_s.F`, `kpp_transport_t.F`
  - L89 Called from outside the package — `KPP_CALC_DIFF_PTR`, `KPP_CALC_DIFF_S`, `KPP_CALC_DIFF_T`, `KPP_CALC_VISC`, `KPP_CALC`, `KPP_CALC_DUMMY`, `KPP_DO_EXCH`, `KPP_OUTPUT`, `KPP_CHECK`, `KPP_INIT_FIXED`, `KPP_INIT_VARIA`, `KPP_READPARMS`, `KPP_TRANSPORT_PTR`, `KPP_TRANSPORT_S`, `KPP_TRANSPORT_T`
  - L106 Verification experiments compiling it (5) — `lab_sea`, `seaice_obcs`, `tutorial_tracer_adjsens`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/kpp@wad.md`
- L1 pkg/kpp  (wad: ~/Documents/research/ECCO/wetting_drying/MITgcm) — `useKPP`
  - L13 Namelist parameters
    - L14 KPP_PARM01 — `kpp_freq`, `kpp_dumpFreq`, `kpp_taveFreq`, `KPPwriteState`, `KPP_ghatUseTotalDiffus`, `KPPuseDoubleDiff`, `LimitHblStable`, `KPPuseSWfrac3D`, `minKPPhbl`, `dB_dz`, `num_v_smooth_Ri`, `KPPmixingMaps`, `num_v_smooth_BV`, `num_z_smooth_sh`, `num_m_smooth_sh`
  - L64 CPP options (defaults as shipped) — `KPP_SMOOTH_SHSQ`, `KPP_SMOOTH_DVSQ`, `KPP_SMOOTH_DBLOC`, `KPP_SMOOTH_DENS`, `KPP_SMOOTH_VISC`, `KPP_SMOOTH_DIFF`, `KPP_ESTIMATE_UREF`, `KPP_DO_NOT_MATCH_DIFFUSIVITIES`, `KPP_DO_NOT_MATCH_DERIVATIVES`, `KPP_SMOOTH_REGULARISATION`, `KPP_SCALE_SHEARMIXING`, `KPP_GHAT`, `EXCLUDE_KPP_SHEAR_MIX`, `EXCLUDE_KPP_DOUBLEDIFF`, `KPP_AUTODIFF_EXCESSIVE_STORE`, `ALLOW_KPP_VERTICALLY_SMOOTH`
  - L83 Headers — `KPP.h`, `KPP_OPTIONS.h`, `KPP_PARAMS.h`
  - L88 Routines (27) — `kpp_calc.F`, `kpp_calc_diff_ptr.F`, `kpp_calc_diff_s.F`, `kpp_calc_diff_t.F`, `kpp_calc_visc.F`, `kpp_check.F`, `kpp_diagnostics_init.F`, `kpp_do_exch.F`, `kpp_forcing_surf.F`, `kpp_init_fixed.F`, `kpp_init_varia.F`, `kpp_output.F`, `kpp_readparms.F`, `kpp_routines.F`, `kpp_transport_ptr.F`, `kpp_transport_s.F`, `kpp_transport_t.F`
  - L91 Called from outside the package — `KPP_CALC_DIFF_PTR`, `KPP_CALC_DIFF_S`, `KPP_CALC_DIFF_T`, `KPP_CALC_VISC`, `KPP_CALC`, `KPP_CALC_DUMMY`, `KPP_DO_EXCH`, `KPP_OUTPUT`, `KPP_CHECK`, `KPP_INIT_FIXED`, `KPP_INIT_VARIA`, `KPP_READPARMS`, `KPP_TRANSPORT_PTR`, `KPP_TRANSPORT_S`, `KPP_TRANSPORT_T`
  - L108 Verification experiments compiling it (4)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/kpp@wadcheckin.md`
- L1 pkg/kpp  (wadcheckin: ~/Documents/research/ECCO/MITgcm_wad_checkin) — `useKPP`
  - L13 Namelist parameters
    - L14 KPP_PARM01 — `kpp_freq`, `kpp_dumpFreq`, `kpp_taveFreq`, `KPPwriteState`, `KPP_ghatUseTotalDiffus`, `KPPuseDoubleDiff`, `LimitHblStable`, `KPPuseSWfrac3D`, `minKPPhbl`, `dB_dz`, `num_v_smooth_Ri`, `KPPmixingMaps`, `num_v_smooth_BV`, `num_z_smooth_sh`, `num_m_smooth_sh`
  - L64 CPP options (defaults as shipped) — `KPP_SMOOTH_SHSQ`, `KPP_SMOOTH_DVSQ`, `KPP_SMOOTH_DBLOC`, `KPP_SMOOTH_DENS`, `KPP_SMOOTH_VISC`, `KPP_SMOOTH_DIFF`, `KPP_ESTIMATE_UREF`, `KPP_DO_NOT_MATCH_DIFFUSIVITIES`, `KPP_DO_NOT_MATCH_DERIVATIVES`, `KPP_SMOOTH_REGULARISATION`, `KPP_SCALE_SHEARMIXING`, `KPP_GHAT`, `EXCLUDE_KPP_SHEAR_MIX`, `EXCLUDE_KPP_DOUBLEDIFF`, `KPP_AUTODIFF_EXCESSIVE_STORE`, `ALLOW_KPP_VERTICALLY_SMOOTH`
  - L83 Headers — `KPP.h`, `KPP_OPTIONS.h`, `KPP_PARAMS.h`
  - L88 Routines (27) — `kpp_calc.F`, `kpp_calc_diff_ptr.F`, `kpp_calc_diff_s.F`, `kpp_calc_diff_t.F`, `kpp_calc_visc.F`, `kpp_check.F`, `kpp_diagnostics_init.F`, `kpp_do_exch.F`, `kpp_forcing_surf.F`, `kpp_init_fixed.F`, `kpp_init_varia.F`, `kpp_output.F`, `kpp_readparms.F`, `kpp_routines.F`, `kpp_transport_ptr.F`, `kpp_transport_s.F`, `kpp_transport_t.F`
  - L91 Called from outside the package — `KPP_CALC_DIFF_PTR`, `KPP_CALC_DIFF_S`, `KPP_CALC_DIFF_T`, `KPP_CALC_VISC`, `KPP_CALC`, `KPP_CALC_DUMMY`, `KPP_DO_EXCH`, `KPP_OUTPUT`, `KPP_CHECK`, `KPP_INIT_FIXED`, `KPP_INIT_VARIA`, `KPP_READPARMS`, `KPP_TRANSPORT_PTR`, `KPP_TRANSPORT_S`, `KPP_TRANSPORT_T`
  - L108 Verification experiments compiling it (3)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/land.md`
- L1 pkg/land — `useLAND`
  - L9 Namelist parameters
    - L10 LAND_MODEL_PAR — `land_calc_grT`, `land_calc_grW`, `land_impl_grT`, `land_calc_snow`, `land_calc_alb`, `land_oldPickup`, `land_grT_iniFile`, `land_grW_iniFile`, `land_snow_iniFile`, `land_deltaT`, `land_taveFreq`, `land_diagFreq`, `land_monFreq`, `land_dzF`, `land_timeave_mnc`, `land_snapshot_mnc`, `land_mon_mnc`, `land_pickup_write_mnc`, `land_pickup_read_mnc`
    - L30 LAND_PHYS_PAR — `land_grdLambda`, `land_heatCs`, `land_CpWater`, `land_wTauDiff`, `land_waterCap`, `land_fractRunOff`, `land_rhoLiqW`, `land_rhoSnow`, `land_Lfreez`, `land_hMaxSnow`, `diffKsnow`, `timeSnowAge`, `hNewSnowAge`, `albColdSnow`, `albWarmSnow`, `tempSnowAlbL`, `albOldSnow`, `hAlbSnow`
  - L50 CPP options (defaults as shipped) — `LAND_DEBUG`, `LAND_OLD_VERSION`
  - L54 Headers — `LAND_OPTIONS.h`, `LAND_PARAMS.h`, `LAND_SIZE.h`, `LAND_VARS.h`
  - L60 Routines (15) — `land_albedo.F`, `land_check.F`, `land_diagnostics_init.F`, `land_diagnostics_state.F`, `land_do_diags.F`, `land_impl_temp.F`, `land_ini_vars.F`, `land_init_fixed.F`, `land_mnc_init.F`, `land_monitor.F`, `land_output.F`, `land_read_pickup.F`, `land_readparms.F`, `land_stepfwd.F`, `land_write_pickup.F`
  - L63 Called from outside the package — `LAND_DIAGNOSTICS_STATE`, `LAND_OUTPUT`, `LAND_CHECK`, `LAND_INIT_FIXED`, `LAND_INI_VARS`, `LAND_READPARMS`, `LAND_WRITE_PICKUP`, `LAND_DO_DIAGS`, `LAND_STEPFWD`, `LAND_ALBEDO`, `LAND_IMPL_TEMP`
  - L76 Verification experiments compiling it (2) — `aim.5l_cs`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/layers.md`
- L1 pkg/layers — `useLAYERS`
  - L10 Namelist parameters
    - L11 LAYERS_PARM01 — `layers_G`, `layers_taveFreq`, `layers_diagFreq`, `LAYER_nb`, `layers_kref`, `useBOLUS`, `layers_bolus`, `layers_name`, `layers_bounds`, `layers_krho`
  - L23 CPP options (defaults as shipped) — `LAYERS_UFLUX`, `LAYERS_VFLUX`, `LAYERS_THICKNESS`, `LAYERS_THERMODYNAMICS`, `LAYERS_FINEGRID_DIAPYCNAL`, `LAYERS_MNC`, `LAYERS_PRHO_REF`, `LAYERS_MSE`
  - L33 Headers — `LAYERS.h`, `LAYERS_OPTIONS.h`, `LAYERS_P2SHARE.h`, `LAYERS_SIZE.h`
  - L39 Routines (17) — `layers_calc.F`, `layers_calc_divergence.F`, `layers_check.F`, `layers_diagnostics_init.F`, `layers_diapycnal.F`, `layers_fill.F`, `layers_fluxcalc.F`, `layers_init_fixed.F`, `layers_init_varia.F`, `layers_locate.F`, `layers_mnc_init.F`, `layers_output.F`, `layers_readparms.F`, `layers_thermodynamics.F`, `layers_wsurf_tr.F`
  - L42 Called from outside the package — `LAYERS_FILL`, `LAYERS_CALC`, `LAYERS_CHECK`, `LAYERS_INIT_FIXED`, `LAYERS_INIT_VARIA`, `LAYERS_READPARMS`, `LAYERS_WSURF_TR`
  - L57 Verification experiments compiling it (3) — `cfc_example`, `tutorial_reentrant_channel`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/longstep.md`
- L1 pkg/longstep — `useLONGSTEP`
  - L9 README
  - L43 Namelist parameters
    - L44 LONGSTEP_PARM01 — `LS_nIter`, `LS_whenToSample`
  - L48 Headers — `LONGSTEP.h`, `LONGSTEP_OPTIONS.h`, `LONGSTEP_PARAMS.h`
  - L53 Routines (17) — `longstep_average.F`, `longstep_average_3d.F`, `longstep_average_3d_fac.F`, `longstep_check.F`, `longstep_check_iters.F`, `longstep_correction_step.F`, `longstep_diagnostics_init.F`, `longstep_fill_3d.F`, `longstep_fill_3d_fac.F`, `longstep_fill_3d_rs.F`, `longstep_forcing_surf.F`, `longstep_init_fixed.F`, `longstep_init_varia.F`, `longstep_readparms.F`, `longstep_reset_3d.F`, `longstep_residual_flow.F`, `longstep_thermodynamics.F`
  - L56 Called from outside the package — `LONGSTEP_AVERAGE`, `LONGSTEP_THERMODYNAMICS`, `LONGSTEP_CHECK`, `LONGSTEP_INIT_FIXED`, `LONGSTEP_INIT_VARIA`, `LONGSTEP_READPARMS`, `LONGSTEP_CHECK_ITERS`
  - L65 Verification experiments compiling it (1) — `lab_sea`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/mangrove@wad.md`
- L1 pkg/mangrove  (wad: ~/Documents/research/ECCO/wetting_drying/MITgcm) — `useMANGROVE`
  - L10 Namelist parameters
    - L11 MANGROVE_PARM01 — `mgNsp`, `mgCd`, `mgTrunkN`, `mgTrunkD`, `mgRootA0`, `mgRootH`, `mgRootP`, `mgRootD`, `mgCoverFile`, `mgNepfAlpha`, `mgMaxKz`, `mgBedShelter`, `mgWaveAtt`, `mgDistFile`
  - L27 Headers — `MANGROVE.h`, `MANGROVE_OPTIONS.h`
  - L31 Routines (10) — `mangrove_calc_visc.F`, `mangrove_check.F`, `mangrove_diagnostics_init.F`, `mangrove_drag.F`, `mangrove_init_fixed.F`, `mangrove_init_varia.F`, `mangrove_readparms.F`, `mangrove_wave_att.F`
  - L34 Called from outside the package — `MANGROVE_CALC_DIFF`, `MANGROVE_CALC_VISC`, `MANGROVE_DRAG`, `MANGROVE_CHECK`, `MANGROVE_INIT_FIXED`, `MANGROVE_INIT_VARIA`, `MANGROVE_READPARMS`, `MANGROVE_WAVE_ATT`
  - L44 Verification experiments compiling it (1)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/matrix.md`
- L1 pkg/matrix — `useMATRIX`
  - L9 Namelist parameters
    - L10 MATRIX_PARM01 — `expMatrixWriteTime`, `impMatrixWriteTime`
  - L14 Headers — `MATRIX.h`, `MATRIX_OPTIONS.h`
  - L18 Routines (7) — `matrix_init_varia.F`, `matrix_output.F`, `matrix_readparms.F`, `matrix_store_tendency.F`, `matrix_write_grid.F`, `matrix_write_tendency.F`
  - L21 Called from outside the package — `MATRIX_OUTPUT`, `MATRIX_INIT_VARIA`, `MATRIX_READPARMS`, `MATRIX_STORE_TENDENCY_IMP`, `MATRIX_STORE_TENDENCY_EXP`
  - L28 Verification experiments compiling it (1) — `matrix_example`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/mdsio.md`
- L1 pkg/mdsio
  - L9 CPP options (defaults as shipped) — `SAFE_IO`, `_NEW_STATUS`, `_OLD_STATUS`, `ALLOW_WHIO`, `ALLOW_WHIO_3D`, `EXCLUDE_WHIO_GLOBUFF_2D`, `INCLUDE_WHIO_GLOBUFF_3D`
  - L20 Headers — `MDSIO_BUFF_3D.h`, `MDSIO_BUFF_WH.h`, `MDSIO_OPTIONS.h`
  - L25 Routines (56) — `mdsio_buffertorl.F`, `mdsio_buffertors.F`, `mdsio_check4file.F`, `mdsio_facef_read.F`, `mdsio_gl.F`, `mdsio_gl_slice.F`, `mdsio_pass_r4torl.F`, `mdsio_pass_r4tors.F`, `mdsio_pass_r8torl.F`, `mdsio_pass_r8tors.F`, `mdsio_rd_rec_rl.F`, `mdsio_rd_rec_rs.F`, `mdsio_read_field.F`, `mdsio_read_meta.F`, `mdsio_read_section.F`, `mdsio_read_tape.F`, `mdsio_read_whalos.F`, `mdsio_readvec_loc.F`, `mdsio_rw_field.F`, `mdsio_rw_slice.F`, `mdsio_seg4torl.F`, `mdsio_seg4tors.F`, `mdsio_seg8torl.F`, `mdsio_seg8tors.F`, `mdsio_segxtorx_2d.F`, `mdsio_wr_metafiles.F`, `mdsio_wr_rec_rl.F`, `mdsio_wr_rec_rs.F`, `mdsio_write_field.F`, `mdsio_write_meta.F` (+5)
  - L28 Called from outside the package — `MDS_FACEF_READ_RS`, `MDS_WR_METAFILES`, `MDS_READVEC_LOC`, `MDS_WRITEVEC_LOC`, `MDS_READ_FIELD`, `MDS_WRITE_FIELD`, `MDS_READ_SEC_XZ`, `MDS_READ_SEC_YZ`, `MDS_WRITE_SEC_XZ`, `MDS_WRITE_SEC_YZ`, `MDS_READ_TAPE`, `MDS_READ_WHALOS`, `MDS_WRITE_TAPE`, `MDS_WRITE_WHALOS`, `MDS_CHECK4FILE`, `MDSREADFIELD_2D_GL`, `MDSREADFIELD_3D_GL`, `MDS_PASS_R8TORL`, `MDSREADFIELD_XZ_GL`, `MDSREADFIELD_YZ_GL`, `MDSWRITEFIELD_2D_GL`, `MDSWRITEFIELD_3D_GL`, `MDSWRITEFIELD_XZ_GL`, `MDSWRITEFIELD_YZ_GL`, `MDS_READ_META`, `MDS_WRITELOCAL`
  - L130 Verification experiments compiling it (58) — `advect_cs`, `advect_xz`, `aim.5l_Equatorial_Channel`, `aim.5l_LatLon`, `aim.5l_cs`, `atm_gray`, `bottom_ctrl_5x5`, `cfc_example`, `cheapAML_box`, `deep_anelastic`, `front_relax`, `global_oce_biogeo_bling`, `global_oce_latlon`, `global_ocean.90x40x15`, `global_ocean.cs32x15`, `halfpipe_streamice`, `ideal_2D_oce`, `internal_wave`, `inverted_barometer`, `lab_sea`, `matrix_example`, `obcs_ctrl`, `offline_exf_seaice`, `seaice_itd`, `seaice_obcs`, `shelfice_2d_remesh`, `short_surf_wave`, `so_box_biogeo`, `tutorial_advection_in_gyre`, `tutorial_baroclinic_gyre` (+12)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/mnc.md`
- L1 pkg/mnc — `useMNC`
  - L8 README
  - L42 Namelist parameters
    - L43 MNC_01 — `mnc_use_indir`, `mnc_use_outdir`, `mnc_outdir_date`, `mnc_outdir_num`, `mnc_use_name_ni0`, `mnc_echo_gvtypes`, `pickup_write_mnc`, `pickup_read_mnc`, `timeave_mnc`, `snapshot_mnc`, `monitor_mnc`, `autodiff_mnc`, `writegrid_mnc`, `readgrid_mnc`, `mnc_outdir_str`, `mnc_indir_str`, `mnc_max_fsize`, `mnc_filefreq`, `mnc_read_bathy`, `mnc_read_salt`, `mnc_read_theta`
  - L66 CPP options (defaults as shipped) — `MNC_DEF_FMNC`, `MNC_DEF_TMNC`
  - L70 Headers — `MNC_BUFF.h`, `MNC_COMMON.h`, `MNC_OPTIONS.h`, `MNC_PARAMS.h`, `MNC_SIZE.h`
  - L77 Routines (93) — `MNC_CW_READWRITE_I.F`, `MNC_CW_READWRITE_RL.F`, `MNC_CW_READWRITE_RS.F`, `mnc_cw_citer.F`, `mnc_cw_cvars.F`, `mnc_cw_init.F`, `mnc_cw_missingvals.F`, `mnc_cw_model_attr.F`, `mnc_cw_udim.F`, `mnc_cw_write_grid_info.F`, `mnc_cwrapper.F`, `mnc_dim.F`, `mnc_dump.F`, `mnc_file.F`, `mnc_grid.F`, `mnc_init.F`, `mnc_readparms.F`, `mnc_update_time.F`, `mnc_utils.F`, `mnc_var.F`
  - L80 Called from outside the package — `MNC_UPDATE_TIME`, `MNC_CW_SET_UDIM`, `MNC_CW_RS_R`, `MNC_CW_SET_CITER`, `MNC_FILE_CLOSE_ALL_MATCHING`, `MNC_CW_ADD_VNAME`, `MNC_CW_DEL_VNAME`, `MNC_CW_ADD_VATTR_TEXT`, `MNC_CW_INIT`, `MNC_INIT`, `MNC_CW_RL_R`, `MNC_READPARMS`, `MNC_FILE_CLOSE_ALL`, `MNC_CW_RS_W`, `MNC_CW_WRITE_GRID_COORD`, `MNC_CW_I_W_S`, `MNC_CW_RL_W`, `MNC_CW_RL_W_S`, `MNC_CW_ADD_GNAME`, `MNC_CW_DEL_GNAME`, `MNC_CW_VATTR_MISSING`, `MNC_CW_I_W`, `MNC_CW_APPEND_VNAME`, `MNC_CW_RL_W_OFFSET`
  - L334 Verification experiments compiling it (21) — `aim.5l_cs`, `cheapAML_box`, `global_oce_biogeo_bling`, `global_oce_latlon`, `global_ocean.90x40x15`, `global_ocean.cs32x15`, `internal_wave`, `inverted_barometer`, `lab_sea`, `tutorial_advection_in_gyre`, `tutorial_baroclinic_gyre`, `tutorial_dic_adjoffline`, `tutorial_global_oce_biogeo`, `tutorial_global_oce_latlon`, `tutorial_held_suarez_cs`, `tutorial_rotating_tank`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/mnc@bbl.md`
- L1 pkg/mnc  (bbl: ~/Documents/research/ECCO/BBL/MITgcm) — `useMNC`
  - L10 README
  - L44 Namelist parameters
    - L45 MNC_01 — `mnc_use_indir`, `mnc_use_outdir`, `mnc_outdir_date`, `mnc_outdir_num`, `mnc_use_name_ni0`, `mnc_echo_gvtypes`, `pickup_write_mnc`, `pickup_read_mnc`, `timeave_mnc`, `snapshot_mnc`, `monitor_mnc`, `autodiff_mnc`, `writegrid_mnc`, `readgrid_mnc`, `mnc_outdir_str`, `mnc_indir_str`, `mnc_max_fsize`, `mnc_filefreq`, `mnc_read_bathy`, `mnc_read_salt`, `mnc_read_theta`
  - L68 CPP options (defaults as shipped) — `MNC_DEF_FMNC`, `MNC_DEF_TMNC`
  - L72 Headers — `MNC_BUFF.h`, `MNC_COMMON.h`, `MNC_ID_HEADER.h`, `MNC_OPTIONS.h`, `MNC_PARAMS.h`, `MNC_SIZE.h`
  - L80 Routines (93) — `MNC_CW_READWRITE_I.F`, `MNC_CW_READWRITE_RL.F`, `MNC_CW_READWRITE_RS.F`, `mnc_cw_citer.F`, `mnc_cw_cvars.F`, `mnc_cw_init.F`, `mnc_cw_missingvals.F`, `mnc_cw_model_attr.F`, `mnc_cw_udim.F`, `mnc_cw_write_grid_info.F`, `mnc_cwrapper.F`, `mnc_dim.F`, `mnc_dump.F`, `mnc_file.F`, `mnc_grid.F`, `mnc_init.F`, `mnc_readparms.F`, `mnc_update_time.F`, `mnc_utils.F`, `mnc_var.F`
  - L83 Called from outside the package — `MNC_UPDATE_TIME`, `MNC_CW_SET_UDIM`, `MNC_CW_RS_R`, `MNC_CW_SET_CITER`, `MNC_FILE_CLOSE_ALL_MATCHING`, `MNC_CW_ADD_VNAME`, `MNC_CW_DEL_VNAME`, `MNC_CW_ADD_VATTR_TEXT`, `MNC_CW_INIT`, `MNC_INIT`, `MNC_CW_RL_R`, `MNC_READPARMS`, `MNC_FILE_CLOSE_ALL`, `MNC_CW_RS_W`, `MNC_CW_WRITE_GRID_COORD`, `MNC_CW_I_W_S`, `MNC_CW_RL_W`, `MNC_CW_RL_W_S`, `MNC_CW_ADD_GNAME`, `MNC_CW_DEL_GNAME`, `MNC_CW_VATTR_MISSING`, `MNC_CW_I_W`, `MNC_CW_APPEND_VNAME`, `MNC_CW_RL_W_OFFSET`
  - L351 Verification experiments compiling it (1)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/mom_common.md`
- L1 pkg/mom_common
  - L10 CPP options (defaults as shipped) — `COSINEMETH_III`, `ISOTROPIC_COS_SCALING`, `ALLOW_LEITH_QG`, `ALLOW_SMAG_3D`, `ALLOW_3D_VISCAH`, `ALLOW_3D_VISCA4`, `ALLOW_BOTTOMDRAG_ROUGHNESS`, `MOM_USE_OLD_DEEP_VERT_ADV`, `ALLOW_MOM_TEND_EXTRA_DIAGS`
  - L21 Headers — `MOM_COMMON_OPTIONS.h`, `MOM_VISC.h`
  - L25 Routines (33) — `mom_calc_3d_strain.F`, `mom_calc_absvort3.F`, `mom_calc_hdiv.F`, `mom_calc_hfacz.F`, `mom_calc_ke.F`, `mom_calc_relvort3.F`, `mom_calc_smag_3d.F`, `mom_calc_strain.F`, `mom_calc_tension.F`, `mom_calc_visc.F`, `mom_diagnostics_init.F`, `mom_hdissip.F`, `mom_init_fixed.F`, `mom_quasihydrostatic.F`, `mom_u_botdrag_coeff.F`, `mom_u_coriolis_nh.F`, `mom_u_implicit_r.F`, `mom_u_metric_nh.F`, `mom_u_rviscflux.F`, `mom_u_sidedrag.F`, `mom_uv_smag_3d.F`, `mom_v_botdrag_coeff.F`, `mom_v_coriolis_nh.F`, `mom_v_implicit_r.F`, `mom_v_metric_nh.F`, `mom_v_rviscflux.F`, `mom_v_sidedrag.F`, `mom_visc_qgl_limit.F`, `mom_visc_qgl_stretch.F`, `mom_w_coriolis_nh.F` (+3)
  - L28 Called from outside the package — `MOM_W_CORIOLIS_NH`, `MOM_W_METRIC_NH`, `MOM_W_SIDEDRAG`, `MOM_W_SMAG_3D`, `MOM_QUASIHYDROSTATIC`, `MOM_CALC_3D_STRAIN`, `MOM_CALC_SMAG_3D`, `MOM_UV_SMAG_3D`, `MOM_U_IMPLICIT_R`, `MOM_V_IMPLICIT_R`, `MOM_INIT_FIXED`, `MOM_CALC_HDIV`, `MOM_CALC_HFACZ`, `MOM_CALC_RELVORT3`, `MOM_VISC_QGL_LIMIT`, `MOM_VISC_QGL_STRETCH`, `MOM_CALC_KE`, `MOM_CALC_STRAIN`, `MOM_CALC_TENSION`, `MOM_CALC_VISC`, `MOM_U_BOTDRAG_COEFF`, `MOM_U_CORIOLIS_NH`, `MOM_U_METRIC_NH`, `MOM_U_RVISCFLUX`, `MOM_U_SIDEDRAG`, `MOM_V_BOTDRAG_COEFF`, `MOM_V_CORIOLIS_NH`, `MOM_V_METRIC_NH`, `MOM_V_RVISCFLUX`, `MOM_V_SIDEDRAG` (+2)
  - L92 Verification experiments compiling it (54) — `advect_cs`, `aim.5l_Equatorial_Channel`, `aim.5l_LatLon`, `aim.5l_cs`, `atm_gray`, `bottom_ctrl_5x5`, `cfc_example`, `cheapAML_box`, `deep_anelastic`, `front_relax`, `global_oce_biogeo_bling`, `global_oce_latlon`, `global_ocean.90x40x15`, `global_ocean.cs32x15`, `ideal_2D_oce`, `internal_wave`, `inverted_barometer`, `lab_sea`, `matrix_example`, `obcs_ctrl`, `offline_exf_seaice`, `seaice_itd`, `seaice_obcs`, `shelfice_2d_remesh`, `short_surf_wave`, `so_box_biogeo`, `tutorial_advection_in_gyre`, `tutorial_baroclinic_gyre`, `tutorial_deep_convection`, `tutorial_global_oce_biogeo` (+8)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/mom_fluxform.md`
- L1 pkg/mom_fluxform — `useMOM_FLUXFORM`
  - L11 CPP options (defaults as shipped) — `MOM_BOUNDARY_CONSERVE`
  - L14 Headers — `MOM_FLUXFORM.h`, `MOM_FLUXFORM_OPTIONS.h`
  - L18 Routines (21) — `mom_calc_rtrans.F`, `mom_fluxform.F`, `mom_u_adv_uu.F`, `mom_u_adv_vu.F`, `mom_u_adv_wu.F`, `mom_u_coriolis.F`, `mom_u_del2u.F`, `mom_u_metric_cylinder.F`, `mom_u_metric_sphere.F`, `mom_u_xviscflux.F`, `mom_u_yviscflux.F`, `mom_uv_boundary.F`, `mom_v_adv_uv.F`, `mom_v_adv_vv.F`, `mom_v_adv_wv.F`, `mom_v_coriolis.F`, `mom_v_del2v.F`, `mom_v_metric_cylinder.F`, `mom_v_metric_sphere.F`, `mom_v_xviscflux.F`, `mom_v_yviscflux.F`
  - L21 Called from outside the package — `MOM_FLUXFORM`
  - L24 Verification experiments compiling it (48) — `advect_cs`, `aim.5l_Equatorial_Channel`, `aim.5l_LatLon`, `atm_gray`, `bottom_ctrl_5x5`, `cfc_example`, `cheapAML_box`, `deep_anelastic`, `front_relax`, `global_oce_biogeo_bling`, `global_oce_latlon`, `global_ocean.90x40x15`, `global_ocean.cs32x15`, `ideal_2D_oce`, `internal_wave`, `inverted_barometer`, `lab_sea`, `matrix_example`, `offline_exf_seaice`, `seaice_itd`, `seaice_obcs`, `shelfice_2d_remesh`, `short_surf_wave`, `so_box_biogeo`, `tutorial_advection_in_gyre`, `tutorial_baroclinic_gyre`, `tutorial_deep_convection`, `tutorial_global_oce_biogeo`, `tutorial_global_oce_in_p`, `tutorial_global_oce_latlon` (+6)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/mom_vecinv.md`
- L1 pkg/mom_vecinv — `useMOM_VECINV`
  - L11 CPP options (defaults as shipped) — `MOM_VI_ORIGINAL_VISCA4`
  - L14 Headers — `MOM_VECINV_OPTIONS.h`
  - L17 Routines (12) — `mom_vecinv.F`, `mom_vi_coriolis.F`, `mom_vi_del2uv.F`, `mom_vi_hdissip.F`, `mom_vi_u_coriolis.F`, `mom_vi_u_coriolis_c4.F`, `mom_vi_u_grad_ke.F`, `mom_vi_u_vertshear.F`, `mom_vi_v_coriolis.F`, `mom_vi_v_coriolis_c4.F`, `mom_vi_v_grad_ke.F`, `mom_vi_v_vertshear.F`
  - L20 Called from outside the package — `MOM_VECINV`, `MOM_VI_U_CORIOLIS`, `MOM_VI_U_CORIOLIS_C4`, `MOM_VI_U_GRAD_KE`, `MOM_VI_V_CORIOLIS`, `MOM_VI_V_CORIOLIS_C4`, `MOM_VI_V_GRAD_KE`, `MOM_VI_DEL2UV`
  - L30 Verification experiments compiling it (53) — `advect_cs`, `aim.5l_Equatorial_Channel`, `aim.5l_LatLon`, `aim.5l_cs`, `atm_gray`, `bottom_ctrl_5x5`, `cfc_example`, `cheapAML_box`, `deep_anelastic`, `front_relax`, `global_oce_biogeo_bling`, `global_oce_latlon`, `global_ocean.90x40x15`, `global_ocean.cs32x15`, `ideal_2D_oce`, `internal_wave`, `inverted_barometer`, `lab_sea`, `matrix_example`, `obcs_ctrl`, `offline_exf_seaice`, `seaice_itd`, `seaice_obcs`, `shelfice_2d_remesh`, `short_surf_wave`, `so_box_biogeo`, `tutorial_advection_in_gyre`, `tutorial_baroclinic_gyre`, `tutorial_deep_convection`, `tutorial_global_oce_biogeo` (+7)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/mom_vecinv@wad.md`
- L1 pkg/mom_vecinv  (wad: ~/Documents/research/ECCO/wetting_drying/MITgcm) — `useMOM_VECINV`
  - L13 CPP options (defaults as shipped) — `MOM_VI_ORIGINAL_VISCA4`
  - L16 Headers — `MOM_VECINV_OPTIONS.h`
  - L19 Routines (12) — `mom_vecinv.F`, `mom_vi_coriolis.F`, `mom_vi_del2uv.F`, `mom_vi_hdissip.F`, `mom_vi_u_coriolis.F`, `mom_vi_u_coriolis_c4.F`, `mom_vi_u_grad_ke.F`, `mom_vi_u_vertshear.F`, `mom_vi_v_coriolis.F`, `mom_vi_v_coriolis_c4.F`, `mom_vi_v_grad_ke.F`, `mom_vi_v_vertshear.F`
  - L22 Called from outside the package — `MOM_VECINV`, `MOM_VI_U_CORIOLIS`, `MOM_VI_U_CORIOLIS_C4`, `MOM_VI_U_GRAD_KE`, `MOM_VI_V_CORIOLIS`, `MOM_VI_V_CORIOLIS_C4`, `MOM_VI_V_GRAD_KE`, `MOM_VI_DEL2UV`
  - L32 Verification experiments compiling it (8)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/mom_vecinv@wadcheckin.md`
- L1 pkg/mom_vecinv  (wadcheckin: ~/Documents/research/ECCO/MITgcm_wad_checkin) — `useMOM_VECINV`
  - L13 CPP options (defaults as shipped) — `MOM_VI_ORIGINAL_VISCA4`
  - L16 Headers — `MOM_VECINV_OPTIONS.h`
  - L19 Routines (12) — `mom_vecinv.F`, `mom_vi_coriolis.F`, `mom_vi_del2uv.F`, `mom_vi_hdissip.F`, `mom_vi_u_coriolis.F`, `mom_vi_u_coriolis_c4.F`, `mom_vi_u_grad_ke.F`, `mom_vi_u_vertshear.F`, `mom_vi_v_coriolis.F`, `mom_vi_v_coriolis_c4.F`, `mom_vi_v_grad_ke.F`, `mom_vi_v_vertshear.F`
  - L22 Called from outside the package — `MOM_VECINV`, `MOM_VI_U_CORIOLIS`, `MOM_VI_U_CORIOLIS_C4`, `MOM_VI_U_GRAD_KE`, `MOM_VI_V_CORIOLIS`, `MOM_VI_V_CORIOLIS_C4`, `MOM_VI_V_GRAD_KE`, `MOM_VI_DEL2UV`
  - L32 Verification experiments compiling it (5)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/monitor.md`
- L1 pkg/monitor
  - L9 CPP options (defaults as shipped) — `MONITOR_TEST_HFACZ`
  - L12 Headers — `MONITOR.h`, `MONITOR_OPTIONS.h`
  - L16 Routines (29) — `mon_advcfl.F`, `mon_advcflw.F`, `mon_advcflw2.F`, `mon_calc_advcfl.F`, `mon_calc_stats_rl.F`, `mon_calc_stats_rs.F`, `mon_init.F`, `mon_ke.F`, `mon_out.F`, `mon_printstats_rl.F`, `mon_printstats_rs.F`, `mon_set_iounit.F`, `mon_set_pref.F`, `mon_solution.F`, `mon_stats_latbnd_rl.F`, `mon_stats_rl.F`, `mon_stats_rs.F`, `mon_surfcor.F`, `mon_vort3.F`, `mon_writestats_rl.F`, `mon_writestats_rs.F`, `monitor.F`, `monitor_ad.F`, `monitor_g.F`
  - L19 Called from outside the package — `MON_PRINTSTATS_RS`, `MON_SET_PREF`, `MON_INIT`, `MON_CALC_ADVCFL_GLOB`, `MON_CALC_ADVCFL_TILE`, `MON_OUT_I`, `MON_OUT_RL`, `MON_WRITESTATS_RL`, `MON_WRITESTATS_RS`, `MON_STATS_LATBND_RL`, `MON_CALC_STATS_RL`
  - L63 Verification experiments compiling it (58) — `advect_cs`, `advect_xz`, `aim.5l_Equatorial_Channel`, `aim.5l_LatLon`, `aim.5l_cs`, `atm_gray`, `bottom_ctrl_5x5`, `cfc_example`, `cheapAML_box`, `deep_anelastic`, `front_relax`, `global_oce_biogeo_bling`, `global_oce_latlon`, `global_ocean.90x40x15`, `global_ocean.cs32x15`, `halfpipe_streamice`, `ideal_2D_oce`, `internal_wave`, `inverted_barometer`, `lab_sea`, `matrix_example`, `obcs_ctrl`, `offline_exf_seaice`, `seaice_itd`, `seaice_obcs`, `shelfice_2d_remesh`, `short_surf_wave`, `so_box_biogeo`, `tutorial_advection_in_gyre`, `tutorial_baroclinic_gyre` (+12)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/my82.md`
- L1 pkg/my82 — `useMY82`
  - L9 Namelist parameters
    - L10 MY_PARM01 — `MYviscMax`, `MYdiffMax`, `MYhblScale`, `MYdumpFreq`, `MYMixingMaps`, `MYwriteState`
  - L18 CPP options (defaults as shipped) — `MY82_SMOOTH_RI`
  - L21 Headers — `MY82.h`, `MY82_OPTIONS.h`
  - L25 Routines (10) — `my82_calc.F`, `my82_calc_diff.F`, `my82_calc_visc.F`, `my82_check.F`, `my82_diagnostics_init.F`, `my82_init_fixed.F`, `my82_init_varia.F`, `my82_output.F`, `my82_readparms.F`, `my82_ri_number.F`
  - L28 Called from outside the package — `MY82_CALC_DIFF`, `MY82_CALC_VISC`, `MY82_CALC`, `MY82_OUTPUT`, `MY82_CHECK`, `MY82_INIT_FIXED`, `MY82_INIT_VARIA`, `MY82_READPARMS`
  - L38 Verification experiments compiling it (1)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/mypackage.md`
- L1 pkg/mypackage — `useMYPACKAGE`
  - L10 Namelist parameters
    - L11 MYPACKAGE_PARM01 — `myPa_MNC`, `myPa_StaV_Cgrid`, `myPa_Tend_Cgrid`, `myPa_applyTendT`, `myPa_applyTendS`, `myPa_applyTendU`, `myPa_applyTendV`, `myPa_doSwitch1`, `myPa_doSwitch2`, `myPa_index1`, `myPa_index2`, `myPa_param1`, `myPa_param2`, `myPa_string1`, `myPa_string2`, `myPa_Scal1File`, `myPa_Scal2File`, `myPa_VelUFile`, `myPa_VelVFile`, `myPa_Surf1File`, `myPa_Surf2File`
  - L34 CPP options (defaults as shipped) — `MYPACKAGE_3D_STATE`, `MYPACKAGE_2D_STATE`, `MYPACKAGE_TENDENCY`, `MYPA_SPECIAL_COMPILE_OPTION1`, `MYPA_SPECIAL_COMPILE_OPTION2`
  - L41 Headers — `MYPACKAGE.h`, `MYPACKAGE_OPTIONS.h`
  - L45 Routines (14) — `mypackage_calc_rhs.F`, `mypackage_check.F`, `mypackage_diagnostics_init.F`, `mypackage_diagnostics_state.F`, `mypackage_init_fixed.F`, `mypackage_init_varia.F`, `mypackage_mnc_init.F`, `mypackage_read_pickup.F`, `mypackage_readparms.F`, `mypackage_tendency_apply.F`, `mypackage_write_pickup.F`
  - L48 Called from outside the package — `MYPACKAGE_TENDENCY_APPLY_S`, `MYPACKAGE_TENDENCY_APPLY_T`, `MYPACKAGE_TENDENCY_APPLY_U`, `MYPACKAGE_TENDENCY_APPLY_V`, `MYPACKAGE_CALC_RHS`, `MYPACKAGE_DIAGNOSTICS_STATE`, `MYPACKAGE_CHECK`, `MYPACKAGE_INIT_FIXED`, `MYPACKAGE_INIT_VARIA`, `MYPACKAGE_READPARMS`, `MYPACKAGE_WRITE_PICKUP`
  - L65 Verification experiments compiling it (1)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/obcs.md`
- L1 pkg/obcs — `useOBCS`
  - L10 Namelist parameters
    - L11 OBCS_PARM01 — `insideOBmaskFile`, `OBNconnectFile`, `OBSconnectFile`, `OBEconnectFile`, `OBWconnectFile`, `OB_Jnorth`, `OB_Jsouth`, `OB_Ieast`, `OB_Iwest`, `OB_singleJnorth`, `OB_singleJsouth`, `OB_singleIeast`, `OB_singleIwest`, `useOrlanskiNorth`, `useOrlanskiSouth`, `useOrlanskiEast`, `useOrlanskiWest`, `useStevensNorth`, `useStevensSouth`, `useStevensEast`, `useStevensWest`, `useOBCSprescribe`, `useOBCStides`, `OBCS_tidalPeriod`, `OBCS_u1_adv_T`, `OBCS_u1_adv_S`, `OBNuFile`, `OBNvFile`, `OBNtFile`, `OBNsFile` (+100)
    - L142 OBCS_PARM02 — `cvelTimeScale`, `useFixedCEast`, `useFixedCWest`
    - L148 OBCS_PARM03 — `spongeThickness`
    - L154 OBCS_PARM04 — `TrelaxStevens`, `SrelaxStevens`, `useStevensPhaseVel`, `useStevensAdvection`
    - L159 OBCS_PARM05 — `seaiceSpongeThickness`
  - L170 CPP options (defaults as shipped) — `ALLOW_OBCS_NORTH`, `ALLOW_OBCS_SOUTH`, `ALLOW_OBCS_EAST`, `ALLOW_OBCS_WEST`, `ALLOW_ORLANSKI`, `ALLOW_OBCS_PRESCRIBE`, `ALLOW_OBCS_STEVENS`, `ALLOW_OBCS_SPONGE`, `ALLOW_OBCS_SEAICE_SPONGE`, `ALLOW_OBCS_BALANCE`, `ALLOW_OBCS_TIDES`, `OBCS_UVICE_OLD`, `OBCS_SEAICE_AVOID_CONVERGENCE`, `OBCS_SEAICE_SMOOTH_UVICE_PERP`, `OBCS_SEAICE_SMOOTH_UVICE_PAR`, `OBCS_SEAICE_COMPUTE_UVICE`, `OBCS_SEAICE_SMOOTH_EDGE`, `OBCS_AGEOS_COST_CONTRIBUTION`, `OBCS_VOLFLUX_COST_CONTRIBUTION`
  - L191 Headers — `OBCS_FIELDS.h`, `OBCS_GRID.h`, `OBCS_OPTIONS.h`, `OBCS_PARAMS.h`, `OBCS_PTRACERS.h`, `OBCS_SEAICE.h`, `ORLANSKI.h`, `obcs_ad_check_lev1_dir.h`, `obcs_ad_check_lev2_dir.h`, `obcs_ad_check_lev3_dir.h`, `obcs_ad_check_lev4_dir.h`
  - L204 Routines (70) — `obcs_add_tides.F`, `obcs_adjust.F`, `obcs_adjust_uvice.F`, `obcs_apply_eta.F`, `obcs_apply_ptracer.F`, `obcs_apply_r_star.F`, `obcs_apply_seaice.F`, `obcs_apply_surf_dr.F`, `obcs_apply_ts.F`, `obcs_apply_uv.F`, `obcs_apply_uvice.F`, `obcs_apply_w.F`, `obcs_balance_flow.F`, `obcs_calc.F`, `obcs_calc_stevens.F`, `obcs_check.F`, `obcs_check_depths.F`, `obcs_copy_tracer.F`, `obcs_copy_uv_n.F`, `obcs_cost_ageos.F`, `obcs_cost_driver.F`, `obcs_cost_final.F`, `obcs_cost_ob_e.F`, `obcs_cost_ob_n.F`, `obcs_cost_ob_s.F`, `obcs_cost_ob_w.F`, `obcs_cost_vol.F`, `obcs_cost_weights.F`, `obcs_diag_balance.F`, `obcs_exchanges.F` (+23)
  - L207 Called from outside the package — `OBCS_SPONGE_S`, `OBCS_SPONGE_T`, `OBCS_SPONGE_U`, `OBCS_SPONGE_V`, `OBCS_APPLY_R_STAR`, `OBCS_APPLY_SURF_DR`, `OBCS_ADJUST`, `OBCS_CALC`, `OBCS_OUTPUT`, `OBCS_APPLY_UV`, `OBCS_COPY_UV_N`, `OBCS_EXCHANGES`, `OBCS_SAVE_UV_N`, `OBCS_CHECK_DEPTHS`, `OBCS_APPLY_W`, `OBCS_CHECK`, `OBCS_INIT_FIXED`, `OBCS_INIT_VARIABLES`, `OBCS_READPARMS`, `OBCS_WRITE_PICKUP`, `OBCS_APPLY_TS`, `OBCS_APPLY_ETA`, `OBCS_COST_DRIVER`, `OBCS_COST_FINAL`, `OBCS_DIAG_BALANCE`, `OBCS_APPLY_PTRACER`, `OBCS_U1_ADV_TRACER`, `OBCS_APPLY_UVICE`, `OBCS_ADJUST_UVICE`, `OBCS_APPLY_SEAICE`
  - L246 Verification experiments compiling it (10) — `internal_wave`, `obcs_ctrl`, `offline_exf_seaice`, `seaice_obcs`, `shelfice_2d_remesh`, `so_box_biogeo`, `tutorial_plume_on_slope`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/obcs@seaicebc.md`
- L1 pkg/obcs  (seaicebc: ~/Documents/research/ECCO/sea_ice_BCs/MITgcm) — `useOBCS`
  - L12 Namelist parameters
    - L13 OBCS_PARM01 — `insideOBmaskFile`, `OBNconnectFile`, `OBSconnectFile`, `OBEconnectFile`, `OBWconnectFile`, `OB_Jnorth`, `OB_Jsouth`, `OB_Ieast`, `OB_Iwest`, `OB_singleJnorth`, `OB_singleJsouth`, `OB_singleIeast`, `OB_singleIwest`, `useOrlanskiNorth`, `useOrlanskiSouth`, `useOrlanskiEast`, `useOrlanskiWest`, `useStevensNorth`, `useStevensSouth`, `useStevensEast`, `useStevensWest`, `useOBCSprescribe`, `useOBCStides`, `OBCS_tidalPeriod`, `OBCS_u1_adv_T`, `OBCS_u1_adv_S`, `OBNuFile`, `OBNvFile`, `OBNtFile`, `OBNsFile` (+100)
    - L144 OBCS_PARM02 — `cvelTimeScale`, `useFixedCEast`, `useFixedCWest`
    - L150 OBCS_PARM03 — `spongeThickness`
    - L156 OBCS_PARM04 — `TrelaxStevens`, `SrelaxStevens`, `useStevensPhaseVel`, `useStevensAdvection`
    - L161 OBCS_PARM05 — `seaiceSpongeThickness`, `seaiceUVspongeWidth`, `seaiceUVspongeMax`, `seaiceBufferWidth`, `seaiceBufferDeltaMin`, `seaiceOBuvAnomTime`, `seaiceUVspongeInterp`, `seaiceSoftOB`, `seaiceSoftOBwgt`, `seaiceUrelaxFile`, `seaiceVrelaxFile`, `seaiceUVrelaxPeriod`, `seaiceUVrelaxStartTime`
  - L184 CPP options (defaults as shipped) — `ALLOW_OBCS_NORTH`, `ALLOW_OBCS_SOUTH`, `ALLOW_OBCS_EAST`, `ALLOW_OBCS_WEST`, `ALLOW_ORLANSKI`, `ALLOW_OBCS_PRESCRIBE`, `ALLOW_OBCS_STEVENS`, `ALLOW_OBCS_SPONGE`, `ALLOW_OBCS_SEAICE_SPONGE`, `ALLOW_OBCS_BALANCE`, `ALLOW_OBCS_TIDES`, `OBCS_UVICE_OLD`, `OBCS_SEAICE_AVOID_CONVERGENCE`, `OBCS_SEAICE_SMOOTH_UVICE_PERP`, `OBCS_SEAICE_SMOOTH_UVICE_PAR`, `OBCS_SEAICE_COMPUTE_UVICE`, `OBCS_SEAICE_SMOOTH_EDGE`, `OBCS_AGEOS_COST_CONTRIBUTION`, `OBCS_VOLFLUX_COST_CONTRIBUTION`
  - L205 Headers — `OBCS_FIELDS.h`, `OBCS_GRID.h`, `OBCS_OPTIONS.h`, `OBCS_PARAMS.h`, `OBCS_PTRACERS.h`, `OBCS_SEAICE.h`, `ORLANSKI.h`, `obcs_ad_check_lev1_dir.h`, `obcs_ad_check_lev2_dir.h`, `obcs_ad_check_lev3_dir.h`, `obcs_ad_check_lev4_dir.h`
  - L218 Routines (70) — `obcs_add_tides.F`, `obcs_adjust.F`, `obcs_adjust_uvice.F`, `obcs_apply_eta.F`, `obcs_apply_ptracer.F`, `obcs_apply_r_star.F`, `obcs_apply_seaice.F`, `obcs_apply_surf_dr.F`, `obcs_apply_ts.F`, `obcs_apply_uv.F`, `obcs_apply_uvice.F`, `obcs_apply_w.F`, `obcs_balance_flow.F`, `obcs_calc.F`, `obcs_calc_stevens.F`, `obcs_check.F`, `obcs_check_depths.F`, `obcs_copy_tracer.F`, `obcs_copy_uv_n.F`, `obcs_cost_ageos.F`, `obcs_cost_driver.F`, `obcs_cost_final.F`, `obcs_cost_ob_e.F`, `obcs_cost_ob_n.F`, `obcs_cost_ob_s.F`, `obcs_cost_ob_w.F`, `obcs_cost_vol.F`, `obcs_cost_weights.F`, `obcs_diag_balance.F`, `obcs_exchanges.F` (+23)
  - L221 Called from outside the package — `OBCS_SPONGE_S`, `OBCS_SPONGE_T`, `OBCS_SPONGE_U`, `OBCS_SPONGE_V`, `OBCS_APPLY_R_STAR`, `OBCS_APPLY_SURF_DR`, `OBCS_ADJUST`, `OBCS_CALC`, `OBCS_OUTPUT`, `OBCS_APPLY_UV`, `OBCS_COPY_UV_N`, `OBCS_EXCHANGES`, `OBCS_SAVE_UV_N`, `OBCS_CHECK_DEPTHS`, `OBCS_APPLY_W`, `OBCS_CHECK`, `OBCS_INIT_FIXED`, `OBCS_INIT_VARIABLES`, `OBCS_READPARMS`, `OBCS_WRITE_PICKUP`, `OBCS_APPLY_TS`, `OBCS_APPLY_ETA`, `OBCS_COST_DRIVER`, `OBCS_COST_FINAL`, `OBCS_DIAG_BALANCE`, `OBCS_APPLY_PTRACER`, `OBCS_U1_ADV_TRACER`, `OBCS_APPLY_UVICE`, `OBCS_SEAICE_SPONGE_UVICE`, `OBCS_ADJUST_UVICE` (+1)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/obsfit.md`
- L1 pkg/obsfit — `useOBSFIT`
  - L11 Namelist parameters
    - L12 OBSFIT_NML — `obsfitDir`, `obsfitFiles`, `mult_obsfit`, `obsfit_facmod`, `obsfitDoNcOutput`, `obsfitDoGenGrid`
  - L20 CPP options (defaults as shipped) — `OBSFIT_USE_MDSFINDUNITS`, `ALLOW_OBSFIT_EXCLUDE_CORNERS`
  - L24 Headers — `OBSFIT.h`, `OBSFIT_OPTIONS.h`, `OBSFIT_SIZE.h`
  - L29 Routines (29) — `obsfit_active_file.F`, `obsfit_active_file_ad.F`, `obsfit_active_file_control.F`, `obsfit_active_file_g.F`, `obsfit_cost.F`, `obsfit_cost_final.F`, `obsfit_findunit.F`, `obsfit_ini_io.F`, `obsfit_init_equifiles.F`, `obsfit_init_fixed.F`, `obsfit_init_varia.F`, `obsfit_inloop.F`, `obsfit_nc_utils.F`, `obsfit_read_obs.F`, `obsfit_readparms.F`, `obsfit_sampling.F`
  - L32 Called from outside the package — `OBSFIT_INI_IO`, `OBSFIT_INIT_FIXED`, `OBSFIT_INIT_VARIA`, `OBSFIT_READPARMS`, `OBSFIT_COST`, `OBSFIT_INLOOP`, `OBSFIT_NC_CLOSE`, `OBSFIT_COST_FINAL`
  - L42 Verification experiments compiling it (1) — `global_oce_biogeo_bling`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/ocn_compon_interf.md`
- L1 pkg/ocn_compon_interf — `useOCN_COMPON_INTERF`
  - L8 Namelist parameters
    - L9 CPL_OCN_PARAM — `cpl_earlyExpImpCall`, `useImportHFlx`, `useImportFW`, `useImportTau`, `useImportSLP`, `useImportRunOff`, `useImportSIce`, `useImportThSIce`, `useImportSltPlm`, `useImportFice`, `useImportCO2`, `useImportWSpd`, `cpl_taveFreq`, `cpl_snapshot_mnc`, `cpl_timeave_mnc`
  - L26 Headers — `CPL_PARAMS.h`, `OCNCPL.h`, `OCN_CPL_OPTIONS.h`
  - L31 Routines (21) — `cpl_diagnostics_fill.F`, `cpl_diagnostics_init.F`, `cpl_exch_configs.F`, `cpl_export_import_data.F`, `cpl_import_cplparms.F`, `cpl_ini_vars.F`, `cpl_init.F`, `cpl_init_fixed.F`, `cpl_readparms.F`, `cpl_register.F`, `cpl_write_pickup.F`, `ocn_apply_import.F`, `ocn_check_cplconfig.F`, `ocn_cpl_diags.F`, `ocn_cpl_read_pickup.F`, `ocn_export_data.F`, `ocn_export_fields.F`, `ocn_export_ocnconfig.F`, `ocn_import_atmconfig.F`, `ocn_import_fields.F`, `ocn_store_my_data.F`
  - L34 Called from outside the package — `CPL_REGISTER`, `CPL_INIT`, `OCN_APPLY_IMPORT`, `OCN_EXPORT_DATA`, `CPL_EXPORT_IMPORT_DATA`, `CPL_EXCH_CONFIGS`, `CPL_IMPORT_CPLPARMS`, `CPL_INIT_FIXED`, `CPL_INI_VARS`, `CPL_READPARMS`, `CPL_WRITE_PICKUP`, `CPL_DIAGNOSTICS_FILL`, `CPL_DIAGNOSTICS_INIT`
  - L49 Verification experiments compiling it (1)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/offline.md`
- L1 pkg/offline — `useOFFLINE`
  - L10 Namelist parameters
    - L11 OFFLINE_PARM01 — `UvelFile`, `VvelFile`, `WvelFile`, `ThetFile`, `SaltFile`, `GMwxFile`, `GMwyFile`, `GMwzFile`, `ConvFile`, `KPP_DiffSFile`, `KPP_ghatKFile`, `HFluxFile`, `SFluxFile`, `IceFile`, `KPP_ghatFile`
    - L27 OFFLINE_PARM02 — `offlineIter0`, `deltaToffline`, `offlineTimeOffset`, `offlineForcingPeriod`, `offlineForcingCycle`, `offlineLoadPrec`, `offlineOffsetIter`
  - L36 CPP options (defaults as shipped) — `NOT_MODEL_FILES`
  - L39 Headers — `OFFLINE.h`, `OFFLINE_OPTIONS.h`, `OFFLINE_SWITCH.h`, `offline_ad_check_lev1_dir.h`, `offline_ad_check_lev2_dir.h`, `offline_ad_check_lev3_dir.h`, `offline_ad_check_lev4_dir.h`
  - L48 Routines (6) — `offline_check.F`, `offline_fields_load.F`, `offline_get_diffus.F`, `offline_init_varia.F`, `offline_readparms.F`, `offline_reset_parms.F`
  - L51 Called from outside the package — `OFFLINE_GET_DIFFUS`, `OFFLINE_FIELDS_LOAD`, `OFFLINE_CHECK`, `OFFLINE_INIT_VARIA`, `OFFLINE_READPARMS`, `OFFLINE_RESET_PARMS`
  - L60 Verification experiments compiling it (2) — `tutorial_cfc_offline`, `tutorial_dic_adjoffline`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/openad.md`
- L1 pkg/openad — `useOPENAD`
  - L8 CPP options (defaults as shipped) — `ALLOW_OPENAD_ACTIVE_READ_XY`, `ALLOW_OPENAD_ACTIVE_READ_XYZ`, `ALLOW_OPENAD_ACTIVE_WRITE`, `ALLOW_OPENAD_DIVA`, `OAD_DEBUG`
  - L15 Headers — `OPENAD_OPTIONS.h`
  - L18 Routines (11) — `ad_s_different_multiple.F`, `ad_s_ifnblnk.F`, `ad_s_ilnblnk.F`, `ad_s_master_cpu_thread.F`, `ad_s_mds_reclen.F`, `externalDummies.F`, `inner_do_loop.F`, `openad_checkpointInit.F`, `openad_dumpAdjoint.F`
  - L21 Called from outside the package — `EXCH1_RL`, `GLOBAL_SUM_TILE_RL`, `INNER_DO_LOOP`, `STREAMICE_SMOOTH_ADJOINT_FIELD`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/opps.md`
- L1 pkg/opps — `useOPPS`
  - L9 Namelist parameters
    - L10 OPPS_PARM01 — `MAX_ABE_ITERATIONS`, `OPPSdebugLevel`, `PlumeRadius`, `STABILITY_THRESHOLD`, `FRACTIONAL_AREA`, `MAX_FRACTIONAL_AREA`, `VERTICAL_VELOCITY`, `ENTRAINMENT_RATE`, `useGCMwVel`
  - L21 CPP options (defaults as shipped) — `ALLOW_OPPS_SNAPSHOT`, `ALLOW_OPPS_DEBUG`
  - L25 Headers — `OPPS.h`, `OPPS_OPTIONS.h`
  - L29 Routines (7) — `opps_calc.F`, `opps_check.F`, `opps_init.F`, `opps_interface.F`, `opps_readparms.F`
  - L32 Called from outside the package — `OPPS_CHECK`, `OPPS_INIT`, `OPPS_READPARMS`, `OPPS_INTERFACE`
  - L38 Verification experiments compiling it (1)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/overflood@wad.md`
- L1 pkg/overflood  (wad: ~/Documents/research/ECCO/wetting_drying/MITgcm) — `useOVERFLOOD`
  - L10 Namelist parameters
    - L11 OVERFLOOD_PARM01 — `ovfIceFile`, `ovfHoleFile`, `ovfFreeboard`, `ovfManningN`, `ovfHoleCd`, `ovfMinDepth`, `ovfMaxFroude`, `ovfCFL`, `ovfOceanKeep`, `ovfMonFreq`, `ovfIceCh`, `ovfIceUmin`, `ovfUseSeaice`, `ovfAreaMin`, `ovfHeffMin`, `ovfLeadFrac`
  - L29 Headers — `OVERFLOOD.h`, `OVERFLOOD_OPTIONS.h`
  - L33 Routines (12) — `overflood_check.F`, `overflood_diagnostics_init.F`, `overflood_init_fixed.F`, `overflood_init_varia.F`, `overflood_read_pickup.F`, `overflood_readparms.F`, `overflood_step.F`, `overflood_tendency_apply.F`, `overflood_write_pickup.F`
  - L36 Called from outside the package — `OVERFLOOD_TENDENCY_APPLY_S`, `OVERFLOOD_TENDENCY_APPLY_T`, `OVERFLOOD_STEP`, `OVERFLOOD_CHECK`, `OVERFLOOD_INIT_FIXED`, `OVERFLOOD_INIT_VARIA`, `OVERFLOOD_READPARMS`, `OVERFLOOD_WRITE_PICKUP`
  - L46 Verification experiments compiling it (1)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/pp81.md`
- L1 pkg/pp81 — `usePP81`
  - L9 Namelist parameters
    - L10 PP81_PARM01 — `PPnRi`, `PPviscMin`, `PPdiffMin`, `PPviscMax`, `PPdumpFreq`, `PPMixingMaps`, `PPwriteState`
  - L21 CPP options (defaults as shipped) — `ALLOW_PP81_LOWERBOUND`, `PP81_SMOOTH_RI`
  - L25 Headers — `PP81.h`, `PP81_OPTIONS.h`
  - L29 Routines (10) — `pp81_calc.F`, `pp81_calc_diff.F`, `pp81_calc_visc.F`, `pp81_check.F`, `pp81_diagnostics_init.F`, `pp81_init_fixed.F`, `pp81_init_varia.F`, `pp81_output.F`, `pp81_readparms.F`, `pp81_ri_number.F`
  - L32 Called from outside the package — `PP81_CALC_DIFF`, `PP81_CALC_VISC`, `PP81_CALC`, `PP81_OUTPUT`, `PP81_CHECK`, `PP81_INIT_FIXED`, `PP81_INIT_VARIA`, `PP81_READPARMS`
  - L42 Verification experiments compiling it (1)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/profiles.md`
- L1 pkg/profiles — `usePROFILES`
  - L11 Namelist parameters
    - L12 PROFILES_NML — `profilesDir`, `mult_profiles`, `mult_profiles_mean`, `profiles_mean_indsamples`, `prof_facmod`, `prof_names`, `prof_namesmod`, `prof_namesclim`, `prof_itracer`, `profilesDoNcOutput`, `profilesDoGenGrid`, `prof_make_nc`, `prof_dBugLevel`
  - L28 CPP options (defaults as shipped) — `ALLOW_PROFILES_CLIMMASK`, `ALLOW_PROFILES_EXCLUDE_CORNERS`, `ALLOW_PROFILES_SAMPLESPLIT_COST`
  - L33 Headers — `PROFILES_OPTIONS.h`, `PROFILES_SIZE.h`, `profiles.h`
  - L38 Routines (25) — `active_file_control_profiles.F`, `active_file_profiles.F`, `active_file_profiles_ad.F`, `active_file_profiles_g.F`, `active_file_profiles_tap_adj.F`, `cost_profiles.F`, `profiles_cost.F`, `profiles_cost_final.F`, `profiles_findunit.F`, `profiles_ini_io.F`, `profiles_init_fixed.F`, `profiles_init_ncfile.F`, `profiles_init_varia.F`, `profiles_inloop.F`, `profiles_interp.F`, `profiles_make_ncfile.F`, `profiles_nc_utils.F`, `profiles_readparms.F`, `profiles_readvector.F`
  - L41 Called from outside the package — `PROFILES_INI_IO`, `PROFILES_INIT_FIXED`, `PROFILES_INIT_VARIA`, `PROFILES_READPARMS`, `PROFILES_COST`, `PROFILES_INLOOP`, `PROFILES_NC_CLOSE`, `PROFILES_COST_FINAL`
  - L51 Verification experiments compiling it (2) — `global_oce_biogeo_bling`, `global_oce_latlon`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/ptracers.md`
- L1 pkg/ptracers — `usePTRACERS`
  - L11 Namelist parameters
    - L12 PTRACERS_PARM01 — `tauTr1ClimRelax`, `PTRACERS_numInUse`, `PTRACERS_Iter0`, `PTRACERS_doAB_onGpTr`, `PTRACERS_addSrelax2EmP`, `PTRACERS_startStepFwd`, `PTRACERS_resetFreq`, `PTRACERS_resetPhase`, `PTRACERS_advScheme`, `PTRACERS_ImplVertAdv`, `PTRACERS_diffKh`, `PTRACERS_diffK4`, `PTRACERS_diffKr`, `PTRACERS_diffKrNr`, `PTRACERS_ref`, `PTRACERS_EvPrRn`, `PTRACERS_useGMRedi`, `PTRACERS_useDWNSLP`, `PTRACERS_useKPP`, `PTRACERS_linFSConserve`, `PTRACERS_stayPositive`, `PTRACERS_initialFile`, `PTRACERS_names`, `PTRACERS_long_names`, `PTRACERS_units`, `PTRACERS_useRecords`, `PTRACERS_dumpFreq`, `PTRACERS_taveFreq`, `PTRACERS_monitorFreq`, `PTRACERS_timeave_mnc` (+4)
  - L48 CPP options (defaults as shipped) — `PTRACERS_ALLOW_DYN_STATE`
  - L51 Headers — `PTRACERS_FIELDS.h`, `PTRACERS_MOD.h`, `PTRACERS_OPTIONS.h`, `PTRACERS_PARAMS.h`, `PTRACERS_SIZE.h`, `PTRACERS_START.h`, `ptracers_ad_check_lev1_dir.h`, `ptracers_ad_check_lev2_dir.h`, `ptracers_ad_check_lev3_dir.h`, `ptracers_ad_check_lev4_dir.h`, `ptracers_adcommon.h`
  - L64 Routines (32) — `ptracers_ad_dump.F`, `ptracers_apply_forcing.F`, `ptracers_calc_wsurf_tr.F`, `ptracers_check.F`, `ptracers_check_pickup.F`, `ptracers_convect.F`, `ptracers_debug.F`, `ptracers_diagnostics_init.F`, `ptracers_diagnostics_state.F`, `ptracers_dyn_state_data_mod.F`, `ptracers_dyn_state_mod.F`, `ptracers_fields_blocking_exch.F`, `ptracers_forcing_surf.F`, `ptracers_init_fixed.F`, `ptracers_init_varia.F`, `ptracers_integrate.F`, `ptracers_mnc_init.F`, `ptracers_monitor.F`, `ptracers_monitor_ad.F`, `ptracers_output.F`, `ptracers_read_pickup.F`, `ptracers_readparms.F`, `ptracers_reset.F`, `ptracers_set_iolabel.F`, `ptracers_switch_onoff.F`, `ptracers_timeave.F`, `ptracers_turnoff_io.F`, `ptracers_write_pickup.F`, `ptracers_write_state.F`, `ptracers_write_timeave.F` (+1)
  - L67 Called from outside the package — `PTRACERS_CONVECT`, `PTRACERS_FIELDS_BLOCKING_EXCH`, `PTRACERS_DIAGNOSTICS_STATE`, `PTRACERS_OUTPUT`, `PTRACERS_FORCING_SURF`, `PTRACERS_RESET`, `PTRACERS_SWITCH_ONOFF`, `PTRACERS_CHECK`, `PTRACERS_INIT_FIXED`, `PTRACERS_INIT_VARIA`, `PTRACERS_READPARMS`, `PTRACERS_WRITE_PICKUP`, `PTRACERS_CALC_WSURF_TR`, `PTRACERS_DEBUG`, `PTRACERS_INTEGRATE`, `PTRACERS_ZONAL_FILT_APPLY`, `PTRACERS_TURNOFF_IO`, `PTRACERS_AD_DUMP`, `ADPTRACERS_MONITOR`
  - L92 Verification experiments compiling it (12) — `cfc_example`, `global_oce_biogeo_bling`, `global_ocean.90x40x15`, `lab_sea`, `so_box_biogeo`, `tutorial_advection_in_gyre`, `tutorial_cfc_offline`, `tutorial_dic_adjoffline`, `tutorial_global_oce_biogeo`, `tutorial_global_oce_latlon`, `tutorial_tracer_adjsens`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/ptracers@bbl.md`
- L1 pkg/ptracers  (bbl: ~/Documents/research/ECCO/BBL/MITgcm) — `usePTRACERS`
  - L13 Namelist parameters
    - L14 PTRACERS_PARM01 — `tauTr1ClimRelax`, `PTRACERS_numInUse`, `PTRACERS_Iter0`, `PTRACERS_doAB_onGpTr`, `PTRACERS_addSrelax2EmP`, `PTRACERS_startStepFwd`, `PTRACERS_resetFreq`, `PTRACERS_resetPhase`, `PTRACERS_advScheme`, `PTRACERS_ImplVertAdv`, `PTRACERS_diffKh`, `PTRACERS_diffK4`, `PTRACERS_diffKr`, `PTRACERS_diffKrNr`, `PTRACERS_ref`, `PTRACERS_EvPrRn`, `PTRACERS_useGMRedi`, `PTRACERS_useDWNSLP`, `PTRACERS_useKPP`, `PTRACERS_linFSConserve`, `PTRACERS_stayPositive`, `PTRACERS_initialFile`, `PTRACERS_names`, `PTRACERS_long_names`, `PTRACERS_units`, `PTRACERS_useRecords`, `PTRACERS_dumpFreq`, `PTRACERS_taveFreq`, `PTRACERS_monitorFreq`, `PTRACERS_timeave_mnc` (+4)
  - L50 CPP options (defaults as shipped) — `PTRACERS_ALLOW_DYN_STATE`
  - L53 Headers — `PTRACERS_FIELDS.h`, `PTRACERS_MOD.h`, `PTRACERS_OPTIONS.h`, `PTRACERS_PARAMS.h`, `PTRACERS_SIZE.h`, `PTRACERS_START.h`, `ptracers_ad_check_lev1_dir.h`, `ptracers_ad_check_lev2_dir.h`, `ptracers_ad_check_lev3_dir.h`, `ptracers_ad_check_lev4_dir.h`, `ptracers_adcommon.h`
  - L66 Routines (30) — `ptracers_ad_dump.F`, `ptracers_apply_forcing.F`, `ptracers_calc_wsurf_tr.F`, `ptracers_check.F`, `ptracers_check_pickup.F`, `ptracers_convect.F`, `ptracers_debug.F`, `ptracers_diagnostics_init.F`, `ptracers_diagnostics_state.F`, `ptracers_dyn_state_data_mod.F`, `ptracers_dyn_state_mod.F`, `ptracers_fields_blocking_exch.F`, `ptracers_forcing_surf.F`, `ptracers_init_fixed.F`, `ptracers_init_varia.F`, `ptracers_integrate.F`, `ptracers_mnc_init.F`, `ptracers_monitor.F`, `ptracers_monitor_ad.F`, `ptracers_output.F`, `ptracers_read_pickup.F`, `ptracers_readparms.F`, `ptracers_reset.F`, `ptracers_set_iolabel.F`, `ptracers_switch_onoff.F`, `ptracers_turnoff_io.F`, `ptracers_write_pickup.F`, `ptracers_write_state.F`, `ptracers_zonal_filt_apply.F`
  - L69 Called from outside the package — `PTRACERS_CONVECT`, `PTRACERS_FIELDS_BLOCKING_EXCH`, `PTRACERS_DIAGNOSTICS_STATE`, `PTRACERS_OUTPUT`, `PTRACERS_FORCING_SURF`, `PTRACERS_RESET`, `PTRACERS_SWITCH_ONOFF`, `PTRACERS_CHECK`, `PTRACERS_INIT_FIXED`, `PTRACERS_INIT_VARIA`, `PTRACERS_READPARMS`, `PTRACERS_WRITE_PICKUP`, `PTRACERS_CALC_WSURF_TR`, `PTRACERS_DEBUG`, `PTRACERS_INTEGRATE`, `PTRACERS_ZONAL_FILT_APPLY`, `PTRACERS_TURNOFF_IO`, `PTRACERS_AD_DUMP`, `ADPTRACERS_MONITOR`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/ptracers@bbl_c68g.md`
- L1 pkg/ptracers  (bbl_c68g: ~/Documents/research/ECCO/BBL/MITgcm_c68g) — `usePTRACERS`
  - L13 Namelist parameters
    - L14 PTRACERS_PARM01 — `tauTr1ClimRelax`, `PTRACERS_numInUse`, `PTRACERS_Iter0`, `PTRACERS_doAB_onGpTr`, `PTRACERS_addSrelax2EmP`, `PTRACERS_startStepFwd`, `PTRACERS_resetFreq`, `PTRACERS_resetPhase`, `PTRACERS_advScheme`, `PTRACERS_ImplVertAdv`, `PTRACERS_diffKh`, `PTRACERS_diffK4`, `PTRACERS_diffKr`, `PTRACERS_diffKrNr`, `PTRACERS_ref`, `PTRACERS_EvPrRn`, `PTRACERS_useGMRedi`, `PTRACERS_useDWNSLP`, `PTRACERS_useKPP`, `PTRACERS_linFSConserve`, `PTRACERS_stayPositive`, `PTRACERS_initialFile`, `PTRACERS_names`, `PTRACERS_long_names`, `PTRACERS_units`, `PTRACERS_useRecords`, `PTRACERS_dumpFreq`, `PTRACERS_taveFreq`, `PTRACERS_monitorFreq`, `PTRACERS_timeave_mnc` (+4)
  - L50 CPP options (defaults as shipped) — `PTRACERS_ALLOW_DYN_STATE`
  - L53 Headers — `PTRACERS_FIELDS.h`, `PTRACERS_MOD.h`, `PTRACERS_OPTIONS.h`, `PTRACERS_PARAMS.h`, `PTRACERS_SIZE.h`, `PTRACERS_START.h`, `PTRACERS_TAVE.h`, `ptracers_ad_check_lev1_dir.h`, `ptracers_ad_check_lev2_dir.h`, `ptracers_ad_check_lev3_dir.h`, `ptracers_ad_check_lev4_dir.h`, `ptracers_adcommon.h`
  - L67 Routines (32) — `ptracers_ad_dump.F`, `ptracers_apply_forcing.F`, `ptracers_calc_wsurf_tr.F`, `ptracers_check.F`, `ptracers_check_pickup.F`, `ptracers_convect.F`, `ptracers_debug.F`, `ptracers_diagnostics_init.F`, `ptracers_diagnostics_state.F`, `ptracers_dyn_state_data_mod.F`, `ptracers_dyn_state_mod.F`, `ptracers_fields_blocking_exch.F`, `ptracers_forcing_surf.F`, `ptracers_init_fixed.F`, `ptracers_init_varia.F`, `ptracers_integrate.F`, `ptracers_mnc_init.F`, `ptracers_monitor.F`, `ptracers_monitor_ad.F`, `ptracers_output.F`, `ptracers_read_pickup.F`, `ptracers_readparms.F`, `ptracers_reset.F`, `ptracers_set_iolabel.F`, `ptracers_switch_onoff.F`, `ptracers_timeave.F`, `ptracers_turnoff_io.F`, `ptracers_write_pickup.F`, `ptracers_write_state.F`, `ptracers_write_timeave.F` (+1)
  - L70 Called from outside the package — `PTRACERS_CONVECT`, `PTRACERS_FIELDS_BLOCKING_EXCH`, `PTRACERS_DIAGNOSTICS_STATE`, `PTRACERS_TIMEAVE`, `PTRACERS_OUTPUT`, `PTRACERS_FORCING_SURF`, `PTRACERS_RESET`, `PTRACERS_SWITCH_ONOFF`, `PTRACERS_CHECK`, `PTRACERS_INIT_FIXED`, `PTRACERS_INIT_VARIA`, `PTRACERS_READPARMS`, `PTRACERS_WRITE_PICKUP`, `PTRACERS_CALC_WSURF_TR`, `PTRACERS_DEBUG`, `PTRACERS_INTEGRATE`, `PTRACERS_ZONAL_FILT_APPLY`, `PTRACERS_TURNOFF_IO`, `PTRACERS_AD_DUMP`, `ADPTRACERS_MONITOR`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/ptracers@wad.md`
- L1 pkg/ptracers  (wad: ~/Documents/research/ECCO/wetting_drying/MITgcm) — `usePTRACERS`
  - L13 Namelist parameters
    - L14 PTRACERS_PARM01 — `tauTr1ClimRelax`, `PTRACERS_numInUse`, `PTRACERS_Iter0`, `PTRACERS_doAB_onGpTr`, `PTRACERS_addSrelax2EmP`, `PTRACERS_startStepFwd`, `PTRACERS_resetFreq`, `PTRACERS_resetPhase`, `PTRACERS_advScheme`, `PTRACERS_ImplVertAdv`, `PTRACERS_diffKh`, `PTRACERS_diffK4`, `PTRACERS_diffKr`, `PTRACERS_diffKrNr`, `PTRACERS_ref`, `PTRACERS_EvPrRn`, `PTRACERS_useGMRedi`, `PTRACERS_useDWNSLP`, `PTRACERS_useKPP`, `PTRACERS_linFSConserve`, `PTRACERS_stayPositive`, `PTRACERS_initialFile`, `PTRACERS_names`, `PTRACERS_long_names`, `PTRACERS_units`, `PTRACERS_useRecords`, `PTRACERS_dumpFreq`, `PTRACERS_taveFreq`, `PTRACERS_monitorFreq`, `PTRACERS_timeave_mnc` (+4)
  - L50 CPP options (defaults as shipped) — `PTRACERS_ALLOW_DYN_STATE`
  - L53 Headers — `PTRACERS_FIELDS.h`, `PTRACERS_MOD.h`, `PTRACERS_OPTIONS.h`, `PTRACERS_PARAMS.h`, `PTRACERS_SIZE.h`, `PTRACERS_START.h`, `ptracers_ad_check_lev1_dir.h`, `ptracers_ad_check_lev2_dir.h`, `ptracers_ad_check_lev3_dir.h`, `ptracers_ad_check_lev4_dir.h`, `ptracers_adcommon.h`
  - L66 Routines (30) — `ptracers_ad_dump.F`, `ptracers_apply_forcing.F`, `ptracers_calc_wsurf_tr.F`, `ptracers_check.F`, `ptracers_check_pickup.F`, `ptracers_convect.F`, `ptracers_debug.F`, `ptracers_diagnostics_init.F`, `ptracers_diagnostics_state.F`, `ptracers_dyn_state_data_mod.F`, `ptracers_dyn_state_mod.F`, `ptracers_fields_blocking_exch.F`, `ptracers_forcing_surf.F`, `ptracers_init_fixed.F`, `ptracers_init_varia.F`, `ptracers_integrate.F`, `ptracers_mnc_init.F`, `ptracers_monitor.F`, `ptracers_monitor_ad.F`, `ptracers_output.F`, `ptracers_read_pickup.F`, `ptracers_readparms.F`, `ptracers_reset.F`, `ptracers_set_iolabel.F`, `ptracers_switch_onoff.F`, `ptracers_turnoff_io.F`, `ptracers_write_pickup.F`, `ptracers_write_state.F`, `ptracers_zonal_filt_apply.F`
  - L69 Called from outside the package — `PTRACERS_CONVECT`, `PTRACERS_FIELDS_BLOCKING_EXCH`, `PTRACERS_DIAGNOSTICS_STATE`, `PTRACERS_OUTPUT`, `PTRACERS_FORCING_SURF`, `PTRACERS_RESET`, `PTRACERS_SWITCH_ONOFF`, `PTRACERS_CHECK`, `PTRACERS_INIT_FIXED`, `PTRACERS_INIT_VARIA`, `PTRACERS_READPARMS`, `PTRACERS_WRITE_PICKUP`, `PTRACERS_CALC_WSURF_TR`, `PTRACERS_DEBUG`, `PTRACERS_INTEGRATE`, `PTRACERS_ZONAL_FILT_APPLY`, `PTRACERS_TURNOFF_IO`, `PTRACERS_AD_DUMP`, `ADPTRACERS_MONITOR`
  - L94 Verification experiments compiling it (3)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/radtrans.md`
- L1 pkg/radtrans  (darwin3: ~/Documents/GitHub/darwin3) — `useRADTRANS`
  - L8 Namelist parameters
    - L9 RADTRANS_FORCING_PARAMS — `RT_Edfile`, `RT_Esfile`, `RT_E_mask`, `RT_E_period`, `RT_E_RepCycle`, `RT_E_startTime`, `RT_E_startdate1`, `RT_E_startdate2`, `RT_Ed_const`, `RT_Ed_exfremo_intercept`, `RT_Ed_exfremo_slope`, `RT_inscal_Ed`, `RT_Es_const`, `RT_Es_exfremo_intercept`, `RT_Es_exfremo_slope`, `RT_inscal_Es`, `RT_icefile`, `RT_iceperiod`, `RT_iceRepCycle`, `RT_iceStartTime`, `RT_icestartdate1`, `RT_icestartdate2`, `RT_iceconst`, `RT_ice_exfremo_intercept`, `RT_ice_exfremo_slope`, `RT_icemask`, `RT_inscal_ice`
    - L37 RADTRANS_INTERP_PARAMS — `RT_E_lon0`, `RT_E_lat0`, `RT_E_nlon`, `RT_E_nlat`, `RT_E_lon_inc`, `RT_E_interpMethod`, `RT_E_lat_inc`, `RT_ice_lon0`, `RT_ice_lat0`, `RT_ice_nlon`, `RT_ice_nlat`, `RT_ice_lon_inc`, `RT_ice_interpMethod`, `RT_ice_lat_inc`
    - L52 RADTRANS_PARAMS — `RT_refract_water`, `RT_rmud_max`, `RT_wbEdges`, `RT_wbRefWLs`, `RT_kmax`, `RT_useOASIMrmud`, `RT_useMeanCosSolz`, `RT_useNoonSolz`, `RT_sfcIrrThresh`, `RT_oasimWgt`
    - L63 RADTRANS_DEPENDENT — `RT_wbWidths`
  - L66 CPP options (defaults as shipped) — `RADTRANS_DIAG_SOLUTION`
  - L69 Headers — `RADTRANS_EXF_PARAMS.h`, `RADTRANS_FIELDS.h`, `RADTRANS_OPTIONS.h`, `RADTRANS_PARAMS.h`, `RADTRANS_SIZE.h`
  - L76 Routines (13) — `radtrans_calc.F`, `radtrans_check.F`, `radtrans_declination_spencer.F`, `radtrans_diagnostics_init.F`, `radtrans_fields_load.F`, `radtrans_init_fixed.F`, `radtrans_init_varia.F`, `radtrans_monitor.F`, `radtrans_readparms.F`, `radtrans_rmud_below.F`, `radtrans_solve.F`, `radtrans_solve_tridiag.F`, `radtrans_solz_daytime.F`
  - L79 Called from outside the package — `RADTRANS_FIELDS_LOAD`, `RADTRANS_CHECK`, `RADTRANS_INIT_FIXED`, `RADTRANS_INIT_VARIA`, `RADTRANS_READPARMS`, `RADTRANS_CALC`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/radtrans@d3backport.md`
- L1 pkg/radtrans  (d3backport: ~/Documents/research/debug/darwin3) — `useRADTRANS`
  - L10 Namelist parameters
    - L11 RADTRANS_FORCING_PARAMS — `RT_Edfile`, `RT_Esfile`, `RT_E_mask`, `RT_E_period`, `RT_E_RepCycle`, `RT_E_startTime`, `RT_E_startdate1`, `RT_E_startdate2`, `RT_Ed_const`, `RT_Ed_exfremo_intercept`, `RT_Ed_exfremo_slope`, `RT_inscal_Ed`, `RT_Es_const`, `RT_Es_exfremo_intercept`, `RT_Es_exfremo_slope`, `RT_inscal_Es`, `RT_icefile`, `RT_iceperiod`, `RT_iceRepCycle`, `RT_iceStartTime`, `RT_icestartdate1`, `RT_icestartdate2`, `RT_iceconst`, `RT_ice_exfremo_intercept`, `RT_ice_exfremo_slope`, `RT_icemask`, `RT_inscal_ice`
    - L39 RADTRANS_INTERP_PARAMS — `RT_E_lon0`, `RT_E_lat0`, `RT_E_nlon`, `RT_E_nlat`, `RT_E_lon_inc`, `RT_E_interpMethod`, `RT_E_lat_inc`, `RT_ice_lon0`, `RT_ice_lat0`, `RT_ice_nlon`, `RT_ice_nlat`, `RT_ice_lon_inc`, `RT_ice_interpMethod`, `RT_ice_lat_inc`
    - L54 RADTRANS_PARAMS — `RT_refract_water`, `RT_rmud_max`, `RT_wbEdges`, `RT_wbRefWLs`, `RT_kmax`, `RT_useOASIMrmud`, `RT_useMeanCosSolz`, `RT_useNoonSolz`, `RT_sfcIrrThresh`, `RT_oasimWgt`
    - L65 RADTRANS_DEPENDENT — `RT_wbWidths`
  - L68 CPP options (defaults as shipped) — `RADTRANS_DIAG_SOLUTION`
  - L71 Headers — `RADTRANS_EXF_PARAMS.h`, `RADTRANS_FIELDS.h`, `RADTRANS_OPTIONS.h`, `RADTRANS_PARAMS.h`, `RADTRANS_SIZE.h`
  - L78 Routines (13) — `radtrans_calc.F`, `radtrans_check.F`, `radtrans_declination_spencer.F`, `radtrans_diagnostics_init.F`, `radtrans_fields_load.F`, `radtrans_init_fixed.F`, `radtrans_init_varia.F`, `radtrans_monitor.F`, `radtrans_readparms.F`, `radtrans_rmud_below.F`, `radtrans_solve.F`, `radtrans_solve_tridiag.F`, `radtrans_solz_daytime.F`
  - L81 Called from outside the package — `RADTRANS_FIELDS_LOAD`, `RADTRANS_CHECK`, `RADTRANS_INIT_FIXED`, `RADTRANS_INIT_VARIA`, `RADTRANS_READPARMS`, `RADTRANS_CALC`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/rbcs.md`
- L1 pkg/rbcs — `useRBCS`
  - L10 Namelist parameters
    - L11 RBCS_PARM01 — `tauRelaxU`, `tauRelaxV`, `tauRelaxT`, `tauRelaxS`, `relaxMaskUFile`, `relaxMaskVFile`, `relaxMaskFile`, `relaxUFile`, `relaxVFile`, `relaxTFile`, `relaxSFile`, `useRBCuVel`, `useRBCvVel`, `useRBCtemp`, `useRBCsalt`, `useRBCptracers`, `rbcsIniter`, `rbcsForcingPeriod`, `rbcsForcingCycle`, `rbcsForcingOffset`, `rbcsVanishingTime`, `rbcsSingleTimeFiles`, `deltaTrbcs`, `rbcsIter0`
    - L36 RBCS_PARM02 — `useRBCpTrNum`, `tauRelaxPTR`, `relaxPtracerFile`
  - L41 CPP options (defaults as shipped) — `DISABLE_RBCS_MOM`
  - L44 Headers — `RBCS_FIELDS.h`, `RBCS_OPTIONS.h`, `RBCS_PARAMS.h`, `RBCS_SIZE.h`, `rbcs_ad_check_lev1_dir.h`, `rbcs_ad_check_lev2_dir.h`, `rbcs_ad_check_lev3_dir.h`, `rbcs_ad_check_lev4_dir.h`
  - L54 Routines (5) — `rbcs_add_tendency.F`, `rbcs_fields_load.F`, `rbcs_init_fixed.F`, `rbcs_init_varia.F`, `rbcs_readparms.F`
  - L57 Called from outside the package — `RBCS_ADD_TENDENCY`, `RBCS_FIELDS_LOAD`, `RBCS_INIT_FIXED`, `RBCS_INIT_VARIA`, `RBCS_READPARMS`
  - L66 Verification experiments compiling it (2) — `tutorial_reentrant_channel`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/regrid.md`
- L1 pkg/regrid — `useREGRID`
  - L9 Namelist parameters
    - L10 REGRID_PARM01 — `regrid_MNC`, `regrid_MDSIO`, `regrid_ngrids`, `regrid_fbname_in`, `regrid_nout`
  - L17 CPP options (defaults as shipped) — `RL_IS_REAL4`
  - L20 Headers — `REGRID.h`, `REGRID_OPTIONS.h`, `REGRID_SIZE.h`
  - L25 Routines (7) — `regrid_check.F`, `regrid_init_fixed.F`, `regrid_init_varia.F`, `regrid_mnc_init.F`, `regrid_readparms.F`, `regrid_scalar_out_RL.F`, `regrid_scalar_out_RS.F`
  - L28 Called from outside the package — `REGRID_CHECK`, `REGRID_INIT_FIXED`, `REGRID_INIT_VARIA`, `REGRID_READPARMS`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/regrid@bbl.md`
- L1 pkg/regrid  (bbl: ~/Documents/research/ECCO/BBL/MITgcm) — `useREGRID`
  - L11 Namelist parameters
    - L12 REGRID_PARM01 — `regrid_MNC`, `regrid_MDSIO`, `regrid_ngrids`, `regrid_fbname_in`, `regrid_nout`
  - L19 CPP options (defaults as shipped) — `RL_IS_REAL4`
  - L22 Headers — `REGRID.h`, `REGRID_OPTIONS.h`, `REGRID_SIZE.h`
  - L27 Routines (7) — `regrid_check.F`, `regrid_init_fixed.F`, `regrid_init_varia.F`, `regrid_mnc_init.F`, `regrid_readparms.F`, `regrid_scalar_out_RL.F`, `regrid_scalar_out_RS.F`
  - L30 Called from outside the package — `REGRID_CHECK`, `REGRID_INIT_FIXED`, `REGRID_INIT_VARIA`, `REGRID_READPARMS`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/runclock.md`
- L1 pkg/runclock — `useRUNCLOCK`
  - L8 Namelist parameters
    - L9 RUNCLOCK — `RC_maxtime_hr`, `RC_maxtime_mi`, `RC_maxtime_sc`
  - L14 CPP options (defaults as shipped) — `RUNCLOCK_USES_DATE_AND_TIME`
  - L17 Headers — `RUNCLOCK.h`, `RUNCLOCK_OPTIONS.h`
  - L21 Routines (5) — `runclock_check.F`, `runclock_continue.F`, `runclock_gettime.F`, `runclock_init.F`, `runclock_readparms.F`
  - L24 Called from outside the package — `RUNCLOCK_CHECK`, `RUNCLOCK_INIT`, `RUNCLOCK_READPARMS`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/rw.md`
- L1 pkg/rw
  - L10 CPP options (defaults as shipped) — `RW_SAFE_MFLDS`, `RW_DISABLE_SMALL_OVERLAP`
  - L14 Headers — `RW_MFLDS.h`, `RW_OPTIONS.h`
  - L18 Routines (54) — `get_write_global_fld.F`, `read_fld_xy_rl.F`, `read_fld_xy_rs.F`, `read_fld_xyz_rl.F`, `read_fld_xyz_rs.F`, `read_glvec_rl.F`, `read_glvec_rs.F`, `read_mflds.F`, `read_rec.F`, `rw_get_suffix.F`, `set_write_global_fld.F`, `write_fld_3d_rl.F`, `write_fld_3d_rs.F`, `write_fld_s3d_rl.F`, `write_fld_s3d_rs.F`, `write_fld_xy_rl.F`, `write_fld_xy_rs.F`, `write_fld_xyz_rl.F`, `write_fld_xyz_rs.F`, `write_glvec_rl.F`, `write_glvec_rs.F`, `write_local_rl.F`, `write_local_rs.F`, `write_rec.F`
  - L21 Called from outside the package — `READ_FLD_XY_RS`, `WRITE_FLD_S3D_RL`, `RW_GET_SUFFIX`, `READ_REC_XY_RS`, `WRITE_FLD_XYZ_RS`, `READ_REC_3D_RS`, `READ_FLD_XYZ_RL`, `WRITE_FLD_XY_RL`, `WRITE_FLD_XY_RS`, `READ_MFLDS_INIT`, `SET_WRITE_GLOBAL_FLD`, `SET_WRITE_GLOBAL_REC`, `SET_WRITE_GLOBAL_SEC`, `READ_FLD_XY_RL`, `READ_GLVEC_RL`, `READ_GLVEC_RS`, `READ_MFLDS_3D_RL`, `READ_MFLDS_CHECK`, `READ_MFLDS_LEV_RS`, `READ_MFLDS_SET`, `READ_REC_3D_RL`, `WRITE_GLVEC_RL`, `WRITE_FLD_XYZ_RL`, `READ_FLD_XYZ_RS`, `WRITE_GLVEC_RS`, `WRITE_REC_3D_RL`, `WRITE_REC_3D_RS`, `WRITE_REC_XYZ_RL`, `WRITE_REC_XY_RL`, `WRITE_LOCAL_RL` (+15)
  - L440 Verification experiments compiling it (58) — `advect_cs`, `advect_xz`, `aim.5l_Equatorial_Channel`, `aim.5l_LatLon`, `aim.5l_cs`, `atm_gray`, `bottom_ctrl_5x5`, `cfc_example`, `cheapAML_box`, `deep_anelastic`, `front_relax`, `global_oce_biogeo_bling`, `global_oce_latlon`, `global_ocean.90x40x15`, `global_ocean.cs32x15`, `halfpipe_streamice`, `ideal_2D_oce`, `internal_wave`, `inverted_barometer`, `lab_sea`, `matrix_example`, `obcs_ctrl`, `offline_exf_seaice`, `seaice_itd`, `seaice_obcs`, `shelfice_2d_remesh`, `short_surf_wave`, `so_box_biogeo`, `tutorial_advection_in_gyre`, `tutorial_baroclinic_gyre` (+12)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/salt_plume.md`
- L1 pkg/salt_plume — `useSALT_PLUME`, `data.salt_plume`
  - L10 Namelist parameters
    - L11 SALT_PLUME_PARM01 — `SaltPlumeSouthernOcean`, `CriterionType`, `PlumeMethod`, `SaltPlumeCriterion`, `SPsalFRAC`, `SPinflectionPoint`, `SaltPlumeSplitBasin`, `SPbrineSconst`, `SPbrineSaltmax`
  - L24 CPP options (defaults as shipped) — `SALT_PLUME_IN_LEADS`, `SALT_PLUME_SPLIT_BASIN`, `SALT_PLUME_VOLUME`
  - L29 Headers — `SALT_PLUME.h`, `SALT_PLUME_OPTIONS.h`, `salt_plume_ad_check_lev1_dir.h`, `salt_plume_ad_check_lev2_dir.h`, `salt_plume_ad_check_lev3_dir.h`, `salt_plume_ad_check_lev4_dir.h`
  - L37 Routines (15) — `salt_plume_apply.F`, `salt_plume_calc_depth.F`, `salt_plume_check.F`, `salt_plume_diagnostics_fill.F`, `salt_plume_diagnostics_init.F`, `salt_plume_do_exch.F`, `salt_plume_forcing_surf.F`, `salt_plume_frac.F`, `salt_plume_init_fixed.F`, `salt_plume_init_varia.F`, `salt_plume_mnc_init.F`, `salt_plume_readparms.F`, `salt_plume_tendency_apply_s.F`, `salt_plume_tendency_apply_t.F`, `salt_plume_volfrac.F`
  - L40 Called from outside the package — `SALT_PLUME_TENDENCY_APPLY_S`, `SALT_PLUME_TENDENCY_APPLY_T`, `SALT_PLUME_APPLY`, `SALT_PLUME_CALC_DEPTH`, `SALT_PLUME_DIAGNOSTICS_FILL`, `SALT_PLUME_DO_EXCH`, `SALT_PLUME_FORCING_SURF`, `SALT_PLUME_VOLFRAC`, `SALT_PLUME_CHECK`, `SALT_PLUME_INIT_FIXED`, `SALT_PLUME_INIT_VARIA`, `SALT_PLUME_READPARMS`, `SALT_PLUME_FRAC`
  - L59 Verification experiments compiling it (3) — `lab_sea`, `seaice_obcs`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/sbo.md`
- L1 pkg/sbo — `useSBO`
  - L9 Namelist parameters
    - L10 SBO_PARM01 — `sbo_taveFreq`, `sbo_monFreq`
  - L14 Headers — `SBO.h`, `SBO_OPTIONS.h`
  - L18 Routines (5) — `sbo_calc.F`, `sbo_check.F`, `sbo_output.F`, `sbo_readparms.F`, `sbo_rho.F`
  - L21 Called from outside the package — `SBO_CALC`, `SBO_OUTPUT`, `SBO_CHECK`, `SBO_READPARMS`
  - L27 Verification experiments compiling it (1) — `global_ocean.90x40x15`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/seaice.md`
- L1 pkg/seaice — `useSEAICE`
  - L10 Namelist parameters
    - L11 SEAICE_PARM01 — `SEAICEuseDYNAMICS`, `SEAICEuseFREEDRIFT`, `SEAICEuseStrImpCpl`, `SEAICEuseMCS`, `SEAICEuseMCE`, `SEAICEuseTD`, `SEAICEusePL`, `SEAICEuseTEM`, `SEAICEuseTilt`, `SEAICEuseMetricTerms`, `SEAICEselectMetricTerms`, `useHB87stressCoupling`, `SEAICEupdateOceanStress`, `usePW79thermodynamics`, `useMaykutSatVapPoly`, `SEAICEuseFlooding`, `SEAICErestoreUnderIce`, `SEAICE_growMeltByConv`, `SEAICE_salinityTracer`, `SEAICE_ageTracer`, `SEAICEadvHeff`, `SEAICEadvArea`, `SEAICEadvSnow`, `SEAICEadvSalt`, `SEAICEadvAge`, `SEAICEaddSnowMass`, `SEAICEmomAdvection`, `SEAICEselectKEscheme`, `SEAICEselectVortScheme`, `SEAICEhighOrderVorticity` (+208)
    - L258 SEAICE_PARM02 — `mult_ice_export`, `mult_ice`, `cost_ice_flag`, `costIceStart1`, `costIceStart2`, `costIceEnd1`, `costIceEnd2`, `SEAICE_cutoff_area`, `SEAICE_cutoff_heff`, `SEAICE_clamp_salt`, `SEAICE_clamp_theta`, `mult_smrsst`, `mult_smrsss`, `mult_smrarea`, `wmean_smrarea`, `smrarea_errfile`
    - L283 SEAICE_PARM03 — `SItrFile`, `SItrName`, `SItrNameLong`, `SItrUnit`, `SItrMate`, `SItrFromOcean0`, `SItrFromOceanFrac`, `SItrFromFlood0`, `SItrFromFloodFrac`, `SItrExpand0`, `IceAgeTrFile`, `SItrNumInUse`
  - L297 CPP options (defaults as shipped) — `SEAICE_DEBUG`, `SEAICE_EXTERNAL_FLUXES`, `SEAICE_ITD`, `SEAICE_DISABLE_SUBLIM`, `SEAICE_DISABLE_HEATCONSFIX`, `SEAICE_VARIABLE_SALINITY`, `SEAICE_GREASE`, `ALLOW_SITRACER`, `ALLOW_SITRACER_ADVCAP`, `ALLOW_SITRACER_DEBUG_DIAG`, `SEAICE_BGRID_DYNAMICS`, `SEAICE_CGRID`, `SEAICE_ALLOW_MOM_ADVECTION`, `SEAICE_ALLOW_BOTTOMDRAG`, `SEAICE_ALLOW_SIDEDRAG`, `SEAICE_ALLOW_JFNK`, `SEAICE_ALLOW_KRYLOV`, `SEAICE_JFNK_MAP_REORDER`, `SEAICE_PRECOND_EXTRA_EXCHANGE`, `SEAICE_GLOBAL_3DIAG_SOLVER`, `SEAICE_ALLOW_EVP`, `SEAICE_ALLOW_CLIPZETA`, `SEAICE_EVP_ELIMINATE_UNDERFLOWS`, `ALLOW_SEAICE_EVP_RESIDUAL`, `SEAICE_DELTA_SMOOTHREG`, `SEAICE_ZETA_SMOOTHREG`, `SEAICE_ALLOW_TEM`, `SEAICE_ALLOW_MCS`, `SEAICE_ALLOW_MCE`, `SEAICE_ALLOW_TEARDROP` (+22)
  - L352 Headers — `SEAICE.h`, `SEAICE_COST.h`, `SEAICE_GRID.h`, `SEAICE_OPTIONS.h`, `SEAICE_PARAMS.h`, `SEAICE_SIZE.h`, `SEAICE_TRACER.h`, `seaice_ad_check_lev1_dir.h`, `seaice_ad_check_lev2_dir.h`, `seaice_ad_check_lev3_dir.h`, `seaice_ad_check_lev4_dir.h`
  - L365 Routines (92) — `advect.F`, `cost_ice_test.F`, `diffus.F`, `dynsolver.F`, `lsr.F`, `ostres.F`, `seaice_ad_dump.F`, `seaice_advdiff.F`, `seaice_advection.F`, `seaice_bottomdrag_coeffs.F`, `seaice_budget_ocean.F`, `seaice_calc_ice_strength.F`, `seaice_calc_lhs.F`, `seaice_calc_residual.F`, `seaice_calc_rhs.F`, `seaice_calc_strainrates.F`, `seaice_calc_stressdiv.F`, `seaice_calc_viscosities.F`, `seaice_check.F`, `seaice_check_pickup.F`, `seaice_cost_accumulate_mean.F`, `seaice_cost_areasst.F`, `seaice_cost_concentration.F`, `seaice_cost_driver.F`, `seaice_cost_export.F`, `seaice_cost_final.F`, `seaice_cost_init_fixed.F`, `seaice_cost_init_varia.F`, `seaice_cost_sensi.F`, `seaice_cost_sss.F` (+47)
  - L368 Called from outside the package — `SEAICE_COST_SENSI`, `SEAICE_FAKE`, `SEAICE_MODEL`, `SEAICE_DIAGNOSTICS_STATE`, `SEAICE_OUTPUT`, `SEAICE_CHECK`, `SEAICE_COST_INIT_FIXED`, `SEAICE_INIT_FIXED`, `SEAICE_INIT_VARIA`, `SEAICE_READPARMS`, `SEAICE_WRITE_PICKUP`, `SEAICE_TURNOFF_IO`, `SEAICE_AD_DUMP`, `SEAICE_COST_FINAL`, `SEAICE_COST_INIT_VARIA`, `ADSEAICE_MONITOR`
  - L386 Verification experiments compiling it (7) — `global_ocean.cs32x15`, `lab_sea`, `offline_exf_seaice`, `seaice_itd`, `seaice_obcs`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/seaice@seaicebc.md`
- L1 pkg/seaice  (seaicebc: ~/Documents/research/ECCO/sea_ice_BCs/MITgcm) — `useSEAICE`
  - L12 Namelist parameters
    - L13 SEAICE_PARM01 — `SEAICEuseDYNAMICS`, `SEAICEuseFREEDRIFT`, `SEAICEuseStrImpCpl`, `SEAICEuseMCS`, `SEAICEuseMCE`, `SEAICEuseTD`, `SEAICEusePL`, `SEAICEuseTEM`, `SEAICEuseTilt`, `SEAICEuseMetricTerms`, `SEAICEselectMetricTerms`, `useHB87stressCoupling`, `SEAICEupdateOceanStress`, `usePW79thermodynamics`, `useMaykutSatVapPoly`, `SEAICEuseFlooding`, `SEAICErestoreUnderIce`, `SEAICE_growMeltByConv`, `SEAICE_salinityTracer`, `SEAICE_ageTracer`, `SEAICEadvHeff`, `SEAICEadvArea`, `SEAICEadvSnow`, `SEAICEadvSalt`, `SEAICEadvAge`, `SEAICEaddSnowMass`, `SEAICEmomAdvection`, `SEAICEselectKEscheme`, `SEAICEselectVortScheme`, `SEAICEhighOrderVorticity` (+208)
    - L260 SEAICE_PARM02 — `mult_ice_export`, `mult_ice`, `cost_ice_flag`, `costIceStart1`, `costIceStart2`, `costIceEnd1`, `costIceEnd2`, `SEAICE_cutoff_area`, `SEAICE_cutoff_heff`, `SEAICE_clamp_salt`, `SEAICE_clamp_theta`, `mult_smrsst`, `mult_smrsss`, `mult_smrarea`, `wmean_smrarea`, `smrarea_errfile`
    - L285 SEAICE_PARM03 — `SItrFile`, `SItrName`, `SItrNameLong`, `SItrUnit`, `SItrMate`, `SItrFromOcean0`, `SItrFromOceanFrac`, `SItrFromFlood0`, `SItrFromFloodFrac`, `SItrExpand0`, `IceAgeTrFile`, `SItrNumInUse`
  - L299 CPP options (defaults as shipped) — `SEAICE_DEBUG`, `SEAICE_EXTERNAL_FLUXES`, `SEAICE_ITD`, `SEAICE_DISABLE_SUBLIM`, `SEAICE_DISABLE_HEATCONSFIX`, `SEAICE_VARIABLE_SALINITY`, `SEAICE_GREASE`, `ALLOW_SITRACER`, `ALLOW_SITRACER_ADVCAP`, `ALLOW_SITRACER_DEBUG_DIAG`, `SEAICE_BGRID_DYNAMICS`, `SEAICE_CGRID`, `SEAICE_ALLOW_MOM_ADVECTION`, `SEAICE_ALLOW_BOTTOMDRAG`, `SEAICE_ALLOW_SIDEDRAG`, `SEAICE_ALLOW_JFNK`, `SEAICE_ALLOW_KRYLOV`, `SEAICE_JFNK_MAP_REORDER`, `SEAICE_PRECOND_EXTRA_EXCHANGE`, `SEAICE_GLOBAL_3DIAG_SOLVER`, `SEAICE_ALLOW_EVP`, `SEAICE_ALLOW_CLIPZETA`, `SEAICE_EVP_ELIMINATE_UNDERFLOWS`, `ALLOW_SEAICE_EVP_RESIDUAL`, `SEAICE_DELTA_SMOOTHREG`, `SEAICE_ZETA_SMOOTHREG`, `SEAICE_ALLOW_TEM`, `SEAICE_ALLOW_MCS`, `SEAICE_ALLOW_MCE`, `SEAICE_ALLOW_TEARDROP` (+22)
  - L354 Headers — `SEAICE.h`, `SEAICE_COST.h`, `SEAICE_GRID.h`, `SEAICE_OPTIONS.h`, `SEAICE_PARAMS.h`, `SEAICE_SIZE.h`, `SEAICE_TRACER.h`, `seaice_ad_check_lev1_dir.h`, `seaice_ad_check_lev2_dir.h`, `seaice_ad_check_lev3_dir.h`, `seaice_ad_check_lev4_dir.h`
  - L367 Routines (84) — `advect.F`, `diffus.F`, `dynsolver.F`, `lsr.F`, `ostres.F`, `seaice_ad_dump.F`, `seaice_advdiff.F`, `seaice_advection.F`, `seaice_bottomdrag_coeffs.F`, `seaice_budget_ocean.F`, `seaice_calc_ice_strength.F`, `seaice_calc_lhs.F`, `seaice_calc_residual.F`, `seaice_calc_rhs.F`, `seaice_calc_strainrates.F`, `seaice_calc_stressdiv.F`, `seaice_calc_viscosities.F`, `seaice_check.F`, `seaice_check_pickup.F`, `seaice_cost_accumulate_mean.F`, `seaice_cost_export.F`, `seaice_cost_final.F`, `seaice_cost_init_fixed.F`, `seaice_cost_init_varia.F`, `seaice_cost_sensi.F`, `seaice_cost_test.F`, `seaice_diagnostics_init.F`, `seaice_diagnostics_state.F`, `seaice_diffusion.F`, `seaice_do_ridging.F` (+39)
  - L370 Called from outside the package — `SEAICE_COST_SENSI`, `SEAICE_FAKE`, `SEAICE_MODEL`, `SEAICE_DIAGNOSTICS_STATE`, `SEAICE_OUTPUT`, `SEAICE_CHECK`, `SEAICE_COST_INIT_FIXED`, `SEAICE_INIT_FIXED`, `SEAICE_INIT_VARIA`, `SEAICE_READPARMS`, `SEAICE_WRITE_PICKUP`, `SEAICE_TURNOFF_IO`, `SEAICE_AD_DUMP`, `SEAICE_COST_FINAL`, `SEAICE_COST_INIT_VARIA`, `ADSEAICE_MONITOR`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/seaice@wad.md`
- L1 pkg/seaice  (wad: ~/Documents/research/ECCO/wetting_drying/MITgcm) — `useSEAICE`
  - L12 Namelist parameters
    - L13 SEAICE_PARM01 — `SEAICEuseDYNAMICS`, `SEAICEuseFREEDRIFT`, `SEAICEuseStrImpCpl`, `SEAICEuseMCS`, `SEAICEuseMCE`, `SEAICEuseTD`, `SEAICEusePL`, `SEAICEuseTEM`, `SEAICEuseTilt`, `SEAICEuseMetricTerms`, `SEAICEselectMetricTerms`, `useHB87stressCoupling`, `SEAICEupdateOceanStress`, `SEAICEimmobileDrag`, `usePW79thermodynamics`, `useMaykutSatVapPoly`, `SEAICEuseFlooding`, `SEAICErestoreUnderIce`, `SEAICE_growMeltByConv`, `SEAICE_salinityTracer`, `SEAICE_ageTracer`, `SEAICEadvHeff`, `SEAICEadvArea`, `SEAICEadvSnow`, `SEAICEadvSalt`, `SEAICEadvAge`, `SEAICEaddSnowMass`, `SEAICEmomAdvection`, `SEAICEselectKEscheme`, `SEAICEselectVortScheme` (+209)
    - L261 SEAICE_PARM02 — `mult_ice_export`, `mult_ice`, `cost_ice_flag`, `costIceStart1`, `costIceStart2`, `costIceEnd1`, `costIceEnd2`, `SEAICE_cutoff_area`, `SEAICE_cutoff_heff`, `SEAICE_clamp_salt`, `SEAICE_clamp_theta`, `mult_smrsst`, `mult_smrsss`, `mult_smrarea`, `wmean_smrarea`, `smrarea_errfile`
    - L286 SEAICE_PARM03 — `SItrFile`, `SItrName`, `SItrNameLong`, `SItrUnit`, `SItrMate`, `SItrFromOcean0`, `SItrFromOceanFrac`, `SItrFromFlood0`, `SItrFromFloodFrac`, `SItrExpand0`, `IceAgeTrFile`, `SItrNumInUse`
  - L300 CPP options (defaults as shipped) — `SEAICE_DEBUG`, `SEAICE_EXTERNAL_FLUXES`, `SEAICE_ITD`, `SEAICE_DISABLE_SUBLIM`, `SEAICE_DISABLE_HEATCONSFIX`, `SEAICE_VARIABLE_SALINITY`, `SEAICE_GREASE`, `ALLOW_SITRACER`, `ALLOW_SITRACER_ADVCAP`, `ALLOW_SITRACER_DEBUG_DIAG`, `SEAICE_BGRID_DYNAMICS`, `SEAICE_CGRID`, `SEAICE_ALLOW_MOM_ADVECTION`, `SEAICE_ALLOW_BOTTOMDRAG`, `SEAICE_ALLOW_SIDEDRAG`, `SEAICE_ALLOW_JFNK`, `SEAICE_ALLOW_KRYLOV`, `SEAICE_JFNK_MAP_REORDER`, `SEAICE_PRECOND_EXTRA_EXCHANGE`, `SEAICE_GLOBAL_3DIAG_SOLVER`, `SEAICE_ALLOW_EVP`, `SEAICE_ALLOW_CLIPZETA`, `SEAICE_EVP_ELIMINATE_UNDERFLOWS`, `ALLOW_SEAICE_EVP_RESIDUAL`, `SEAICE_DELTA_SMOOTHREG`, `SEAICE_ZETA_SMOOTHREG`, `SEAICE_ALLOW_TEM`, `SEAICE_ALLOW_MCS`, `SEAICE_ALLOW_MCE`, `SEAICE_ALLOW_TEARDROP` (+22)
  - L355 Headers — `SEAICE.h`, `SEAICE_COST.h`, `SEAICE_GRID.h`, `SEAICE_OPTIONS.h`, `SEAICE_PARAMS.h`, `SEAICE_SIZE.h`, `SEAICE_TRACER.h`, `seaice_ad_check_lev1_dir.h`, `seaice_ad_check_lev2_dir.h`, `seaice_ad_check_lev3_dir.h`, `seaice_ad_check_lev4_dir.h`
  - L368 Routines (84) — `advect.F`, `diffus.F`, `dynsolver.F`, `lsr.F`, `ostres.F`, `seaice_ad_dump.F`, `seaice_advdiff.F`, `seaice_advection.F`, `seaice_bottomdrag_coeffs.F`, `seaice_budget_ocean.F`, `seaice_calc_ice_strength.F`, `seaice_calc_lhs.F`, `seaice_calc_residual.F`, `seaice_calc_rhs.F`, `seaice_calc_strainrates.F`, `seaice_calc_stressdiv.F`, `seaice_calc_viscosities.F`, `seaice_check.F`, `seaice_check_pickup.F`, `seaice_cost_accumulate_mean.F`, `seaice_cost_export.F`, `seaice_cost_final.F`, `seaice_cost_init_fixed.F`, `seaice_cost_init_varia.F`, `seaice_cost_sensi.F`, `seaice_cost_test.F`, `seaice_diagnostics_init.F`, `seaice_diagnostics_state.F`, `seaice_diffusion.F`, `seaice_do_ridging.F` (+39)
  - L371 Called from outside the package — `SEAICE_COST_SENSI`, `SEAICE_FAKE`, `SEAICE_MODEL`, `SEAICE_DIAGNOSTICS_STATE`, `SEAICE_OUTPUT`, `SEAICE_CHECK`, `SEAICE_COST_INIT_FIXED`, `SEAICE_INIT_FIXED`, `SEAICE_INIT_VARIA`, `SEAICE_READPARMS`, `SEAICE_WRITE_PICKUP`, `SEAICE_TURNOFF_IO`, `SEAICE_AD_DUMP`, `SEAICE_COST_FINAL`, `SEAICE_COST_INIT_VARIA`, `ADSEAICE_MONITOR`
  - L389 Verification experiments compiling it (4)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/sediment@wad.md`
- L1 pkg/sediment  (wad: ~/Documents/research/ECCO/wetting_drying/MITgcm) — `useSEDIMENT`
  - L10 Namelist parameters
    - L11 SEDIMENT_PARM01 — `sedNcl`, `sedPtrIdx`, `sedWs`, `sedTauCE`, `sedEroM`, `sedTauCD`, `sedRiverC`, `sedRiverCFile`, `sedBulkRho`, `sedBedFile`, `sedDryDepth`, `sedDryDeposit`, `sedRefDepth`, `sedMonFreq`, `sedConsTime`, `sedTauCEf`, `sedEroMf`
    - L29 SEDIMENT_PARM02 — `sedNlay`, `sedCohesive`, `sedPorosity`, `sedFcLo`, `sedFcHi`, `sedTcbNew`, `sedTcbOffset`, `sedTcbSlope`, `sedTcbMin`, `sedTcbMax`, `sedTcons`, `sedTswell`, `sedZaCoef`, `sedZaMin`, `sedLayMass`
    - L45 SEDIMENT_PARM03 — `sedWaveMode`, `sedWaveUbFile`, `sedWavePerFile`, `sedWaveDirFile`, `sedWaveHsFile`, `sedWaveNrec`, `sedWaveDt`, `sedWindSpeed`, `sedWindDir`, `sedFetch`, `sedWaveGamma`, `sedZ0`, `sedWaveIce`
  - L60 Headers — `SEDIMENT.h`, `SEDIMENT_OPTIONS.h`
  - L64 Routines (13) — `sediment_add_tendency.F`, `sediment_bed.F`, `sediment_check.F`, `sediment_diagnostics_init.F`, `sediment_init_fixed.F`, `sediment_init_varia.F`, `sediment_read_pickup.F`, `sediment_readparms.F`, `sediment_step.F`, `sediment_waves.F`, `sediment_write_pickup.F`
  - L67 Called from outside the package — `SEDIMENT_STEP`, `SEDIMENT_CHECK`, `SEDIMENT_INIT_FIXED`, `SEDIMENT_INIT_VARIA`, `SEDIMENT_READPARMS`, `SEDIMENT_WRITE_PICKUP`, `SEDIMENT_ADD_TENDENCY`
  - L76 Verification experiments compiling it (2)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/shap_filt.md`
- L1 pkg/shap_filt — `useSHAP_FILT`
  - L12 Namelist parameters
    - L13 SHAP_PARM01 — `Shap_funct`, `shap_filt_uvStar`, `shap_filt_TrStagg`, `Shap_alwaysExchUV`, `Shap_alwaysExchTr`, `nShapT`, `nShapS`, `nShapTrPhys`, `Shap_Trtau`, `Shap_TrLength`, `nShapUV`, `nShapUVPhys`, `Shap_uvtau`, `Shap_uvLength`, `Shap_noSlip`, `Shap_diagFreq`
  - L31 CPP options (defaults as shipped) — `SEQUENTIAL_2D_SHAP`, `USE_OLD_SHAPIRO_FILTERS`, `NO_SLIP_SHAP`, `USE_SHAP_CALC_VORTICITY`
  - L37 Headers — `SHAP_FILT.h`, `SHAP_FILT_OPTIONS.h`
  - L41 Routines (17) — `shap_filt_apply_ts.F`, `shap_filt_apply_uv.F`, `shap_filt_computvort.F`, `shap_filt_diagnostics_init.F`, `shap_filt_init_fixed.F`, `shap_filt_readparms.F`, `shap_filt_relvort3.F`, `shap_filt_tracer_s1.F`, `shap_filt_tracer_s2.F`, `shap_filt_tracer_s4.F`, `shap_filt_tracerold.F`, `shap_filt_u.F`, `shap_filt_uv_s1.F`, `shap_filt_uv_s2.F`, `shap_filt_uv_s2c.F`, `shap_filt_uv_s4.F`, `shap_filt_v.F`
  - L44 Called from outside the package — `SHAP_FILT_APPLY_UV`, `SHAP_FILT_INIT_FIXED`, `SHAP_FILT_READPARMS`, `SHAP_FILT_APPLY_TS`
  - L51 Verification experiments compiling it (12) — `aim.5l_Equatorial_Channel`, `aim.5l_LatLon`, `aim.5l_cs`, `atm_gray`, `tutorial_held_suarez_cs`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/shelfice.md`
- L1 pkg/shelfice — `useSHELFICE`
  - L10 Namelist parameters
    - L11 SHELFICE_PARM01 — `SHELFICEsaltToHeatRatio`, `SHELFICEheatTransCoeff`, `SHELFICEsaltTransCoeff`, `SHELFICEMassStepping`, `SHI_update_kTopC`, `rhoShelfice`, `SHELFICElatentHeat`, `SHELFICEHeatCapacity_Cp`, `no_slip_shelfice`, `SHELFICEDragLinear`, `SHELFICEDragQuadratic`, `SHELFICEselectDragQuadr`, `SHELFICEthetaSurface`, `useISOMIPTD`, `SHELFICEboundaryLayer`, `SHI_withBL_realFWflux`, `SHI_withBL_uStarTopDz`, `SHELFICEwriteState`, `SHELFICE_dumpFreq`, `SHELFICE_taveFreq`, `SHELFICE_tave_mnc`, `SHELFICE_dump_mnc`, `SHELFICEtopoFile`, `SHELFICEmassFile`, `SHELFICEloadAnomalyFile`, `SHELFICEMassDynTendFile`, `SHELFICETransCoeffTFile`, `SHELFICEDynMassOnly`, `SHELFICEadvDiffHeatFlux`, `SHELFICEuseGammaFrict` (+13)
  - L60 CPP options (defaults as shipped) — `ALLOW_ISOMIP_TD`, `SHI_ALLOW_GAMMAFRICT`, `SHI_SALTBAL_FWFLX`, `ALLOW_SHELFICE_REMESHING`, `SHELFICE_REMESH_PRINT`
  - L67 Headers — `SHELFICE.h`, `SHELFICE_COST.h`, `SHELFICE_OPTIONS.h`, `shelfice_ad_check_lev1_dir.h`, `shelfice_ad_check_lev2_dir.h`, `shelfice_ad_check_lev3_dir.h`, `shelfice_ad_check_lev4_dir.h`
  - L76 Routines (27) — `shelfice_check.F`, `shelfice_cost_accumulate.F`, `shelfice_cost_final.F`, `shelfice_cost_shifwflx.F`, `shelfice_diagnostics_drag.F`, `shelfice_forcing.F`, `shelfice_forcing_surf.F`, `shelfice_init_depths.F`, `shelfice_init_fixed.F`, `shelfice_init_varia.F`, `shelfice_mask_seaice.F`, `shelfice_mnc_init.F`, `shelfice_output.F`, `shelfice_read_pickup.F`, `shelfice_readparms.F`, `shelfice_remesh_c_mask.F`, `shelfice_remesh_calc_w.F`, `shelfice_remesh_state.F`, `shelfice_remesh_uv_mask.F`, `shelfice_remeshing.F`, `shelfice_step_icemass.F`, `shelfice_thermodynamics.F`, `shelfice_u_drag_coeff.F`, `shelfice_v_drag_coeff.F`, `shelfice_write_pickup.F`
  - L79 Called from outside the package — `SHELFICE_FORCING_S`, `SHELFICE_FORCING_T`, `SHELFICE_DIAGNOSTICS_DRAG`, `SHELFICE_THERMODYNAMICS`, `SHELFICE_OUTPUT`, `SHELFICE_FORCING_SURF`, `SHELFICE_REMESHING`, `SHELFICE_INIT_DEPTHS`, `SHELFICE_CHECK`, `SHELFICE_INIT_FIXED`, `SHELFICE_INIT_VARIA`, `SHELFICE_READPARMS`, `SHELFICE_WRITE_PICKUP`, `SHELFICE_COST_FINAL`, `SHELFICE_COST_ACCUMULATE`, `SHELFICE_U_DRAG_COEFF`, `SHELFICE_V_DRAG_COEFF`, `SHELFICE_NETMASSFLUX_SURF`, `SHELFICE_MASK_SEAICE`
  - L109 Verification experiments compiling it (2) — `shelfice_2d_remesh`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/showflops.md`
- L1 pkg/showflops — `useSHOWFLOPS`
  - L8 CPP options (defaults as shipped) — `USE_FLIPS`
  - L11 Headers — `SHOWFLOPS.h`, `SHOWFLOPS_INIT.h`, `SHOWFLOPS_OPTIONS.h`
  - L16 Routines (3) — `showflops_init.F`, `showflops_inloop.F`, `showflops_insolve.F`
  - L19 Called from outside the package — `SHOWFLOPS_INLOOP`, `SHOWFLOPS_INSOLVE`, `SHOWFLOPS_INIT`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/smooth.md`
- L1 pkg/smooth — `useSMOOTH`
  - L10 Namelist parameters
    - L11 SMOOTH_NML — `smooth2D_Lx0`, `smooth2D_Ly0`, `smooth2DmaskName`, `smooth3DtypeH`, `smooth3DsizeH`, `smooth3DtypeZ`, `smooth3DsizeZ`, `smooth3D_Lx0`, `smooth3D_Ly0`, `smooth3D_Lz0`, `smooth3DmaskName`, `smoothDir`
  - L31 Headers — `SMOOTH.h`, `SMOOTH_OPTIONS.h`
  - L35 Routines (19) — `smooth2d.F`, `smooth3d.F`, `smooth_basic2d.F`, `smooth_check.F`, `smooth_correl2d.F`, `smooth_correl2dw.F`, `smooth_correl3d.F`, `smooth_diff2d.F`, `smooth_diff3d.F`, `smooth_filtervar2d.F`, `smooth_filtervar3d.F`, `smooth_hetero2d.F`, `smooth_impldiff.F`, `smooth_init2d.F`, `smooth_init3d.F`, `smooth_init_fixed.F`, `smooth_init_varia.F`, `smooth_readparms.F`, `smooth_rhs.F`
  - L38 Called from outside the package — `SMOOTH_CHECK`, `SMOOTH_INIT_FIXED`, `SMOOTH_INIT_VARIA`, `SMOOTH_READPARMS`, `SMOOTH_CORREL2D`, `SMOOTH_CORREL3D`, `SMOOTH_BASIC2D`, `SMOOTH_HETERO2D`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/sphere.md`
- L1 pkg/sphere — `useSPHERE`
  - L8 Headers — `sphere.h`
  - L11 Routines (13) — `exch_all_2d_rl.F`, `sphere.F`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/steep_icecavity.md`
- L1 pkg/steep_icecavity — `useSTEEP_ICECAVITY`
  - L9 Namelist parameters
    - L10 STIC_PARM01 — `STIClengthFile`, `STICdepthFile`
  - L14 CPP options (defaults as shipped) — `ALLOW_SHITRANSCOEFF_3D`
  - L17 Headers — `STIC.h`, `STIC_OPTIONS.h`
  - L21 Routines (7) — `stic_check.F`, `stic_init_depths.F`, `stic_init_fixed.F`, `stic_init_varia.F`, `stic_readparms.F`, `stic_solve4fluxes.F`, `stic_thermodynamics.F`
  - L24 Called from outside the package — `STIC_THERMODYNAMICS`, `STIC_INIT_DEPTHS`, `STIC_CHECK`, `STIC_INIT_FIXED`, `STIC_INIT_VARIA`, `STIC_READPARMS`
  - L32 Verification experiments compiling it (1)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/streamice.md`
- L1 pkg/streamice — `useSTREAMICE`
  - L10 Namelist parameters
    - L11 STREAMICE_PARM01 — `streamice_density`, `streamice_density_ocean_avg`, `streamice_density_firn`, `B_glen_isothermal`, `n_glen`, `eps_glen_min`, `eps_u_min`, `C_basal_fric_const`, `n_basal_friction`, `streamice_vel_update`, `streamice_cg_tol`, `streamice_nonlin_tol`, `streamice_nonlin_tol_fp`, `streamice_err_norm`, `streamice_max_cg_iter`, `streamice_max_nl_iter`, `streamice_maxcgiter_cpl`, `streamice_maxnliter_cpl`, `STREAMICEthickInit`, `STREAMICEsigcoordInit`, `STREAMICEsigcoordFile`, `STREAMICEthickFile`, `STREAMICEcalveMaskFile`, `STREAMICEcostMaskFile`, `STREAMICE_dump_mdsio`, `STREAMICE_tave_mdsio`, `STREAMICE_dump_mnc`, `STREAMICE_tave_mnc`, `STREAMICE_move_front`, `STREAMICE_calve_to_mask` (+69)
    - L111 STREAMICE_PARMTRACER — `STREAMICETrac2DBCxFile`, `STREAMICETrac2DBCyFile`
    - L115 STREAMICE_PARMPETSC — `PETSC_PRECOND_TYPE`, `PETSC_SOLVER_TYPE`, `streamice_use_petsc`, `streamice_maxnliter_petsc`, `streamice_petsc_pcfactorlevels`
    - L121 STREAMICE_PARMOAD — `streamice_nonlin_tol_adjoint`, `streamice_nonlin_tol_adjoint_rl`, `STREAMICE_OAD_petsc_reuse`, `PETSC_PRECOND_OAD`
    - L126 STREAMICE_COST — `STREAMICEvelOptimSnapBasename`, `STREAMICEvelOptimTCBasename`, `STREAMICEsurfOptimTCBasename`, `STREAMICEBglenCostMaskFile`, `streamice_wgt_drift`, `streamice_wgt_vel`, `streamice_wgt_vel_norm`, `streamice_wgt_surf`, `streamice_wgt_tikh_beta`, `streamice_wgt_tikh_bglen`, `streamice_wgt_tikh_gen`, `streamice_wgt_prior_bglen`, `streamice_wgt_prior_gen`, `STREAMICE_do_snapshot_cost`, `STREAMICE_do_timedep_cost`, `STREAMICE_do_verification_cost`, `STREAMICE_do_vaf_cost`, `streamice_vel_cost_timesteps`, `streamice_surf_cost_timesteps`
    - L146 STREAMICE_PARM02 — `shelf_max_draft`, `shelf_min_draft`, `shelf_edge_pos`, `shelf_slope_scale`, `shelf_flat_width`, `flow_dir`
    - L153 STREAMICE_PARM03 — `min_x_noflow_NORTH`, `max_x_noflow_NORTH`, `min_x_noflow_SOUTH`, `max_x_noflow_SOUTH`, `min_y_noflow_WEST`, `max_y_noflow_WEST`, `min_y_noflow_EAST`, `max_y_noflow_EAST`, `min_x_noStress_NORTH`, `max_x_noStress_NORTH`, `min_x_noStress_SOUTH`, `max_x_noStress_SOUTH`, `min_y_noStress_WEST`, `max_y_noStress_WEST`, `min_y_noStress_EAST`, `max_y_noStress_EAST`, `min_x_FluxBdry_NORTH`, `max_x_FluxBdry_NORTH`, `min_x_FluxBdry_SOUTH`, `max_x_FluxBdry_SOUTH`, `min_y_FluxBdry_WEST`, `max_y_FluxBdry_WEST`, `min_y_FluxBdry_EAST`, `max_y_FluxBdry_EAST`, `min_x_Dirich_NORTH`, `max_x_Dirich_NORTH`, `min_x_Dirich_SOUTH`, `max_x_Dirich_SOUTH`, `min_y_Dirich_WEST`, `max_y_Dirich_WEST` (+16)
  - L201 CPP options (defaults as shipped) — `STREAMICE_CONSTRUCT_MATRIX`, `STREAMICE_HYBRID_STRESS`, `STREAMICE_FLOWLINE_BUTTRESS`, `USE_ALT_RLOW`, `STREAMICE_GEOM_FILE_SETUP`, `STREAMICE_SMOOTH_FLOATATION`, `STREAMICE_SMOOTH_FLOATATION2`, `ALLOW_PETSC`, `ALLOW_STREAMICE_2DTRACER`, `STREAMICE_TRACER_AB`, `STREAMICE_SERIAL_TRISOLVE`, `STREAMICE_3D_GLEN_CONST`, `STREAMICE_COULOMB_SLIDING`, `STREAMICE_ECSECRYO_DOSUM`, `STREAMICE_FALSE`, `STREAMICE_FIRN_CORRECTION`, `STREAMICE_PETSC_3_8`, `STREAMICE_STRESS_BOUNDARY_CONTROL`, `ALLOW_STREAMICE_TIMEDEP_FORCING`, `ALLOW_STREAMICE_FLUX_CONTROL`, `ALLOW_STREAMICE_TC_COST`, `ALLOW_STREAMICE_FP_ADJ`
  - L225 Headers — `STREAMICE.h`, `STREAMICE_ADV.h`, `STREAMICE_BDRY.h`, `STREAMICE_CG.h`, `STREAMICE_COST_SIZE.h`, `STREAMICE_CTRL_FLUX.h`, `STREAMICE_FP.h`, `STREAMICE_OPTIONS.h`, `STREAMICE_PETSC_MOD.h`, `streamice_ad_check_lev1_dir.h`, `streamice_ad_check_lev2_dir.h`, `streamice_ad_check_lev3_dir.h`, `streamice_ad_check_lev4_dir.h`
  - L240 Routines (66) — `adstreamice_cg_solve.F`, `adstreamice_invert_surf_forthick.F`, `eta_gl_prime_streamice.F`, `eta_gl_streamice.F`, `eta_gl_streamice_prime.F`, `phi_gl_streamice.F`, `slope_limiter.F`, `streamice_adv_flux_fl_x.F`, `streamice_adv_flux_fl_y.F`, `streamice_adv_front.F`, `streamice_advect_2dtracer.F`, `streamice_advect_thickness.F`, `streamice_apply_flux_ctrl.F`, `streamice_bstress_exponent.F`, `streamice_cg_functions.F`, `streamice_cg_solve.F`, `streamice_cg_solve_matfree.F`, `streamice_cg_solve_petsc.F`, `streamice_cg_wrapper.F`, `streamice_check.F`, `streamice_cost_accum.F`, `streamice_cost_final.F`, `streamice_cost_reg_accum.F`, `streamice_cost_surf_accum.F`, `streamice_cost_vel_accum.F`, `streamice_diagnostics_state.F`, `streamice_driving_stress.F`, `streamice_driving_stress_fd.F`, `streamice_driving_stress_ppm.F`, `streamice_dump.F` (+32)
  - L243 Called from outside the package — `STREAMICE_DIAGNOSTICS_STATE`, `STREAMICE_TIMESTEP`, `STREAMICE_CHECK`, `STREAMICE_INIT_FIXED`, `STREAMICE_INIT_VARIA`, `STREAMICE_READPARMS`, `STREAMICE_WRITE_PICKUP`, `STREAMICE_FINALIZE_PETSC`, `STREAMICE_COST_FINAL`, `STREAMICE_COST_ACCUM`, `STREAMICE_INITIALIZE_PETSC`, `ADSTREAMICE_CG_SOLVE`, `STREAMICE_CG_SOLVE`
  - L259 Verification experiments compiling it (1) — `halfpipe_streamice`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/tapenade.md`
- L1 pkg/tapenade — `useTAPENADE`
  - L8 CPP options (defaults as shipped) — `ALLOW_TAPENADE_ACTIVE_READ_XYZ`, `ALLOW_TAPENADE_ACTIVE_READ_XY`, `ALLOW_TAPENADE_ACTIVE_WRITE`
  - L13 Headers — `COST_TAP_TLM.h`, `TAPENADE_OPTIONS.h`, `adBinomial.h`, `adComplex.h`, `adFixedPoint.h`, `adStack.h`
  - L21 Routines (35) — `active_read_tap.F`, `active_write_tap.F`, `dummy_tap.F`, `stubs_tap_adj.F`, `stubs_tap_tlm.F`
  - L24 Called from outside the package — `ADEXCH_3D_RL`, `ADEXCH_UV_3D_RL`, `ADEXCH_UV_XY_RS`, `ADEXCH_XY_RL`, `ADEXCH_XY_RS`
  - L39 Verification experiments compiling it (8) — `global_oce_biogeo_bling`, `global_oce_latlon`, `global_ocean.cs32x15`, `halfpipe_streamice`, `lab_sea`, `tutorial_global_oce_biogeo`, `tutorial_tracer_adjsens`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/thsice.md`
- L1 pkg/thsice — `useTHSICE`
  - L10 Namelist parameters
    - L11 THSICE_CONST — `cpIce`, `cpWater`, `kIce`, `kSnow`, `bMeltCoef`, `albColdSnow`, `albWarmSnow`, `tempSnowAlb`, `albOldSnow`, `hNewSnowAge`, `snowAgTime`, `albIceMax`, `albIceMin`, `hAlbIce`, `hAlbSnow`, `i0swFrac`, `dhSnowLin`, `saltIce`, `S_winton`, `mu_Tf`, `nitMaxTsf`, `hIceMin`, `hiMax`, `hsMax`, `iceMaskMax`, `iceMaskMin`, `fracEnMelt`, `fracEnFreez`, `hThinIce`, `hThickIce` (+1)
    - L52 THSICE_PARM01 — `startIceModel`, `stepFwd_oceMxL`, `thSIce_calc_albNIR`, `thSIce_skipThermo`, `thSIce_deltaT`, `thSIce_dtTemp`, `ocean_deltaT`, `tauRelax_MxL`, `tauRelax_MxL_salt`, `hMxL_default`, `sMxL_default`, `vMxL_default`, `thSIce_diffK`, `thSIceAdvScheme`, `stressReduction`, `thSIceBalanceAtmFW`, `thSIce_taveFreq`, `thSIce_diagFreq`, `thSIce_monFreq`, `thSIce_tave_mnc`, `thSIce_snapshot_mnc`, `thSIce_mon_mnc`, `thSIce_pickup_read_mnc`, `thSIce_pickup_write_mnc`, `thSIceFract_InitFile`, `thSIceThick_InitFile`, `thSIceSnowH_InitFile`, `thSIceSnowA_InitFile`, `thSIceEnthp_InitFile`, `thSIceTsurf_InitFile`
    - L83 THSICE_COST — `mult_thsice`, `thsice_cost_ice_flag`
  - L87 CPP options (defaults as shipped) — `THSICE_FRACEN_POWERLAW`, `ALLOW_DBUG_THSICE`, `CHECK_ENERGY_CONSERV`, `THSICE_REGULARIZE_CALC_THICKN`
  - L93 Headers — `THSICE_COST.h`, `THSICE_DEBUG.h`, `THSICE_OPTIONS.h`, `THSICE_PARAMS.h`, `THSICE_SIZE.h`, `THSICE_VARS.h`, `thsice_ad_check_lev1_dir.h`, `thsice_ad_check_lev2_dir.h`, `thsice_ad_check_lev3_dir.h`, `thsice_ad_check_lev4_dir.h`, `thsice_test_addfluid.h`
  - L106 Routines (41) — `thsice_advdiff.F`, `thsice_advection.F`, `thsice_albedo.F`, `thsice_ave.F`, `thsice_balance_frw.F`, `thsice_calc_thickn.F`, `thsice_check.F`, `thsice_check_conserv.F`, `thsice_cost_driver.F`, `thsice_cost_final.F`, `thsice_cost_init_varia.F`, `thsice_cost_test.F`, `thsice_diagnostics_init.F`, `thsice_diagnostics_state.F`, `thsice_diffusion.F`, `thsice_do_advect.F`, `thsice_do_exch.F`, `thsice_extend.F`, `thsice_get_bulkf.F`, `thsice_get_exf.F`, `thsice_get_ocean.F`, `thsice_get_precip.F`, `thsice_get_velocity.F`, `thsice_impl_temp.F`, `thsice_ini_vars.F`, `thsice_init_fixed.F`, `thsice_main.F`, `thsice_map_exf.F`, `thsice_mnc_init.F`, `thsice_monitor.F` (+11)
  - L109 Called from outside the package — `THSICE_MAIN`, `THSICE_DIAGNOSTICS_STATE`, `THSICE_OUTPUT`, `THSICE_CHECK`, `THSICE_INIT_FIXED`, `THSICE_INI_VARS`, `THSICE_READPARMS`, `THSICE_WRITE_PICKUP`, `THSICE_TURNOFF_IO`, `THSICE_DO_ADVECT`, `THSICE_DO_EXCH`, `THSICE_SLAB_OCEAN`, `THSICE_STEP_FWD`, `THSICE_ALBEDO`, `THSICE_IMPL_TEMP`, `THSICE_AVE`, `THSICE_GET_OCEAN`, `THSICE_COST_FINAL`, `THSICE_COST_INIT_VARIA`, `THSICE_COST_DRIVER`
  - L137 Verification experiments compiling it (4) — `aim.5l_cs`, `global_ocean.cs32x15`, `offline_exf_seaice`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/wad@wad.md`
- L1 pkg/wad  (wad: ~/Documents/research/ECCO/wetting_drying/MITgcm) — `useWAD`
  - L11 Namelist parameters
    - L12 WAD_PARM01 — `wadMinDepth`, `wadCritDepth`, `wadDragDepth`, `wadAdvDepth`, `wadGMDepth`, `wadIceDepth`, `wadLoadDepth`, `wadLoadTau`, `wadIceMinHeff`, `wadKPPDepth`, `wadKPPCapDepth`, `wadSmoothWidth`, `wadMonFreq`, `wadMaxSpeed`, `wadManningN`, `wadDragMax`, `wadMaxFroude`, `wadManningFile`, `wadManningT1`, `wadManningT2`, `wadConserveVol`, `wadCarryVel`, `wadIceGround`, `wadUpwindFace`, `wadDryForcing`, `wadIceTopMelt`, `wadIceSeal`, `wadIceOBpack`
  - L42 CPP options (defaults as shipped) — `WAD_DEBUG`, `WAD_SMOOTH_MASK`
  - L46 Headers — `WAD.h`, `WAD_OPTIONS.h`
  - L50 Routines (17) — `wad_apply_floor.F`, `wad_botdrag.F`, `wad_check.F`, `wad_check_speed.F`, `wad_diagnostics_init.F`, `wad_dry_forcing.F`, `wad_gmredi_taper.F`, `wad_init_fixed.F`, `wad_init_varia.F`, `wad_kpp_taper.F`, `wad_limit_froude.F`, `wad_limit_outflow.F`, `wad_read_pickup.F`, `wad_readparms.F`, `wad_seaice_dry.F`, `wad_update_masks.F`, `wad_write_pickup.F`
  - L53 Called from outside the package — `WAD_DRY_FORCING`, `WAD_APPLY_FLOOR`, `WAD_CHECK_SPEED`, `WAD_LIMIT_OUTFLOW`, `WAD_UPDATE_MASKS`, `WAD_CHECK`, `WAD_INIT_FIXED`, `WAD_INIT_VARIA`, `WAD_READPARMS`, `WAD_WRITE_PICKUP`, `WAD_GMREDI_TAPER`, `WAD_KPP_TAPER`, `WAD_SEAICE_DRY`
  - L69 Verification experiments compiling it (8)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/wad@wadcheckin.md`
- L1 pkg/wad  (wadcheckin: ~/Documents/research/ECCO/MITgcm_wad_checkin) — `useWAD`
  - L11 Namelist parameters
    - L12 WAD_PARM01 — `wadMinDepth`, `wadCritDepth`, `wadDragDepth`, `wadAdvDepth`, `wadGMDepth`, `wadKPPDepth`, `wadKPPCapDepth`, `wadSmoothWidth`, `wadMonFreq`, `wadMaxSpeed`, `wadManningN`, `wadDragMax`, `wadMaxFroude`, `wadManningFile`, `wadManningT1`, `wadManningT2`, `wadConserveVol`, `wadCarryVel`, `wadUpwindFace`, `wadDryForcing`
  - L34 CPP options (defaults as shipped) — `WAD_DEBUG`, `WAD_SMOOTH_MASK`
  - L38 Headers — `WAD.h`, `WAD_OPTIONS.h`
  - L42 Routines (16) — `wad_apply_floor.F`, `wad_botdrag.F`, `wad_check.F`, `wad_check_speed.F`, `wad_diagnostics_init.F`, `wad_dry_forcing.F`, `wad_gmredi_taper.F`, `wad_init_fixed.F`, `wad_init_varia.F`, `wad_kpp_taper.F`, `wad_limit_froude.F`, `wad_limit_outflow.F`, `wad_read_pickup.F`, `wad_readparms.F`, `wad_update_masks.F`, `wad_write_pickup.F`
  - L45 Called from outside the package — `WAD_DRY_FORCING`, `WAD_APPLY_FLOOR`, `WAD_CHECK_SPEED`, `WAD_LIMIT_OUTFLOW`, `WAD_UPDATE_MASKS`, `WAD_CHECK`, `WAD_INIT_FIXED`, `WAD_INIT_VARIA`, `WAD_READPARMS`, `WAD_WRITE_PICKUP`, `WAD_GMREDI_TAPER`, `WAD_KPP_TAPER`
  - L60 Verification experiments compiling it (5)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/pkg/zonal_filt.md`
- L1 pkg/zonal_filt — `useZONAL_FILT`
  - L9 Namelist parameters
    - L10 ZONFILT_PARM01 — `zonal_filt_uvStar`, `zonal_filt_TrStagg`, `zonal_filt_lat`, `zonal_filt_cospow`, `zonal_filt_sinpow`, `zonal_filt_mode2dx`
  - L18 Headers — `FFTPACK.h`, `ZONAL_FILT.h`, `ZONAL_FILT_OPTIONS.h`
  - L23 Routines (34) — `fftpack.F`, `zonal_filt_apply_ts.F`, `zonal_filt_apply_uv.F`, `zonal_filt_init.F`, `zonal_filt_nofill.F`, `zonal_filt_postsmooth.F`, `zonal_filt_presmooth.F`, `zonal_filt_readparms.F`, `zonal_filter.F`
  - L26 Called from outside the package — `ZONAL_FILT_APPLY_UV`, `ZONAL_FILT_INIT`, `ZONAL_FILT_READPARMS`, `ZONAL_FILT_APPLY_TS`, `ZONAL_FILTER`
  - L34 Verification experiments compiling it (2) — `aim.5l_LatLon`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/1D_ocean_ice_column.md`
- L1 verification/1D_ocean_ice_column
  - L3 Build variants (code*/)
  - L14 Input variants (input*/)
  - L22 Reference results — `output_adm.txt`, `output_tlm.txt.gz`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/MLAdjust.md`
- L1 verification/MLAdjust
  - L3 README (first 40 lines)
  - L48 Build variants (code*/)
  - L54 Input variants (input*/)
  - L77 Reference results — `output.A4FlxF.txt`, `output.AhFlxF.txt`, `output.AhStTn.txt`, `output.AhVrDv.txt`, `output.QGLthGM.txt`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/adjustment.128x64x1.md`
- L1 verification/adjustment.128x64x1
  - L3 Build variants (code*/)
  - L7 Input variants (input*/)
  - L12 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/adjustment.cs-32x32x1.md`
- L1 verification/adjustment.cs-32x32x1
  - L3 README (first 40 lines)
  - L48 Build variants (code*/)
  - L60 Input variants (input*/)
  - L70 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/advect_cs.md`
- L1 verification/advect_cs
  - L3 Build variants (code*/)
  - L10 Input variants (input*/)
  - L15 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/advect_xy.md`
- L1 verification/advect_xy
  - L3 Build variants (code*/)
  - L9 Input variants (input*/)
  - L17 Reference results — `output.ab3_c4.txt`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/advect_xz.md`
- L1 verification/advect_xz
  - L3 Build variants (code*/)
  - L9 Input variants (input*/)
  - L20 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/aim.5l_Equatorial_Channel.md`
- L1 verification/aim.5l_Equatorial_Channel
  - L3 README (first 40 lines)
  - L48 Build variants (code*/)
  - L54 Input variants (input*/)
  - L59 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/aim.5l_LatLon.md`
- L1 verification/aim.5l_LatLon
  - L3 README (first 40 lines)
  - L39 Build variants (code*/)
  - L44 Input variants (input*/)
  - L49 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/aim.5l_cs.md`
- L1 verification/aim.5l_cs
  - L3 Build variants (code*/)
  - L10 Input variants (input*/)
  - L18 Reference results — `output.thSI.txt`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/atm_gray.md`
- L1 verification/atm_gray
  - L3 README (first 40 lines) — `atm_phys`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/bottom_ctrl_5x5.md`
- L1 verification/bottom_ctrl_5x5
  - L3 Build variants (code*/)
  - L10 Input variants (input*/)
  - L18 Reference results — `output_adm.facg2d.txt`, `output_adm.txt`, `output_tlm.facg2d.txt.gz`, `output_tlm.txt.gz`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/cfc_example.md`
- L1 verification/cfc_example
  - L3 Build variants (code*/)
  - L10 Input variants (input*/)
  - L15 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/cheapAML_box.md`
- L1 verification/cheapAML_box
  - L3 Build variants (code*/)
  - L9 Input variants (input*/)
  - L17 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/cpl_aim+ocn.md`
- L1 verification/cpl_aim+ocn
  - L3 README (first 40 lines) — `code_cpl`, `code_atm`, `code_ocn`, `input_cpl`, `input_atm`, `input_ocn`
  - L48 Build variants (code*/)
  - L61 Input variants (input*/)
  - L77 Reference results — `atmSTDOUT.0000`, `atmSTDOUT.icedyn`, `ocnSTDOUT.0000`, `ocnSTDOUT.icedyn`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/deep_anelastic.md`
- L1 verification/deep_anelastic
  - L3 Build variants (code*/)
  - L9 Input variants (input*/)
  - L17 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/dome.md`
- L1 verification/dome
  - L3 Build variants (code*/)
  - L10 Input variants (input*/)
  - L15 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/exp2.md`
- L1 verification/exp2
  - L3 README (first 40 lines)
  - L15 Build variants (code*/)
  - L21 Input variants (input*/)
  - L29 Reference results — `output.rigidLid.txt`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/exp4.md`
- L1 verification/exp4
  - L3 README (first 40 lines) — `useOBCSprescribe`, `obcs_calc.F`
  - L48 Build variants (code*/)
  - L55 Input variants (input*/)
  - L69 Reference results — `output.with_flt.txt`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/fizhi-cs-32x32x40.md`
- L1 verification/fizhi-cs-32x32x40
  - L3 Build variants (code*/)
  - L9 Input variants (input*/)
  - L14 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/fizhi-cs-aqualev20.md`
- L1 verification/fizhi-cs-aqualev20
  - L3 Build variants (code*/)
  - L10 Input variants (input*/)
  - L15 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/fizhi-gridalt-hs.md`
- L1 verification/fizhi-gridalt-hs
  - L3 Build variants (code*/)
  - L10 Input variants (input*/)
  - L15 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/front_relax.md`
- L1 verification/front_relax
  - L3 README (first 40 lines) — `GM_background_K`, `rhoNil`
  - L48 Build variants (code*/)
  - L54 Input variants (input*/)
  - L71 Reference results — `output.in_p.txt`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/global_oce_biogeo_bling.md`
- L1 verification/global_oce_biogeo_bling
  - L3 Build variants (code*/)
  - L18 Input variants (input*/)
  - L30 Reference results — `output_adm.obsfit.txt`, `output_adm.txt`, `output_tap_adj.txt`, `output_tlm.obsfit.txt.gz`, `output_tlm.txt.gz`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/global_oce_latlon.md`
- L1 verification/global_oce_latlon
  - L3 README (first 40 lines) — `tutorial_global_oce_latlon`, `global_ocean.90x40x15`, `global_ocean.cs32x15`, `input.in_p`
  - L48 Build variants (code*/)
  - L63 Input variants (input*/)
  - L92 Reference results — `output_adm.ebm.txt`, `output_adm.ggl90.txt`, `output_adm.txt`, `output_adm.w_exf.txt`, `output_tap_adj.txt`, `output_tap_adj.w_exf.txt`, `output_tap_tlm.txt`, `output_tap_tlm.w_exf.txt`, `output_tlm.ggl90.txt.gz`, `output_tlm.txt.gz`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/global_oce_latlon@bbl.md`
- L1 verification/global_oce_latlon@bbl
  - L5 README (first 40 lines) — `tutorial_global_oce_latlon`, `global_ocean.90x40x15`, `global_ocean.cs32x15`, `input.in_p`, `bbl_useNofSpeed`
  - L50 Build variants (code*/)
  - L65 Input variants (input*/)
  - L97 Reference results — `output_adm.ebm.txt`, `output_adm.ggl90.txt`, `output_adm.txt`, `output_adm.w_exf.txt`, `output_tap_adj.txt`, `output_tap_adj.w_exf.txt`, `output_tap_tlm.txt`, `output_tap_tlm.w_exf.txt`, `output_tlm.ggl90.txt.gz`, `output_tlm.txt.gz`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/global_ocean.90x40x15.md`
- L1 verification/global_ocean.90x40x15
  - L3 README (first 40 lines)
  - L48 Build variants (code*/)
  - L59 Input variants (input*/)
  - L80 Reference results — `output_adm.bottomdrag.txt`, `output_adm.kapgm.txt`, `output_adm.kapredi.txt`, `output_adm.txt`, `output_tlm.bottomdrag.txt.gz`, `output_tlm.kapgm.txt.gz`, `output_tlm.kapredi.txt.gz`, `output_tlm.txt.gz`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/global_ocean.cs32x15.md`
- L1 verification/global_ocean.cs32x15
  - L3 README (first 40 lines)
  - L44 Build variants (code*/)
  - L62 Input variants (input*/)
  - L96 Reference results — `output.in_p.txt`, `output.viscA4.txt`, `output_adm.seaice.txt`, `output_adm.seaice_dynmix.txt`, `output_adm.thsice.txt`, `output_adm.txt`, `output_tap_adj.txt`, `output_tap_tlm.txt`, `output_tlm.seaice.txt.gz`, `output_tlm.seaice_dynmix.txt.gz`, `output_tlm.thsice.txt.gz`, `output_tlm.txt.gz`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/halfpipe_streamice.md`
- L1 verification/halfpipe_streamice
  - L3 Build variants (code*/)
  - L18 Input variants (input*/)
  - L28 Reference results — `output_adm.txt`, `output_tap_adj.txt.gz`, `output_tap_tlm.txt.gz`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/hs94.128x64x5.md`
- L1 verification/hs94.128x64x5
  - L3 Build variants (code*/)
  - L10 Input variants (input*/)
  - L15 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/hs94.1x64x5.md`
- L1 verification/hs94.1x64x5
  - L3 README (first 40 lines)
  - L48 Build variants (code*/)
  - L59 Input variants (input*/)
  - L67 Reference results — `output_adm.txt`, `output_tlm.txt.gz`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/hs94.cs-32x32x5.md`
- L1 verification/hs94.cs-32x32x5
  - L3 Build variants (code*/)
  - L10 Input variants (input*/)
  - L18 Reference results — `output.impIGW.txt`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/ideal_2D_oce.md`
- L1 verification/ideal_2D_oce
  - L3 README (first 40 lines)
  - L29 Build variants (code*/)
  - L35 Input variants (input*/)
  - L42 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/internal_wave.md`
- L1 verification/internal_wave
  - L3 README (first 40 lines)
  - L34 Build variants (code*/)
  - L41 Input variants (input*/)
  - L49 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/inverted_barometer.md`
- L1 verification/inverted_barometer
  - L3 README (first 40 lines)
  - L20 Build variants (code*/)
  - L26 Input variants (input*/)
  - L31 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/isomip.md`
- L1 verification/isomip
  - L3 Build variants (code*/)
  - L19 Input variants (input*/)
  - L45 Reference results — `output_adm.htd.txt`, `output_adm.stic.txt`, `output_adm.txt`, `output_tap_adj.txt`, `output_tap_tlm.txt`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/lab_sea.md`
- L1 verification/lab_sea
  - L3 README (first 40 lines) — `SEAICE_PARAMS.h`, `output.natl_box.txt`, `output.salt_plume.txt`, `output_adm.noseaice.txt`, `output_adm.noseaicedyn.txt`, `output_adm.txt`, `output_tap_adj.noecco.txt`, `output_tap_adj.txt`, `output_tlm.txt.gz`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/matrix_example.md`
- L1 verification/matrix_example
  - L3 README (first 40 lines)
  - L31 Build variants (code*/)
  - L37 Input variants (input*/)
  - L42 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/obcs_ctrl.md`
- L1 verification/obcs_ctrl
  - L3 Build variants (code*/)
  - L10 Input variants (input*/)
  - L15 Reference results — `output_adm.txt`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/offline_exf_seaice.md`
- L1 verification/offline_exf_seaice
  - L3 README (first 40 lines)
  - L48 Build variants (code*/)
  - L59 Input variants (input*/)
  - L97 Reference results — `output.dyn_ellnnfr.txt`, `output.dyn_jfnk.txt`, `output.dyn_lsr.txt`, `output.dyn_mce.txt`, `output.dyn_paralens.txt`, `output.dyn_teardrop.txt`, `output_adm.obcs.txt`, `output_adm.thsice.txt`, `output_adm.txt`, `output_tlm.obcs.txt.gz`, `output_tlm.thsice.txt.gz`, `output_tlm.txt.gz`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/seaice_itd.md`
- L1 verification/seaice_itd
  - L3 README (first 40 lines)
  - L48 Build variants (code*/)
  - L55 Input variants (input*/)
  - L64 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/seaice_obcs.md`
- L1 verification/seaice_obcs
  - L3 README (first 40 lines) — `mk_tides.m`
  - L29 Build variants (code*/)
  - L35 Input variants (input*/)
  - L47 Reference results — `output.regDenom.txt`, `output.seaiceSponge.txt`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/shelfice_2d_remesh.md`
- L1 verification/shelfice_2d_remesh
  - L3 README (first 40 lines)
  - L40 Build variants (code*/)
  - L46 Input variants (input*/)
  - L51 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/short_surf_wave.md`
- L1 verification/short_surf_wave
  - L3 Build variants (code*/)
  - L9 Input variants (input*/)
  - L14 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/so_box_biogeo.md`
- L1 verification/so_box_biogeo
  - L3 README (first 40 lines)
  - L48 Build variants (code*/)
  - L54 Input variants (input*/)
  - L65 Reference results — `output.caSat0.txt`, `output.caSat3.txt`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/solid-body.cs-32x32x1.md`
- L1 verification/solid-body.cs-32x32x1
  - L3 README (first 40 lines) — `rSphere`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/tutorial_advection_in_gyre.md`
- L1 verification/tutorial_advection_in_gyre
  - L3 Build variants (code*/)
  - L9 Input variants (input*/)
  - L14 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/tutorial_baroclinic_gyre.md`
- L1 verification/tutorial_baroclinic_gyre
  - L3 README (first 40 lines)
  - L35 Build variants (code*/)
  - L41 Input variants (input*/)
  - L46 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/tutorial_barotropic_gyre.md`
- L1 verification/tutorial_barotropic_gyre
  - L3 README (first 40 lines)
  - L37 Build variants (code*/)
  - L41 Input variants (input*/)
  - L46 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/tutorial_cfc_offline.md`
- L1 verification/tutorial_cfc_offline
  - L3 README (first 40 lines)
  - L42 Build variants (code*/)
  - L48 Input variants (input*/)
  - L56 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/tutorial_deep_convection.md`
- L1 verification/tutorial_deep_convection
  - L3 README (first 40 lines)
  - L33 Build variants (code*/)
  - L39 Input variants (input*/)
  - L47 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/tutorial_dic_adjoffline.md`
- L1 verification/tutorial_dic_adjoffline
  - L3 Build variants (code*/)
  - L10 Input variants (input*/)
  - L15 Reference results — `output_adm.txt`, `output_tlm.txt.gz`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/tutorial_global_oce_biogeo.md`
- L1 verification/tutorial_global_oce_biogeo
  - L3 README (first 40 lines)
  - L48 Build variants (code*/)
  - L64 Input variants (input*/)
  - L74 Reference results — `output_adm.txt`, `output_tap_adj.txt`, `output_tap_tlm.txt`, `output_tlm.txt.gz`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/tutorial_global_oce_in_p.md`
- L1 verification/tutorial_global_oce_in_p
  - L3 README (first 40 lines)
  - L48 Build variants (code*/)
  - L54 Input variants (input*/)
  - L59 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/tutorial_global_oce_latlon.md`
- L1 verification/tutorial_global_oce_latlon
  - L3 README (first 40 lines)
  - L30 Build variants (code*/)
  - L36 Input variants (input*/)
  - L41 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/tutorial_global_oce_optim.md`
- L1 verification/tutorial_global_oce_optim
  - L3 README (first 40 lines)
  - L33 Build variants (code*/)
  - L40 Input variants (input*/)
  - L45 Reference results — `output_adm.txt`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/tutorial_held_suarez_cs.md`
- L1 verification/tutorial_held_suarez_cs
  - L3 Build variants (code*/)
  - L10 Input variants (input*/)
  - L15 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/tutorial_plume_on_slope.md`
- L1 verification/tutorial_plume_on_slope
  - L3 README (first 40 lines) — `zRoughBot`, `roughBot`, `output.roughBot.txt`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/tutorial_reentrant_channel.md`
- L1 verification/tutorial_reentrant_channel
  - L3 README (first 40 lines)
  - L34 Build variants (code*/)
  - L41 Input variants (input*/)
  - L46 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/tutorial_rotating_tank.md`
- L1 verification/tutorial_rotating_tank
  - L3 Build variants (code*/)
  - L10 Input variants (input*/)
  - L15 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/tutorial_tracer_adjsens.md`
- L1 verification/tutorial_tracer_adjsens
  - L3 README (first 40 lines)
  - L35 Build variants (code*/)
  - L47 Input variants (input*/)
  - L58 Reference results — `output_adm.som81.txt`, `output_adm.txt`, `output_tap_adj.txt`, `output_tap_tlm.txt`, `output_tlm.som81.txt.gz`, `output_tlm.txt.gz`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/vermix.md`
- L1 verification/vermix
  - L3 Build variants (code*/)
  - L9 Input variants (input*/)
  - L27 Reference results — `output.gglLC.txt`

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/wad_balzano@wad.md`
- L1 verification/wad_balzano@wad
  - L5 Build variants (code*/)
  - L12 Input variants (input*/)
  - L22 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/wad_balzano@wadcheckin.md`
- L1 verification/wad_balzano@wadcheckin
  - L5 README (first 40 lines) — `wadCritDepth`
  - L50 Build variants (code*/)
  - L57 Input variants (input*/)
  - L67 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/wad_estuary_3d@wad.md`
- L1 verification/wad_estuary_3d@wad
  - L5 Build variants (code*/)
  - L12 Input variants (input*/)
  - L25 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/wad_estuary_3d@wadcheckin.md`
- L1 verification/wad_estuary_3d@wadcheckin
  - L5 README (first 40 lines) — `WAD_NX`, `WAD_NY`, `WAD_NR`
  - L50 Build variants (code*/)
  - L57 Input variants (input*/)
  - L70 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/wad_flat_xz@wad.md`
- L1 verification/wad_flat_xz@wad
  - L5 Build variants (code*/)
  - L12 Input variants (input*/)
  - L20 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/wad_flat_xz@wadcheckin.md`
- L1 verification/wad_flat_xz@wadcheckin
  - L5 README (first 40 lines)

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/wad_iceground@wad.md`
- L1 verification/wad_iceground@wad
  - L5 Build variants (code*/)
  - L12 Input variants (input*/)
  - L20 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/wad_mangrove@wad.md`
- L1 verification/wad_mangrove@wad
  - L5 Build variants (code*/)
  - L12 Input variants (input*/)
  - L27 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/wad_mudflat@wad.md`
- L1 verification/wad_mudflat@wad
  - L5 Build variants (code*/)
  - L12 Input variants (input*/)
  - L30 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/wad_mudflat@wadcheckin.md`
- L1 verification/wad_mudflat@wadcheckin
  - L5 README (first 40 lines) — `addMass`, `useRealFreshWaterFlux`, `addMass.bin`
  - L50 Build variants (code*/)
  - L57 Input variants (input*/)
  - L68 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/wad_overflood@wad.md`
- L1 verification/wad_overflood@wad
  - L5 Build variants (code*/)
  - L12 Input variants (input*/)
  - L23 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/wad_thacker_1d@wad.md`
- L1 verification/wad_thacker_1d@wad
  - L5 Build variants (code*/)
  - L11 Input variants (input*/)
  - L21 Reference results

## `code_util/claude_code/mitgcm-ecco/references/mitgcm-index/verification/wad_thacker_1d@wadcheckin.md`
- L1 verification/wad_thacker_1d@wadcheckin
  - L5 README (first 40 lines) — `bathy_thacker.bin`, `eta_thacker.bin`
  - L50 Build variants (code*/)
  - L56 Input variants (input*/)
  - L66 Reference results

## `code_util/claude_code/mitgcm-ecco/references/troubleshooting/README.md`
- L1 MITgcm troubleshooting history
  - L28 How far to trust an entry
  - L40 Raw material and updating

## `code_util/claude_code/mitgcm-ecco/references/troubleshooting/adjoint.md`
- L1 Troubleshooting: adjoint / TAF / Tapenade / ctrl / cost / grdchk (MITgcm support history 2003-2026)
  - L5 TAF errors / warnings and AD build
    - L7 TAF ERROR "keyword REC, KIND, or SHAPE expected" on `CADJ STORE ... kind = isbyte` (TAF -f08 / TAF_FORTRAN_VERS='F08')
    - L13 TAF ERROR "cannot generate correct recomputations for ..." / "unresolvable conflict" then run dies with "!!!!!!! PANIC !!!!!!! in S/R BARRIER myThid = 0"
    - L19 TAF errors only after upgrading TAF: "TAF ERROR ... store directive for non-existing variable" / "storage directive not initialised"
    - L25 TAF ERROR "code of function tbar not seen and no flow information set, cannot generate derivative common block variable" (new cost term)
    - L31 `make adall`: "staf: Command not found" / ad_taf_output.f not generated
    - L37 TAF errors on mixed fixed/free-format source; `-nocat4ad` (-ncad) breaks on *.flowdir / macOS `.for` files — `TAF_FORTRAN_VERS`
    - L43 `make -j N adtaf` intermittently fails comparing ad_config.template / AD_CONFIG.h (race)
    - L49 TAF-generated tangent-linear model fails at run time: "EXF_INTERP_READ: filename ???? File does not exist" (USE_EXF_INTERPOLATION)
    - L55 Link error "undefined reference to adexch_uv_3d_rl_ / adexch_xy_rs_" (and mdthe_main_loop_) in the adjoint
    - L61 TAF "*ERROR* : variable precip must be available" in a customised ctrl_map_ini_genarr.F
    - L67 Link error "relocation truncated to fit: R_X86_64_32S" building mitgcmuv_ad (or "relocation truncated" in ifort libs) — `nchklev_1`
    - L73 Old checkpoints only: "duplicate symbol _adams_bashforth2_" / "undefined option -nonew_arg" / wrong adjoint after a TAF default change
  - L79 Recomputation, tapes and store directives
    - L81 "TAF RECOMPUTATION WARNING" (any): how the core team handles it
    - L87 Odd result: `if (c) call A; if (c) call B` recomputes but `if (c) {call A; call B}` does not
    - L93 Inner (checkpoint level 1) tape computation unstable / wrong with nonlinFreeSurf > 2 (r*), calc_r_star "too SMALL rStarFac"
    - L99 mitgcmuv_ad blows up (CALC_R_STAR too SMALL) 16 days into the REVERSE run, forward run fine; or garbage gradients far from tape start
    - L105 Adjoint run blows up (adxx ~1e16, finite-difference fine, optim "linesearch failed")
    - L111 Wrong-looking gradients with GMREDI_WITH_STABLE_ADJOINT (stableGmAdjTap): tape key in gmredi_slope_limit
    - L117 Run dies "ADWRITE: tapeFileCounter > tapeMaxCounter 691 690" (tapelev3_..._gptrnm1)
    - L123 Wrong adxx (all negative / flipped sign) with divided adjoint + ALLOW_AUTODIFF_WHTAPEIO + ALLOW_ADCTRLBOUND
    - L129 divided adjoint: "end-of-file during read, unit 76, file divided.ctrl" at the end of the run
    - L135 ecco_phys.F recomputation after adding BEGIN_MASTER/END_MASTER (multithreaded safety) with ALLOW_PSBAR_STERIC
    - L141 Major recomputations from pkg/ecco (ecco_phys.F below line ~141) when NONLIN_FRSURF is undefined
    - L147 pkg/bling AD: wrong store keys / recomputations (bling_light.F, mixedlayer, bio_nitrogen), PHYTO_SELF_SHADING crashes
    - L153 Sea-ice AD: seaice_dynsolver called from seaice_model_ad; useHB87stressCoupling + SEAICE_deltaTdyn>deltaT; SEAICE_USE_GROWTH_ADX
    - L159 Convective adjustment (isomip) gradient depends on unrelated code (store of theta/salt at k-1)
    - L165 Adjoint 10x slower with an advection or physics option (GAD_SMOLARKIEWICZ_HACK, RBCS, PCELL_MIX_CODE)
    - L171 SOLVE_DIAGONAL_LOWMEMORY / SOLVE_DIAGONAL_KINNER with AD (implicit solvers)
    - L177 AD arrays real*4 mismatch: floating overflow in update_cg2d_ad.f with `-ur4` and cg2dFullAdjoint=.TRUE.; adread/adwrite type error with `-warn all`
    - L183 SIGSEGV in adjoint only (forward fine): ctrl_map_genarr2d_ad / adthe_main_loop / profiles, "address not mapped"
    - L189 Segfault in adjoint on 2 nodes with ifort+OpenMPI only (tutorial_tracer_adjsens): fixed with AUTODIFF_USE_OLDSTORE_3D
  - L195 Performance / memory / disk
    - L197 Adjoint runs on Lustre (Pleiades) slow and variable; tapes dominate I/O — `adTapeDir`, `tmpDir`
    - L203 mdsioLocalDir + control pack/unpack fails ("attempt to read/write past end of record" in mdsio_gl.f, MDSREADFIELD_3D_GL)
    - L209 optim.x killed "Exit code -5" on large ecco_cost files
  - L215 ctrl / cost set-up
    - L217 Run stops "MDS_READ_FIELD ... Non-existing record number" xx_atemp.effective.NNNN.data when restarting (nIter0 > 0) with time-varying gentim2d controls
    - L223 gentim2d sensitivities appear in adxx_<name>.effective, adxx_<name> is zero
    - L229 Time-varying ctrl interpolated backwards between records when useCAL=.FALSE.
    - L235 exf field with fldPeriod = 0 plus xx_gentim2d control: control has no effect or accumulates
    - L241 Control adjustment non-zero where control weight is zero; mixed scaled/unscaled adxx (genarr3d, doscaling)
    - L247 Compile error ctrl_getobcs[nsew].F when ALLOW_OBCSN_CONTROL defined but ALLOW_OBCS_NORTH is not
    - L253 Mixing exf atmospheric ctrl flags (ALLOW_ATEMP_CONTROL, ...) with ALLOW_GENARR3D_CONTROL: exf controls silently not applied / missing xx_uwind.effective
    - L259 Scalar (0-D) control needed
    - L265 CPP options defined in several header files give order-dependent behaviour (ALLOW_PACKUNPACK_METHOD2, ALLOW_GENCOST_CONTRIBUTION)
    - L271 gencost time mask: "end-of-file during read" at end of forward run; cost term is zero
    - L277 ADJ* files not written; or wrong snapshot sign; adEtaN zero; ADJ values not "per unit" — `adjDumpFreq`
    - L282 viscFacInAd has no effect — `viscFacInAd`
    - L288 Bottom drag / user control gives zero sensitivity
    - L294 ecco gencost: m_horflux_vol_N suffix disables the cost; objf_gencost array order; boxmean wrong weights
    - L300 Observation records misaligned: cost_gencal assumes obs start with run start; single multi-year gencost_datafile treated as climatology
    - L306 Cost = time-average of a field (SSH, box) in the pre-gencost code
  - L311 Gradient check, TLM, profiles
    - L313 grdchk summary RMS of ratios bad (1e-3 ... 1e-1) although AD is right — `grdchk_eps`
    - L319 Gradient check fails with r* (time-varying hFacC) on the 2nd perturbation / central difference
    - L325 grdchk_main: recglo ignored; grdchk output missing with useSingleCpuIO; gradient index depends on data.ctrl order
    - L330 TLM monitor `g__dynstat_*` all zeros; WARNING CTRL_INIT_CTRLVAR "could not find g_xx_theta ... initialize this file with all zeros"
    - L336 profiles: gradient check wrong with NetCDF model equivalents (profilesDoNcOutput); segfault in profiles_interp_ad
    - L342 Salt (or other) sensitivities all zero when nSx > 1 (e.g. nSx=6,nPx=4) in global_oce_cs32 ad.sens, fine with nSx=1,nPx=6
  - L348 Optimisation (optim / lsopt / optim_m1qn3)
    - L350 optim_m1qn3: omode=5/4/6, terminates early or "nsim/niter changed"; line search "mlis3" stalls
    - L356 optim.x EOF / NaN on Cray (ARCHER): "End of file" in optim_readdata.f
    - L362 Bundled lsopt/optim: build errors, "ifail = 4 the search direction is not a descent one"
  - L368 Tapenade (AD tool for pkg/tapenade, genmake2 -tap)
    - L370 genmake2: Tapenade run fails "flow_tap" not found when the build dir is outside the MITgcm tree
    - L376 Tapenade Makefile rebuilds everything / race / unusable on macOS (case-insensitive FS)
    - L382 Tapenade: "System: Not a tree operator" for every file / FortranParser crash
    - L388 Tapenade streamice test/adjoint 50x slower or wrong; fixed-point; PETSC compile errors
    - L394 Tapenade/TAF compile: CTRL_SIZE.h / CTRL_DUMMY.h included twice (ALLOW_GENTIM2D_CONTROL), wrong _RL vs _RS in exch_tap_b.F
    - L400 OpenAD is being removed (issue 971); build questions
  - L406 ECCO-specific notes
    - L408 ECCO v4 (llc90) "rStarFac too SMALL" after altering wind forcing
    - L414 ecco adjoint verification set-ups (what to run)
    - L419 Monthly/averaged adjoint output with ECCO v4

## `code_util/claude_code/mitgcm-ecco/references/troubleshooting/build.md`
- L1 Troubleshooting: Build / compile / link / start-up (genmake2, optfiles, MPI, NetCDF, TAF/Tapenade)
  - L4 genmake2 / Makefile basics
    - L6 "These lines are here to deliberately cause a compile-time error" (chksum_tiled.f, "Non-numeric character in statement label")
    - L12 Edited/new header or source in code/ is ignored; "This name does not have a type" for a new/renamed variable (diagMdsDir, sideDragFactor, viscC4smag, inAdMode, MNC_TAG_ID); "OBCS.h: No such file"
    - L18 Fortran errors: "_d" / "WORDLENGTH" not defined, "Missing kind-parameter" `0.5 _d 0`, "Illegal use of symbol d0 - KIND", unexpanded `_RL`/`_EXCH_XYZ_R8`, "missing terminating ' character" — `_EXCH_XYZ_R8`, `_EXCH_XYZ_RL`
    - L24 "Error: No Fortran compilers were found in your path" / "can't read OPTFILE=.../linux_aarch64_mpif77|darwin_arm64_mpif77|irix64_ip27_g77" / "no options file found that matches this platform" — `MITGCM_OF`
    - L30 "gfortran: error: unrecognized argument in option '-mcmodel=medium'" on aarch64 Linux (Docker/VM on Apple M-series)
    - L36 "Warning: FORTRAN compiler test failed" / "Makefile might be unusable" though everything is fine
  - L42 MPI / NetCDF library detection and linking
    - L44 `makedepend: warning ... cannot find include file "mpif.h"` / `EESUPPORT.h: error: mpif.h: No such file` / `MPI_SUCCESS has no type` / `MPI_REAL has no type` — `MPI_HOME`
    - L50 "MPI_HOME is not set and pkg-config not available, aborting" (darwin_amd64/arm64_gfortran, linux_amd64_gfortran with -mpi)
    - L56 Link step: `undefined reference to mpi_wait_ / mpi_recv_ / mpi_init_` / `_gfortran_st_write` / `multiple definition of main (for_main.o)` / `cannot find -lnetcdf`
    - L62 genmake2: "Can we create NetCDF-enabled binaries... no" / WARNING: the "mnc" package was enabled but tests failed ... now DISABLED — `NETCDF_ROOT`
    - L68 macOS/Homebrew arm64: `ld: library 'netcdf' not found` or NetCDF check silently disabled (nf-config gives `-lnetcdff -lnetcdf -lnetcdf`, two Cellar prefixes)
    - L74 Run time: `error while loading shared libraries: libnetcdff.so.5 | libifport.so.5 | libmpi_f77.so.1` / `symbol lookup error ... undefined symbol: __intel_sse2_strlen` / `Trace/BPT trap` after NetCDF upgrade
  - L80 Large executables
    - L82 `relocation truncated to fit: R_X86_64_PC32 / R_X86_64_32S against symbol ..._ defined in COMMON section` ; `failed to convert GOTPCREL relocation; relink with --no-relax` ; macOS `ld: 32-bit RIP relative reference out of range`
    - L88 Compiles fine, runs out of memory at load ("There is not enough memory", "Out of memory") with diagnostics — `numLevels`, `numDiags`
  - L94 Compiler-specific
    - L96 gfortran >= 10: `Error: Type mismatch between actual argument at (1) and actual argument at (2)` (MPI_RECV, MPI_Allreduce in gather/global_max/cumulsum)
    - L102 clang/gcc >= 10 on macOS: `implicitly declaring library function 'sprintf'` in setdir.c; `clang: error: invalid version number in '-mmacosx-version-min=11.2'`
    - L108 macOS (case-insensitive FS): "Your file system cannot distinguish between *.F and *.f" ; gfortran "linker input file unused because linking not done" for .fr9 free-form files
    - L114 Wrong results / NaN / segfault only with optimization or only with one compiler (ifort 8-13, PGI 13.4, ifort vs g77) — `__fsd_exp`
    - L120 Model blows up within a few steps with ifort (or at certain tile sizes) but not gfortran; results differ with tile size — `GLOBAL_SUM_ORDER_TILES`
  - L126 Namelists, run directory, start-up
    - L128 Fortran runtime: "namelist not terminated with / or &end", "Cannot match namelist object name X", "End of file" in ini_parms / packages_boot / ptracers_readparms / mnc_readparms
    - L134 scratch1.0000#### / scratch2.0000#### files left in run dir; "Permission denied trying to open file /tmp/gfortrantmp..."; "namelist read ... end of file reached without finding group" (many MPI ranks)
    - L140 `INI_PROCS: needs MPI for multi-procs (nPx*nPy=N) ... usingMPI=False` / `EEBOOT_MINIMAL: No. of procs= 1 not equal to nPx*nPy` / `S/R LOAD_GRID_SPACING ABNORMAL END` / identical output on every tile
    - L146 `forrtl: severe (36): attempt to access non-existent record, unit 9, file bathy.bin|topog.*|OB*.bin` / "Non-existing record number" / "do_ud: end of file" (MDS_READ_FIELD) — `externForcingCycle`
    - L152 Binary input/output looks like garbage / NaN / "inf"; python-written files do not match
    - L158 STDOUT/STDERR empty or not flushed after crash or queue-kill; "output.txt empty" under MPI
    - L164 Run script needs to detect failed runs; no non-zero exit code on crash
    - L170 Segmentation fault at/near start (before first step, in a package init, with huge arrays) — `OMP_STACKSIZE`, `GOMP_STACKSIZE`, `KMP_STACKSIZE`
    - L176 "MON_SOLUTION: STOPPED DUE TO EXTREME VALUES OF SOLUTION" / "SOLUTION IS HEADING OUT OF BOUNDS"; "EEDIE: Only 0 threads have completed"
  - L182 Performance
    - L184 MPI time dominated by GLOBAL_SUM_TILE (MPI_Allreduce of a mostly-zero vector) at many tiles (e.g. llc540, 2819 ranks)
    - L190 Run much slower than expected: huge Wall vs User time, monitor/file I/O dominates — `chkptFreq`, `monitorFreq`, `globalFiles`
  - L196 Packages / CPP-option combinations
    - L198 "*** ERROR *** from PACKAGES_CHECK: run-time control flag useX is set but pkg/X was not compiled" ; "CONFIG_CHECK: #undef NONLIN_FRSURF and nonlinFreeSurf is non-zero" ; "DIAGSTATS_SET_REGIONS: #define DIAGSTATS_REGION_MASK missing in DIAG_OPTIONS.h" — `sizRegMsk`
    - L204 Compile error when `ALLOW_CTRL_OBCS{N,S,E,W}` defined but `ALLOW_OBCS_{NORTH,SOUTH,EAST,WEST}` undefined (OBNt/OBNu undeclared in ctrl_getobcsn.F); obcs_cost_driver.F missing useObcsCostContribution
    - L210 Type-mismatch compile errors with `genmake2 -devel` / `-ur4` (use _RS=real*4) when writing or updating packages — `DIAGNOSTICS_FILL`, `DIAGNOSTICS_FILL_RS`, `useDiagnostics`
    - L216 Compiling pkg/seaice without pkg/exf: "seaice_growth.f ... snowPrecip / EVAP / PRECIP has no type"; STOP in seaice_check "need to define pkg/exf ALLOW_ATM_TEMP"
    - L222 F90 module code: "Can't open module file X.mod" / compiled in wrong order with `make -j` / free-format source rejected by a fixed-format compiler
    - L228 KPP_ESTIMATE_UREF: dbloc / work1 not declared in kpp_forcing_surf.F; LOG of negative zref in vermix
    - L234 NetCDF-less build: pkg/obsfit / pkg/profiles / mnc fail without NetCDF — `useMNC`
    - L240 mnc output wrong grid (X,Y = 1..N) with pkg/exch2 and a lat-lon grid
    - L246 How to use blank tiles (exch2) - blankList / dimsFacets — `blankList`, `dimsFacets`
  - L252 Adjoint / TAF / Tapenade / OpenAD builds
    - L254 TAF adjoint link: `undefined reference to mdthe_main_loop_` (or `adthe_main_loop_`) in the_model_main.o — `mdthe_main_loop`, `adthe_main_loop`, `AD_TAF_FLAGS`
    - L260 TAF: `undefined reference to adexch_3d_rl_ / adexch_uv_3d_rl_ / adexch_uv_xy_rs_ / adexch_xy_rs_` (addummy_in_stepping.o, monitor_ad.o, copy_ad_uv_outp.o) — `AUTODIFF_EXCLUDE_ADEXCH_RS`
    - L266 TAF errors: `CADJ STORE ... *ERROR* identifier not defined` (StoreDynVars3D, tapelev4); `keyword REC, KIND, or SHAPE expected` with -f08 — `USE_EXTENDED_SRC`
    - L272 Parallel make of adjoint: race with `make -j N adtaf` (ad_config.template removed by another make instance)
    - L278 Tapenade: `undefined reference to the_main_loop_b_`; "docker: command not found" fallback; hangs/ fails at fortranParser
    - L284 Tapenade makefile: "flow_tap" file not found; TOOLSDIR undefined in adjoint_tap
    - L290 OpenAD (old): `make adAll` redoes everything / dontTransform, dontCompile files ignored
  - L296 Misc build-tool notes
    - L298 `makedepend: error: out of space: increase MAXFILES` (adjoint/ECCO builds with thousands of files)
    - L304 Harmless genmake2 / compiler messages
    - L309 Windows: genmake2 "No Fortran compilers", "C pre-processor failed the test case", sigreg.c ucontext.h missing
    - L315 "make: ... multiple definition" / duplicate objects after cvs/git update; or `ln: ./CVS: Operation not permitted` in testreport

## `code_util/claude_code/mitgcm-ecco/references/troubleshooting/diagnostics.md`
- L1 Troubleshooting: pkg/diagnostics, output formats (MDS/MNC/NetCDF), budget closure
  - L4 Budget closure: tracers (heat / salt / passive)
    - L6 Heat/salt budget does not close at k=1 with linear free surface (nonlinFreeSurf=0, no r*)
    - L12 TOTTTEND / TOTSTEND do not close the budget with nonlinear free surface / r* (ECCO v4 style), or any run where hFac varies
    - L18 Which surface-flux diagnostic closes the surface heat budget (TFLUX vs oceQnet vs surForcT)
    - L24 Shortwave heating missing in the vertical heat balance in the top ~100 m
    - L30 Subsurface budget residual: ADV*/DF* are fluxes (degC m3/s), plus implicit, KPP and GM terms
    - L36 DFrE_TH non-zero although implicitDiffusion=.TRUE. (or 'DFrE_TH has not been filled')
    - L42 WVELTH vs WTHMASS differ when staggerTimeStep=.TRUE.
    - L48 ADVx_TH / ADVy_TH / ADVr_TH off by a factor of the time step with SOM advection (scheme 80/81)
    - L54 ADVr_TH / DFrE_TH are identically zero at k=1 but WVELTH is not
    - L60 Divergence of u not zero / volume budget with nonlinear free surface
    - L66 DIC / passive-tracer budget does not close (ptracers, DIC, darwin)
    - L72 Restoring (RBCS) / geothermal / SW forcing not visible in budget
    - L78 ocean-to-sea-ice heat flux: which SI* diagnostic
  - L84 Budget closure: momentum / energy
    - L86 Momentum budget residual (flux-form or vector-invariant)
    - L92 Momentum budget off by ~10% in non-hydrostatic runs
    - L98 Um_Diss misses bottom drag / Smag-3D, UBotDrag diagnostic absent
    - L104 Kinetic-energy budget / momKE definition
  - L110 Output timing: frequency, timePhase, calendarDumps
    - L112 Snapshots come out at the middle of the interval (e.g. t=200,600,1000 instead of 0,400,800)
    - L118 Output at unexpected times / extra files with frequency that is not a multiple of deltaT
    - L124 Extra outputs at iter-1 / iter+1 around the right one after long runs (e.g. U.0016783199, ...200, ...201), or 3 files per dump
    - L130 Start diagnostics output only after spin-up / for a limited window
    - L136 Snapshot of the last time step missing; averaged output at end of run
    - L142 Calendar-month / calendar-year averages (calendarDumps) misbehave: half-day offset, snapshot ignored, wrong timeInterval
    - L148 Average over a restart shorter than the averaging frequency gives wrong mean / diagnostics pickup segfault
  - L154 data.diagnostics: setup errors and warnings
    - L156 "WARNING - diag.# NN : NAME has not been filled (ndiag= 0) ... write ZEROS instead" (stats: "write UNDEF instead")
    - L162 "*** ERROR *** DIAGNOSTICS_SET_POINTERS: <name> is not a Diagnostic"
    - L168 "DIAGNOSTICS_ADDTOLIST: Exceed Max.Number of diagnostics ndiagMax= 500" / numDiags, numperlist, diagSt_size too small / segmentation or memory problems
    - L174 Mixing 2-D (ETAN) and 3-D (UVEL) fields in one output stream gives only level 1 of the 3-D field
    - L180 How to add your own diagnostic
    - L186 "DIAGNOSTICS_STATUS_ERROR ... DIAGNOSTICS_FILL called from the WRONG place" (after the last DIAGNOSTICS_WRITE in DO_THE_MODEL_IO)
    - L192 Namelist read failures that look like diagnostics bugs (end-of-file, unit 11, crashes only with useDiagnostics=T)
    - L198 Output precision: diagnostics file is 64-bit when you want 32-bit (or reverse); NetCDF in float
  - L204 Statistics diagnostics, regions, sub-domain output
    - L206 Regional / single-point statistics (DIAG_STATIS_PARMS, region masks)
    - L212 Cannot write a sub-region or vertical section with pkg/diagnostics; need only point/section output
    - L218 Snapshot adjoint diagnostics (ADJtheta etc.) have the wrong sign for frequency < 0
  - L224 MDS / MNC / NetCDF output
    - L226 useMNC=.TRUE. but "pkg/mnc has not been compiled (#undef ALLOW_MNC)", or genmake2 "mnc package was enabled but tests failed to compile NetCDF applications ... DISABLED"
    - L232 MNC writes one NetCDF file per tile; useSingleCpuIO / globalFiles have no effect on NetCDF
    - L238 MNC restart/pickup errors: "MNC ERROR: variable 'RC' is already defined ... different grid shape", "dimension 'Z'/'Zp1' does not exist", "cannot open either a per-face or a per-tile file"
    - L244 MNC output slower than MDS (up to 10x) or produces many files
    - L250 MNC output has Xp1/Yp1 one longer than X/Y (halo copy), XG/YG last column duplicates the first
    - L256 MNC missing values (useMissingValue) mask ice-shelf or surface diagnostics wrongly
    - L262 MDS: one file per tile (S.0000000000.001.001.data ...), joinmds fails, restart with different processor count
    - L268 Grid files (grid.*.nc, XC/YC .data) not written
    - L274 rdmds / MITgcmutils errors (ValueError parsing .meta; NaN import; cannot read missing .meta; dimension reshape)
    - L280 Reading MDS output in other tools; offline-run file naming
  - L286 Interpreting specific diagnostics
    - L288 PHIHYD / PH / PHIBOT / PHI_NH units and how to get pressure
    - L294 Buoyancy frequency / stratification from output with a nonlinear EOS
    - L300 Bolus overturning from GM_PsiX/PsiY (and GM_Kwy) is non-zero at depth
    - L306 EXFuwind/EXFvwind are not the 10-m winds with useRelativeWind=.TRUE.
    - L312 EXF diagnostics constant in time / ERA5 forcing sanity (EXFswdn, EXFlwdn, EXFatemp)
    - L318 pkg/layers diagnostics are all zero or have confusing names/units
    - L324 Mixed layer depth diagnostics: MXLDEPTH vs KPPhbl
    - L330 Frictional heating diagnostic (HeatDiss) is all zeros in the ocean
  - L336 Misc
    - L338 Removed: pkg/timeave and CPP EXACT_CONSERV
    - L344 DIC pkg: dic_pCO2 error; surface-only diagnostics; OMEGAC zero; 3-D silicate
    - L350 Mellor-Yamada (pkg/my82) / MY82 diagnostics location and units
    - L356 Huge STDOUT, too many text files, memory and link errors
    - L362 Zero output / zero depth because of wrong input files, grid sign or SIZE.h mismatch
    - L368 Anelastic option and RhoRef.data

## `code_util/claude_code/mitgcm-ecco/references/troubleshooting/forcing.md`
- L1 Troubleshooting: surface forcing (pkg/exf, cal, restoring, runoff)
  - L4 EXF set-up: CPP combinations, namelist errors, unrealistic heat flux
    - L6 Temperature goes wild / hflux far larger than input, e.g. exf_hflux_max 1200 when the file max is 679 (EXF_OPTIONS.h flags inconsistent with forcing files) — `ALLOW_ATM_TEMP`, `ALLOW_DOWNWARD_RADIATION`, `ALLOW_BULKFORMULAE`
    - L12 Ocean cools indefinitely (-20 C coastal, lake, "simple EXF setup") with downward radiation missing or wrong
    - L18 STOP "EXF_CHECK: use u,v_wind components but not wind-stress" (or the reverse) — `ALLOW_ATM_WIND`, `EXF_ATM_WIND`
    - L24 Namelist error "variable not in namelist EXF_NML_0x" / "read unexpected character" / "Fortran runtime error: syntax error in NAMELIST" — `USE_EXF_INTERPOLATION`, `_lat_inc`
    - L30 "EXF_INTERP_UV: input grid must encompass output grid" / striped or distorted forcing / "Non-existing record number" with interpolation — `EXF_INTERP_USE_DYNALLOC`
    - L36 Bicubic vs bilinear in exf_interp: spurious values, rectangular patterns in wind-stress curl
    - L41 EXF values wrong with compiler optimisation on (atemp of constant 300 K shows max 1200)
    - L46 Constant-in-time forcing with exf; stress read, wrong grid — `readStressOnAgrid`, `zonalWindFile`, `periodicExternalForcing`
    - L51 Wrong values in generic (non-exf) forcing: taux0/tauy0 empty, fu/fv constant — `periodicExternalForcing`
    - L57 Forcing in non-exf mode is linearly interpolated, "suddenly changing" wind not possible
    - L62 exf relative wind: is ocean current subtracted for latent/sensible heat too?
    - L67 exf requires pkg cal (build fails with missing cal headers; STOP CAL_FULLDATE called too early)
  - L73 Calendar, start dates, periods, restarts
    - L75 Model reads the wrong record / crashes with "Fortran runtime error: Non-existing record number" (forcing file too short)
    - L81 Monthly/climatological forcing: period -12, mid-month start, leap years, "how long is a year" — `file_YYYY`, `cal_isleap.F`, `nDaysLeap`
    - L86 Starting year/phase issues: restart does not read forcing from the right time; repeatPeriod, yearly fields, very long runs — `repeatPeriod`, `pickupSuff`, `DIC_forcingCycle`
    - L91 calendarDumps / month-end output and pickups — `rn_pickup`, `dumpInitAndLast`, `calendarDumps`
    - L96 EXF fldPeriod=0 with ctrl xx_gentim2d; ctrl time-varying weights backwards with useCAL=F — `CTRL_GET_GEN_REC`
    - L102 Interannual monthly forcing broke OBCS/EXF in checkpoint68d
    - L108 TAF TL code with USE_EXF_INTERPOLATION: "EXF_INTERP_READ: File does not exist" — `ALLOW_AUTODIFF_MONITOR`
  - L113 Heat and freshwater fluxes, runoff, free surface
    - L115 Negative salinity / crazy salinity (-200..90) at coastal cells or big river runoff — `EXACT_CONSERV`, `temp_EvPrRn`
    - L121 Runoff: how to put river water in; at depth; heat content of runoff — `runoffFile`, `runoff_lat_inc`, `ALLOW_ADDFLUID`, `addMassFile`, `runOffMapFile`
    - L126 Prescribing heat/salt flux (hflux, sflux) with seaice on — `EXF_READ_EVAP`
    - L132 Net heat flux sign and components (Qnet, Qsw, shortwave penetration) — `surfQnetFile`, `ALLOW_BULKFORMULAE`, `balanceEmPmR`, `balanceQnet`, `ALLOW_BALANCE_FLUXES`
    - L137 Freshwater flux, sea-ice melt vs exf budget (docs gap) — `useRealFreshWaterFlux`, `convertFW2Salt`, `temp_EvPrRn`, `exactConserv`
    - L142 Evaporation/precip input precision: salt drifts although sflux integrates to zero
    - L148 Wrong friction velocity ustar in exf_wind.F / wrong stable-case term psimh in exf_bulkformulae.F (old code) — `exf_wind.F`, `exf_bulkformulae.F`
    - L153 Frictional heating (ALLOW_FRICTION_HEATING) gives zeros in ocean — `frictionHeating`
  - L159 Restoring, RBCS, sponge, relaxation
    - L161 Surface SST/SSS restoring: set-up, time scale, regional restoring, restarts — `climsstTauRelax`, `climsssTauRelax`, `tauThetaClimRelax`, `thetaClimFile`, `latBandClimRelax`
    - L166 pkg/rbcs: restoring velocities, tracers, ice; maskLEN; time-scale; restart offset — `maskLEN`, `rbcsForcingPeriod`, `rbcsForcingCycle`, `rbcsForcingOffset`, `rbcsIniter`
    - L171 Restoring tracer RBCS: asymmetric solution under symmetric forcing
    - L177 3D forcing/tendency in momentum or tracer equations (body force, drag, tidal body force, bottom heat flux) — `external_forcing.F`, `USE_OLD_EXTERNAL_FORCING`, `MYPACKAGE_TENDENCY`
    - L182 Tidal forcing options — `tidePot`, `useOBCSprescribe`, `useOBCSbalance`, `ATMOSPHERIC_LOADING`
  - L187 OBCS interactions (forcing-related)
    - L189 OBCS + EXF: STOP EXF_CHECK_RANGE after some time, or noise at open boundary
    - L195 Prescribed OBCS time control (with and without exf); yearly OBCS fields — `pickupSuff`
    - L200 OBCS sponge not tile-proof; domain decomposition changes results
    - L206 Periodicity assumptions with walls and OBCS — `maskInC`, `notUsingXPeriodicity`
  - L211 Input files, binary format, tools
    - L213 Input file format: always unblocked IEEE big-endian; precision; orientation
    - L218 'make' after header changes / "relocation truncated to fit" / makedepend out of space
    - L223 Forcing is read but nothing happens / ocean at rest / no response — `zonalWindFile`, `periodicExternalForcing`, `SHELFICEloadAnomalyFile`
  - L228 pkg/cheapaml, aim, thsice forcing (few entries)
    - L230 CheapAML on llc90 / with seaice
    - L235 aim_v23 with thsice: SST from file, orbital parameters, insolation — `aim_useFMsurfBC`, `aim_useMMsurfFc`, `ALLOW_INSOLATION`
  - L240 Misc forcing physics Q&A with clear answers
    - L242 Adding mass/evaporation without salinity change (addmass)
    - L247 Atmospheric pressure: sign/offset in EOS; setting p_top > 0

## `code_util/claude_code/mitgcm-ecco/references/troubleshooting/grid.md`
- L1 Troubleshooting: grids, bathymetry, hFac, vertical coordinates, r*
  - L4 Input files, namelists, rebuilds
    - L6 Custom `CPP_OPTIONS.h` (or any header in `code/`) ignored; CONFIG_CHECK says `#undef NONLIN_FRSURF`
    - L12 `S/R LOAD_GRID_SPACING: No value for delX at i = N` (or delY) after changing resolution
    - L18 `S/R LOAD_GRID_SPACING: delX(i=k)= 0.0 : MUST BE >0` for a block of i — `readBinaryPrec`, `load_grid_spacing.F`
    - L24 `S/R INI_VERTICAL_GRID: neither delR nor delRc are defined` / `S/R INI_PARMS: No value for delZ/delP/delR at k = N` — `ini_parms.F`, `ini_vertical_grid.F`
    - L30 `namelist read: read unexpected character` / `invalid reference to variable in NAMELIST input` / `Cannot match namelist object name` — `useMDSIO`, `debugMode`
    - L36 Input files opened read-write; fails on read-only filesystem — `_READONLY_ACTION`, `EXCLUDE_OPEN_ACTION`
    - L42 Immediate segfault / `Killed` at startup of a large setup
  - L48 Bathymetry, masks, binary I/O
    - L50 Depth.data all zeros / no ocean; or model depth = input depth + constant — `Ro_SeaLevel`, `top_Pres`, `seaLev_Z`
    - L56 Segfault, `forrtl: error (65): floating invalid` in `ini_masks_etc`, or run dies at once, after editing bathymetry / T,S init
    - L62 Bathymetry/forcing map appears transposed/rotated, wrong number of points
    - L68 Initial/forcing fields interpolated to a finer land-sea mask have zero/garbage at coast
    - L74 Single-row (2-D x-z) domain is periodic in y; need varying cross-section width — `no_slip_sides`, `addSwallFile`, `addWwallFile`, `add_walls2masks.F`, `ini_cartesian_grid.F`
    - L80 T,S,U,V = 0 at the bottom level / "bottom" output is zero
  - L86 Horizontal grids: Cartesian, spherical, rotated, curvilinear
    - L88 Wind stress / forcing on a rotated or curvilinear grid points the wrong way; UVEL is not "eastward" — `zonalWindFile`
    - L94 `usingCurvilinearGrid=.TRUE.`: `file not found tile002.mitgrid`, every proc reads tile001, or `attempt to access non-existent record` — `horizGridFile`, `tileNNN.mitgrid`, `ini_curvilinear_grid.F`
    - L100 OLD_GRID_IO file names; dummy/copy metrics (DXF=DXG=DXC...) give noise or blow-up — `ini_spherical_polar_grid.F`, `ini_cartesian_grid.F`
    - L106 `delX/delY/ygOrigin/f0` have no effect once `usingCurvilinearGrid=.TRUE.` — `rotate_spherical_polar_grid.F`, `ini_spherical_polar_grid.F`
    - L112 Where do curvilinear/cs/llc grid files come from? (SPGrid dead; mitgrid with `rAz=0`, `dxC~22 m`)
    - L118 AngleCS/AngleSN (angleCosC/angleSinC) wrong, jump to 1/0, in last row/column (i=Nx, j=Ny) — `CALC_GRID_ANGLES`, `ini_curvilinear_grid.F`, `ini_spherical_polar_grid.F`
    - L124 Locations of U/V points on curvilinear grid; no XU/YU/XV/YV output
    - L130 MNC grid output: last row/column of XG/YG equals the first; gluemnc puts curvilinear tiles side by side — `usingCurvilinearGrid`
  - L136 Cubed-sphere and LLC specifics
    - L138 `forrtl: severe (36): attempt to access non-existent record` reading `grid_csNN.faceNNN.bin` — `ini_curvilinear_grid.F`
    - L144 `EXCH1_RS_CUBE: Wrong Tiling ... works only with sNx=sNy & nSx=6 & nSy=nPx=nPy=1` — `ini_curvilinear_grid.F`
    - L150 Blank tiles / `data.exch2` / "cpus on land points" — `global_ocean.90x40x15`, `SIZE.h`
    - L156 Cubed-sphere velocity / vector output looks broken across tile edges; initializing winds — `uvLatLon2cube.m`, `rotate_uv2uvEN.m`, `cubeCalcAngle.m`
    - L162 Cubed-sphere/LLC unified `.meta/.data` per time step; tiled output merging — `globalFiles`
    - L168 MOC / overturning on cs or llc grid gives nonsense; broken-line file mismatch — `use_isoLat_bkl.m`
    - L174 Non-Earth planet on cubed sphere (grid files assume Earth radius) — `rSphere`, `radius_fromHorizGrid`, `ini_curvilinear_grid.F`
    - L180 `deepAtmosphere`: where is `rSphere`? p* + deepAtmosphere?
  - L186 Vertical grid, partial cells (hFac), cavities
    - L188 Blow-up (`SOLUTION IS HEADING OUT OF BOUNDS`, `EXTREME Pot.Temp`, NaN, cg2d) after enabling partial cells or refining dz — `monitorFreq`, `advcfl_W_hf_max`, `advcfl_wvel_max`, `debugLevel`, `viscAhGrid`, `viscAh`
    - L194 High-resolution horizontal stripes / grid noise; how big is viscAhGrid? — `mom_calc_visc.F`, `viscAhGridMax`, `viscA4Grid`
    - L200 `OBCS_CHECK: Inside Mask and OB locations disagree` with pkg/shelfice (ice draft meets bathymetry at hFacMin*delR) — `ini_masks_etc.F`
    - L206 Thin ice shelf (hFacC(k=1)>0): thsice/seaice forms under shelf, short-wave warms cavity
    - L212 Isolated coastal cells get salinity -200 ... +90 (surface fresh-water/EmPmR forcing) — `nonlinFreeSurf`, `useRealFreshWaterFlux`
    - L218 Coastal currents veer / flow follows shelf wrongly in regional setup
  - L224 r*, nonlinear free surface, SSH
    - L226 `STOP in CALC_R_STAR : too SMALL rStarFac[C,W,S] !` (`WARNING: r*FacC < hFacInf`) — `useRealFreshWaterFlux`, `MAX_HEFF`, `calc_r_star.F`
    - L232 Which free-surface settings are valid? `select_rStar`, `nonlinFreeSurf`, `DISABLE_RSTAR_CODE` — `config_check.F`, `DISABLE_RSTAR_CODE`, `DISABLE_SIGMA_CODE`, `SURF_ADJUSTMENT`
    - L238 Negative / depressed SSH where sea ice is thick; SSH correlates with sea-ice thickness
    - L244 Mean SSH drifts continuously (regional model) — `balanceSaltClimRelax`, `useOBCSbalance`
    - L250 No wetting/drying; SSH lower than the top cell or the whole depth
    - L256 r* with ice-shelf cavities / thick sea ice and JMD95Z EOS
    - L262 r* grad-P ("sloping term") wrong when model top is not at r=0
    - L268 Which tRef/hydrogThetaFile? Kelvin-wave test gives V != 0 at step 1 — `hydrogThetaFile`, `tRef`
  - L274 Budgets, pressure, diagnostics tied to the grid
    - L276 Volume/heat/tracer budget does not close with grid-file hFacW/hFacS
    - L282 PHIHYD / PHL / bottom pressure: values huge, topography imprint, PH != pressure — `totPhiHyd`, `REMOVE_MEAN_RL`
    - L288 Closing the momentum budget in LLC4320 with self-computed hydrostatic pressure (noisy stripes) — `Um_dPhiX`, `Vm_dPhiY`
    - L294 GM_Psi (bolus streamfunction) weighting with hFac and deepAtmosphere — `gmredi_residual_flow.F`
    - L300 pkg/layers: overturning in density coordinates in p-coordinates / wrong reference pressure; budget with LinFSConserveTr=F — `layers_krho`, `LAYERS_THERMODYNAMICS`, `LAYERS_UFLUX`, `LAYERS_VFLUX`, `LAYERS_THICKNESS`
    - L306 SEAICE LSR solver: ice does not move at fine grid spacing (<~500 m), depends on tile edges
    - L312 Tides on a prescribed grid: forcing via ATMOSPHERIC_LOADING instead of editing external_forcing.F — `pLoadFile`
    - L318 Byte-swapped input on little-endian/Pathscale: `-D_BYTESWAPIO`

## `code_util/claude_code/mitgcm-ecco/references/troubleshooting/instability.md`
- L1 Troubleshooting: instability, blow-ups, NaNs, CFL, time step
  - L5 How to diagnose a blow-up — `monitorFreq`, `deltaT`, `debugLevel`, `advcfl_uvel_max`, `advcfl_vvel_max`, `advcfl_wvel_max`, `advcfl_W_hf_max`, `dumpFreq`, `readBinaryPrec`, `useExfCheckRange`
  - L18 CFL, time step, viscosity/diffusion limits
    - L19 SOLUTION IS HEADING OUT OF BOUNDS: tMin,tMax= ...  (ABNORMAL END: S/R MON_SOLUTION, stops due to EXTREME Pot.Temp) — `monSolutionMaxRange`
    - L25 Model fine for N steps then NaN / "tMin,tMax" large, after changing resolution or grid — `viscAhGrid`, `viscA4Grid`
    - L31 viscAhGrid / viscAh details; numerical noise although advective CFL << 1 (e.g. 0.02) — `ivdc_kappa`, `viscAhGridMax`, `viscC2smag`
    - L37 Blow-up with thin partial cells (hFacMin small, hFacMinDr small): vertical CFL — `hFacMin`, `hFacMinDr`, `advcfl_W_hf_max`, `advcfl_wvel_max`
    - L43 Blows up in first hours/days after starting from climatology / interpolated / perturbed initial conditions — `ivdc_kappa`
    - L49 NaNs/blow-up caused by tracer reference profile ordering or EOS parameters (tRef reversed, tAlpha large, extreme salinity) — `tRef`, `sRef`, `tAlpha`, `selectP_inEOS_Zc`
    - L55 STOP ABNORMAL END: S/R INI_THETA / INI_SALT "found N wet grid-pts with theta=0 identically" — `hydrogThetaFile`
    - L61 Immediate NaN / "Initial hydrostatic pressure did not converge after 15 steps" / "Iteration 1, RMS-difference = NaN" / seaice-dynamics NaN after one step
    - L67 Uniform/zero-forcing run still unstable: central advection without diffusion (default tempAdvScheme=2) — `tempAdvScheme`, `diffKhT`, `diffKrT`
  - L73 Advection scheme / time-stepping scheme
    - L74 ** WARNING ** GAD_CHECK: potentially unstable time-stepping (Internal Wave) / need "staggerTimeStep=.TRUE." in "data", nml PARM01
    - L80 Instability/odd 2-D or poles runs fixed by Adams-Bashforth parameter
    - L86 Blow-up / NaN with higher-order scheme and tiny overlap (OLx, OLy too small; sNy=1 2-D)
    - L92 Quasi-hydrostatic + staggerTimeStep fatal instability (large QH terms, weak stratification) — `ALLOW_QHYD_STAGGER_TS`
    - L98 Instability after updating code past checkpoint64e with SOM advection (schemes 80/81)
    - L104 Hydrostatic grid-scale W noise ("columns of alternating w") over slopes/mixed layer — `ivdc_kappa`, `cAdjFreq`
    - L110 Surface-wave / implicit free surface damping or rigid lid blow-ups
    - L116 NH: NaN after minutes with sharp density step / cg3d not converging — `cg2dUseMinResSol`, `cg3dMaxIters`
  - L122 Time-step parameters (deltaT, deltaTmom, deltaTtracer, deltaTClock)
    - L123 Setting deltaTmom/deltaTtracer/deltaTfreesurf/deltaTClock separately: no speed-up, "wrong" results, ice not forming, Kelvin waves too slow — `deltaT`
    - L129 Offline tracer run unstable: deltaTtracer=86400 with daily flow — `deltaToffline`, `deltaTtracer`
    - L135 cal_Set: The time step is less than a second. / contains fractions of a second.
    - L141 Diagnostics/dump output never written at the requested time
  - L147 Surface forcing (EXF), heat/salt/fresh-water
    - L148 ABNORMAL END: S/R EXF_CHECK_RANGE
    - L154 SST 30C+ ("heat trapping") then SURF_ADJUSTMENT / cg2d divergence in summer
    - L160 NaN with cheapaml / negative specific humidity / negative relative humidity
    - L166 exf_interp_uv produces |vector| > 1e6 near the pole
    - L172 Crash with forcing/BC file read: "Non-existing record number" / read past end (looks like a blow-up)
  - L178 Open boundaries (OBCS)
    - L179 Eta/sea level drifts then NaN: unbalanced OBC inflow/outflow (SSH growth, "unrealistic sea surface height growth")
    - L185 Spurious boundary jets / blow-up with OBCS and useCDscheme=.TRUE. — `viscA4`
    - L191 Blow-up at sponge layer / Orlanski boundary / inconsistent lateral boundaries, "horizontal strips in currents"
  - L197 Sea ice
    - L198 STOP in CALC_R_STAR : too SMALL rStarFac[C,W,S] ! / ABNORMAL END: S/R CALC_R_STAR  (very thick ice, big eta excursions) — `hFacInf`
    - L204 MOM_IMPLICIT_R: error when solving 3-Diag problem. (also after sea-ice crash)
    - L210 EVP sea-ice solver explodes (ice velocities blow up at a point, no LSR problem)
    - L216 Runaway ice thickness (HEFF +50 m in one step), large precipitation / snow
    - L222 Huge Qnet (1e4 W/m2) where supercooled ice-shelf water surfaces — `SEAICE_frazilFrac`
    - L228 Sea ice dynamics gives NaN with no ice present (checkpoint65u LSR)
    - L234 Huge heat release to ocean (+25 C) when snow floods, SEAICE_USE_GROWTH_ADX
    - L240 Sea-ice stress divergence missing metric terms (curvilinear grids) — `SEAICEselectMetricTerms`
  - L245 Ice shelf (pkg/shelfice, icefront)
    - L246 shelfice_thermodynamics NaN: freshWaterFlux = NaN when saltFreeze = 0 (sLoc = 0), e.g. llc4320 cavity init
    - L252 Blow-up within minutes after enabling useSHELFICE (regional, OBCS+seaice) — `SHELFICEloadAnomalyFile`, `rhoConst`
  - L258 Adjoint / TAF
    - L259 Adjoint run dies with CALC_R_STAR "too SMALL rStarFac" in tape (forward) computations although forward run is fine — `doSinglePrecTapelev`
  - L265 Input/setup errors that look like instabilities
    - L266 Crash "forrtl: error (65): floating invalid" right after reading an edited pickup / NaN at step 1-2 with new input
    - L272 Segmentation fault / NaN: SIZE.h inconsistent with data
    - L278 Zero ocean velocity/odd values: bathymetry sign, grid latitude — `bathyFile`, `phiMin`
  - L284 Packages
    - L285 aim_v23 / SPEEDY: "AIM_DYN2AIM: Temp out of range 100 400" / "SUFLUX_POST: TS out of range 100 400" after centuries — `Shap_Trtau`
    - L291 NaN in DIC/GCHEM tracers
    - L297 Passive homogeneous salinity blows up after very long run
    - L303 Floats (pkg/flt) segfault when model near blow-up
  - L309 Platform, compiler, NaNs continue
    - L310 Model keeps running writing NaN after blow-up (burning CPU) — `monitorFreq`
    - L316 Note: IEEE_INVALID_FLAG / IEEE_OVERFLOW_FLAG printed at "STOP NORMAL END"
    - L322 Blow-up depends on compiler flags / optimization / tile count (non-reproducible)
    - L328 Build/run: "relocation truncated to fit: R_X86_64_PC32 ... COMMON" / gfortran "-mcmodel=medium" unrecognized (arm64)
    - L334 Job hangs after a rank STOPs (MPI deadlock) / run "dies" but SLURM keeps burning hours — `writePickupAtEnd`
    - L340 Stderr/stdout files: STOP/ABNORMAL END but STDERR empty / "Cannot open file '#dev/null'" with SINGLE_DISK_IO

## `code_util/claude_code/mitgcm-ecco/references/troubleshooting/obcs.md`
- L1 Troubleshooting: pkg/obcs (open boundaries, Orlanski, Stevens, sponge, tides, balance, OB input files)
  - L4 Enabling the package and CPP/namelist mismatches
    - L6 genmake2 "In ../code/CPP_OPTIONS.h there is an illegal line: #define ALLOW_OBCS" / runtime "Run-time control flag useOBCS was used when CPP flag ALLOW_OBCS was unset" — `OBCS_OPTIONS.h`
    - L12 `OBCS_CHECK: Cannot set useOBCSsponge=.TRUE. ... with #undef ALLOW_OBCS_SPONGE` (same for prescribe, tides, Stevens) — `OBCS_OPTIONS.h`, `ALLOW_OBCS_PRESCRIBE`, `ALLOW_OBCS_BALANCE`, `ALLOW_ORLANSKI`, `ALLOW_OBCS_SPONGE`, `ALLOW_OBCS_STEVENS`, `ALLOW_OBCS_TIDES`, `ALLOW_OBCS_SEAICE_SPONGE`, `obcs_check.F`, `obcs_check.f`
    - L18 `Fortran runtime error: End of file` in obcs_readparms / `Cannot match namelist object name` / `"tidalPeriod" no longer allowed in file "data.obcs"` — `tidalPeriod`, `OBCS_tidalPeriod`
  - L24 Where the boundary is: OB_Ieast/Iwest/Jnorth/Jsouth, OBCS_CHECK, corners
    - L26 Meaning and length of OB_Ieast, OB_Iwest, OB_Jnorth, OB_Jsouth
    - L32 `OBCS_CHECK: Inside Mask and OB locations disagree` / `N errors in OB location vs Mask` / `N errors in tile OB set-up` — `ALLOW_OBCS_NORTH`, `insideOBmaskFile`
    - L38 OB point in a domain corner or inside a bay: OBCS_CHECK fails, large w, or SSH drift with four open sides
    - L44 Cannot prescribe a river/"inner" OB inside a domain whose perimeter is all open, or a second OB on the same row — `insideOBmaskFile`
    - L50 Inclined / rotated open boundaries
  - L56 OB input files: dimensions, records, time, precision, defaults
    - L58 `forrtl: severe (36)/ Fortran runtime error: Non-existing record number` in mdsio_read_section (file OBNu...), `MDSREADFIELD_XZ_GL: File does not exist` — `readBinaryPrec`, `mdsio_read_section.F`
    - L64 OB velocities/tracers appear as ±1e32, zeros, or "model never reads OB files" — `readBinaryPrec`, `exf_iprec_obcs`, `exf_iprec`, `useOBCSprescribe`, `ALLOW_OBCS_PRESCRIBE`
    - L70 What are the OB values when you give no files? (excess cooling / relaxation to tRef at the boundary) — `obcs_calc.F`
    - L76 OB time stamps: how records map to time (no EXF vs EXF), monthly/yearly series — `repCycle`, `useOBCSYearlyFields`
    - L82 Which grid point/hFac do OB values refer to? (staggering, inward point, OBCS_uvApplyFac)
    - L88 Bad T/S/u/v at land-ocean transitions along the boundary (e.g. salinity ~0 next to coast)
    - L94 OB inputs on an irregular vertical grid; OB from ECCO/other model output
    - L100 OB eta files: `OB*etaFile(s) only allowed with nonlinFreeSurf`; SSH at the boundary does not drive the interior; OBNhFile
    - L106 Using the internal_wave verification code: constant OB forcing is ignored — `obcs_calc.F`
  - L112 Volume balance and SSH drift
    - L114 Sea level rises/falls by metres to tens of metres in a regional run with prescribed OBs (or Orlanski) — `balanceEmPmR`, `ALLOW_BALANCE_FLUXES`, `OBCS_monitorFreq`
    - L120 What do the OBCS_balanceFac values mean? (`-1`, `0`, `1`, `2`)
    - L126 Balance does not correct non-OBCS volume sources (ice-sheet/shelf volume change, E-P) — `balanceEmPmR`, `inFlow`, `obcs_balance_flow.F`, `OBCSbalanceSurf`
    - L132 `useOBCSbalance=.TRUE.` crash (0/0) with nSx>1 and OBs at domain edge
    - L138 dynstat_eta_max = 0 / eta exactly 0 on open boundary
  - L144 Orlanski radiation
    - L146 `OBCS_CHECK: useOrlanski* OBC not yet implemented for pTracers` / `nonlinFreeSurf not yet implemented in Orlanski OBC` / with seaice / `OBCS not yet implemented in CD-Scheme` — `obcs_check.F`, `useCDscheme`, `obcs_calc.F`
    - L152 Orlanski: reflection of internal waves, no net flow, SSH drift, tides ignored — `cVelTimeScale`, `useOBCSbalance`, `viscC2smag`
    - L158 NaN in one variable (e.g. T) along an Orlanski boundary after long runs (divide by zero) — `tutorial_plume_on_slope`
    - L164 Orlanski restart and OB pickup/section files written per tile — `useSingleCpuIO`
  - L170 Sponge layer
    - L172 How to choose spongeThickness, Urelaxobcsinner/bound; sponge restores to OB values, not to "nothing" — `obcs_sponge.F`, `spongeThickness`, `_UatNS`, `_UatEW`, `_VatNS`, `_VatEW`, `_Theta`, `_Salt`, `useLinearSponge`
    - L178 Sponge bugs fixed upstream: South salinity timescale, T/S relaxation skipped where OB value is 0, east OB `isl-1` — `lambda_obcs_s`, `obcs_apply_ts.F`, `obcs_sponge_v`, `shelfice_2d_remesh`, `obcs_ctrl`
    - L184 Sponge thicker than a tile starts "suddenly" / sponge appears at every tile edge (parallel)
    - L190 Sponge for passive tracers, T-only sponge, "closed" wall + sponge, or sea-ice sponge — `useSeaiceSponge`, `seaiceSpongeThickness`, `ALLOW_OBCS_SEAICE_SPONGE`, `OBCSsponge_Theta`
    - L196 Spurious vertical velocity / waveguide / upwelling along an open boundary
  - L202 Stevens (1990) boundary conditions
    - L204 `useStevens*`: `stops due to EXTREME Pot.Temp` near the boundary; warnings about nonlinFreeSurf/pTracers
  - L210 Tides
    - L212 Setting up barotropic tides at OBs (units, files, nothing happens) — `obcs_add_tides.F`, `advcfl_W`
    - L218 Tides + sponge/daily boundary data; tidal currents too strong in a narrow channel; no SSH forcing — `obcs_add_tides`, `zRoughBot`
    - L224 Internal-tide sponge gives artificial U/T disturbances; OB internal-wave reflection
  - L230 Sea ice at open boundaries
    - L232 Sea ice piles up/gets stuck at OBs; `Orlanski ... useSEAICE`; periodic ice with ocean OBs — `OBCS_SEAICE_AVOID_CONVERGENCE`, `OBCS_SEAICE_COMPUTE_UVICE`, `ALLOW_OBCS_SEAICE_SPONGE`, `seaiceSpongeThickness`, `seaice_obcs.seaiceSponge`, `obcs_adjust_uvice.F`, `obcs_apply_seaice.F`, `obcs_apply_uvice.F`, `obcs_seaice_sponge.F`
  - L238 Passive tracers at open boundaries
    - L240 Passive tracer concentration explodes (1e38, negative) near an OB; what is the default tracer OB? — `so_box_biogeo`
  - L246 Domain, parallel layout, misc
    - L248 2-D slab (sNy=1) with OBCS: T/S go to zero or the run blows up after the first step
    - L254 Default boundaries are doubly periodic; OB values "appear on both sides"; extra Nx+1/Ny+1 column in netCDF output
    - L260 Freshwater (river) through an OB: negative salinity, cg2d diverges at step 1
    - L266 OBCS + EXF: `EXF_CHECK_RANGE` / hflux out of range after enabling OBs; OBCS not combined with CD scheme; time/NaN near boundary
    - L272 Non-hydrostatic + OBCS stop in CONFIG_CHECK / w at the OB — `config_check`, `obcs_apply_w`, `obcs_apply_w.F`
    - L278 One-way nesting / regional ECCO-style downscaling: what pkg/obcs can do — `obcs_ctrl`, `tutorial_plume_on_slope`, `so_box_biogeo`

## `code_util/claude_code/mitgcm-ecco/references/troubleshooting/parallel.md`
- L1 Troubleshooting: MPI / tiling / SIZE.h / exch2 / blank tiles / OpenMP / scaling
  - L4 Launch / "No. of procs" / MPI build
    - L5 EEBOOT_MINIMAL: No. of procs= 1 not equal to nPx*nPy= N  (also "No. of processes not equal to nPx*nPy 12 2", "EEDIE: earlier error in multi-proc/thread setting")
    - L11 SIZE.h edits ignored with -mpi / "SIZE.h_mpi" silently overrides SIZE.h — `SIZE.h_mpi`
    - L17 EXCH-1 useCubedSphereExchange unsafe with usingMPI=.TRUE.
    - L23 EESET_PARMS: Error reading parameter file "eedata" / EEDIE
    - L29 Compile/link: "cannot find include file mpif.h", "mpif77: command not found", "/usr/bin/ld: cannot find -lmpi"
    - L35 Fatal error in PMPI_Comm_rank: Invalid communicator / symbol lookup error libopen-pal.so ... undefined symbol __intel_sse2_strlen
    - L41 MPI_Recv: Invalid tag (coupled runs on Cray XC30 / ARCHER)
  - L47 Hangs, aborts, output flushing
    - L48 Job hangs until walltime after a model STOP (e.g. "MOM_IMPLICIT_R: error when solving 3-Diag problem", missing pickup) — `ALL_PROC_DIE`
    - L54 STDOUT.* / STDERR.* not updating until job ends (Cray compiler, ARCHER)
    - L60 No STDOUT.* files; all ranks print to screen in a jumble — `SINGLE_DISK_IO`
  - L66 Stack / memory / size limits
    - L67 Segmentation fault at start or in first big routine (do_oceanic_phys, mom_calc_visc, dynamics) after raising Nr or tile size — `OMP_STACKSIZE`, `MP_STACKSIZE`
    - L73 Link error "relocation truncated to fit: R_X86_64_PC32 / R_X86_64_32S" or "32-bit RIP relative reference out of range" with big sNx*sNy
    - L79 Reduce memory per process
  - L85 Decomposition-dependent results (reproducibility)
    - L86 Results differ when only nPx/nPy/sNx/sNy (or processor count) change
    - L92 Non-hydrostatic run still not decomposition-reproducible with CG2D_SINGLECPU_SUM
    - L98 OBCS sponge layer differs with tile size (strip of non-zero values at inner edge of sponge)
    - L104 Home-made sponge in external_forcing.F repeats at every tile boundary
    - L110 Spurious standing wave / instability with wavenumber nPx, or blow-up on a tile boundary
    - L116 cg2d converges for one processor layout but hits cg2dMaxIters on another (cg2d_res stuck ~1e-6) — `cg2dTargetResidual`, `cg2dMaxIters`
    - L122 cg2dTargetResWunit / cg3dTargetResWunit tolerance scaled wrongly (off by sqrt(Nxy)) — `cg2dTargetResidual`
  - L128 Large processor counts
    - L129 Tile suffix shows ".***.001.data" or "scratch.****", STDOUT.**** (nPx*nSx > 999 or > 9999 procs)
    - L135 EEBOOT_MINIMAL stops when requesting >MAX_NO_PROCS processes (old: 128 -> 1024 -> 2048)
    - L141 MNC: more than 999 tiles / MNC_MAX_ID too small / gluemncbig "found 999 tiles, need N"
    - L147 useMNC unsafe with multi-threads (OpenMP+mnc writes mdsio instead)
  - L153 exch2, cubed-sphere/llc, blank tiles
    - L154 How to use blank (land-only) tiles to save cores
    - L160 S/R W2_SET_MAP_TILES: Domain Total # of tiles = 468 does not match (SIZE.h+blankList)= 492 (also "ABNORMAL END: S/R W2_SET_MAP_TILES")
    - L166 Patchy values / noise at corners of tiles on cubed-sphere with small tiles
    - L172 Blank tiles: output files contain zeros / are not smaller
    - L178 Blank tiles on a lat-lon grid give wrong YC/rA/Coriolis (all tiles in lowest row), or global file too short when last tile is blank
    - L184 w2_e2setup.f: This name does not have a type (nTiles = W2_maxNbTiles); old matlab driver.m topology files
    - L190 forrtl: severe (36): attempt to access non-existent record ... tile001.mitgrid
    - L196 Using one global tile001.mitgrid with multiple tiles on a Cartesian/horizGridFile grid
    - L202 pkg/flt (floats) with exch2/LLC
  - L208 Changing processor count / restart / I/O layout
    - L209 Restart with different nPx/nPy (pickup from other decomposition)
    - L215 Output split in files like U.0000000010.001.002.data (or "output files not separated into 48 subfiles")
    - L221 rdmds('XG.001.001') returns a larger array than the tile (108x30 instead of 54x30)
    - L227 ZONAL_FILT_INIT: Multi-tiles ( nSx*nPx= N ) in Zonal (X) dir. not implemented in Zonal-Filter code
  - L233 Performance and scaling
    - L234 MPI_Allreduce dominates runtime with thousands of tiles (global_sum_tile / cg2d), e.g. llc540 on Pleiades
    - L240 Same Allreduce problem in a package: thousands of GLOBAL_SUM_RL calls (profiles_make_ncfile.F, obsfit) make MPI adjoint 4-10x slower than serial
    - L246 Need a global sum over only some tiles, or vector (element-wise) global sum, in a new package
    - L252 Adjoint/Pleiades slow, Lustre load: many small file open/close in a tight loop (PLEIADES_LUSTRE_OPT1/2/3) — `adTapeDir`, `mdsioLocalDir`, `profilesDir`, `diagMdsDir`, `SINGLE_DISK_IO`
    - L258 BLOCKING_EXCHANGES takes 80% of time after enabling pkg/ptracers (8 tracers)
    - L264 Implicit vertical solver slow (solve_*diagonal.F) — `SOLVE_DIAGONAL_KINNER`
    - L270 Poor scaling, wall-clock >> user time (CG2D, BLOCKING_EXCHANGES dominate); rules of thumb
    - L276 Mac laptop/iMac MITgcm run slow, CPU 2-35%
  - L282 Domain-size and SIZE.h rules
    - L283 Array out of bounds / wrong results with OLx or OLy too small
    - L289 Nx/Ny, tile terminology (what to set for a given core count)
  - L295 OpenMP / threading
    - L296 Threaded (OpenMP) runs: little gain vs MPI, packages break threading
    - L302 Array bound error only in multi-threaded configs: "bi loop upper bound = myByLo" in coupling-interface init
  - L308 Legacy / obsolete
    - L309 mds_byteswapi4.f: Operands of logical operator '.and.' are INTEGER(4)/BOZ (gfortran >=10, OpenMPI) with -DFAST_BYTESWAP
    - L315 Why MPI_COMM_WORLD still appears in eedie.F

## `code_util/claude_code/mitgcm-ecco/references/troubleshooting/physics.md`
- L1 Troubleshooting: physics (free surface, EOS/TEOS, shelfice, RBCS, OBCS/tides, drag/viscosity/Coriolis, advection, cg2d/cg3d)
  - L4 Solver convergence (cg2d / cg3d)
    - L6 cg2dTargetResWunit / cg3dTargetResWunit gives a tolerance that is too loose (off by sqrt(Nxy)) — `cg2dTargetResidual`
    - L12 cg2d converges "too fast"/imprecisely after shifting initial SSH (cg2dNormaliseRHS=T)
    - L18 cg2d: NaN / "STOP ... MOM_IMPLICIT_R (or GAD_IMPLICIT_R): error when solving 3-Diag problem" / "SOLUTION IS HEADING OUT OF BOUNDS ... exceeds allowed range (monSolutionMaxRange=1.E+03)" / "IEEE_INVALID_FLAG" — `advcfl_`
    - L24 cg2d hits cg2dMaxIters and does not converge (flat-bottom/idealised box, tides, rigid lid) then blows up
    - L30 cg2d is the scaling bottleneck (GLOBAL_SUM_TILE_RL dominates, SOLVE_FOR_PRESSURE fraction grows with nProcs) — `ALLOW_SRCG`
    - L36 Non-hydrostatic run 10x slower than hydrostatic; cg3d dominates (SOLVE_FOR_PRESSURE ~ CG3D) — `cg3dMaxIters`
    - L42 Adjoint: tape/recompute problems and wrong cg2d adjoint with nonlinFreeSurf>2 or _RS=real*4 — `AUTODIFF_ALLOW_VISCFACADJ`
  - L48 Free surface, volume and freshwater conservation
    - L50 "CONFIG_CHECK: nonlinFreeSurf not yet implemented in nonHydrostatic code" (and "nonHydrostatic NOT SAFE with non fully implicit barotropic solver") — `selectNHfreeSurf`
    - L56 useRealFreshWaterFlux with implicDiv2DFlow < 1 (or selectAddFluid=1) gives wrong volume
    - L62 addMassFile / ALLOW_ADDFLUID shows no effect, or surface mass source "does nothing"
    - L68 Sea level / mean salinity drifts in a regional or lake domain (SSH decreasing, dynstat_eta_mean drifting)
    - L74 Very negative SSH (-20 m...), "WARNING: hFac < hFacInf", "STOP in CALC_R_STAR : too SMALL rStarFac", "ABNORMAL END: S/R CALC_R_STAR"
    - L80 Surface-cell tracer not conserved / "layers" surface term wrong (linear FS + useRealFreshWaterFlux)
    - L86 Implicit vertical advection with limiter schemes crashes in CALC_R_STAR (e.g. VertAdvScheme=33 + ImplVertAdv=T on llc90)
    - L92 Eta below column depth / dry-looking cells at ebb tide
    - L98 Tides come out damped / weak surface gravity waves with implicitFreeSurface
  - L104 Equation of state, pressure, TEOS-10
    - L106 Which pressure goes into the EOS (JMD95Z vs MDJWF vs TEOS10): selectP_inEOS_Zc — `selectP_inEOS_Zc`
    - L112 eosType='TEOS10': FIND_ALPHA/FIND_BETA derivative errors; surface heat flux bias from theta=CT
    - L118 Atmospheric loading (pLoad / apressure) enters the EOS as full pressure — `surf_pRef`
    - L124 NaN/blow-up with extreme T/S initial conditions (hypersaline lake, paleo 130 psu, huge T/S range)
    - L130 Linear EOS: RHOAnoma zero / weird, "dRho", tRef semantics
  - L136 Ice shelves (pkg/shelfice, frazil, GM/Redi, mixing)
    - L138 shelfice_thermodynamics.F NaN: freshWaterFlux = shiTransCoeffS*(1 - sLoc/saltFreeze) with saltFreeze = 0
    - L144 Supercooled water (theta ~ -100 degC) in cavity cells, sea ice growing next to cavity
    - L150 Ice-shelf masking wrong where draft is thinner than top cell / Ro_surf ~ 1e-13 (kTopC = kSurfC in open ocean); remeshing cannot merge thin ice — `SHI_update_kTopC`
    - L156 pkg/gmredi under ice shelves loses/gains tracer; 'linear' GM taper extrapolates through ice; KPP/GGL90 under shelves
    - L162 Turn off / prescribe ice-shelf thermodynamics; remove conduction; prescribe melt rate
    - L168 Ice-shelf cavity blows up at start / pressure at ice base (SHELFICEloadAnomalyFile, rhoConst/rhoNil)
    - L174 SHI_withBL_realFWFlux / SHELFICEboundaryLayer with autodiff -> STOP
    - L180 pkg/frazil: not coupled to shelfice; freezing point differs from shelfice/seaice
    - L186 Ice-shelf diagnostics confusion (SHIForcS vs SHIfwFlx, tendencies not where flux applied)
  - L192 RBCS (restoring) checklist
    - L194 RBCS appears to do nothing / "pkg/rbcs was not compiled (ALLOW_RBCS undef)" / namelist error (exit 19) — `maskLEN`
    - L200 RBCS with tauRelax < deltaT: noise or no restoring
    - L206 RBCS time axis: forcing offset, single-time files, cycles, error "startTime before rbcsForcingOffset+0.5*rbcsForcingPeriod"
    - L212 RBCS velocities blow up / asymmetry when fields are on the wrong grid
    - L218 Model much slower with RBCS (3-D relax fields)
  - L224 OBCS, tides and sponges
    - L226 SSH does not follow the parent model / drifts after enabling tides (useOBCSbalance, useOBCStides)
    - L232 Ways to force tides; units; division by zero in obcs_add_tides — `tidePot`
    - L237 Spurious boundary jets / noise at open boundaries (CD scheme, Orlanski, prescribed-only, sponge) — `spongeThickness`
    - L243 Model blows up / strange fresh river inflow through OBCS (negative salinity, NaN in cg2d)
    - L249 Unexpected periodicity: same tide at E/N and W/S, wrap-around (MITgcm is periodic by default)
  - L255 Drag, viscosity, Coriolis
    - L257 Bottom drag confusion: no_slip_bottom vs bottomDragLinear/Quadratic vs sideDragFactor; too-strong tidal channel currents
    - L263 Grid-scale noise / waves near boundaries on a lat-lon beta-plane (tutorial_baroclinic_gyre)
    - L269 Viscosity/biharmonic limits and cosPower on isotropic grids — `ISOTROPIC_COS_SCALING`
    - L275 vectorInvariantMomentum + highOrderVorticity reverses flow over topography (ACC channel)
    - L281 Momentum budget will not close (TOTUTEND vs terms) — `ALLOW_MOM_TEND_EXTRA_DIAGS`
    - L287 Pressure diagnostics: PHIHYD, PHL, PH, PHRefC, z* differences — `PHIHYDcR`
  - L293 Advection schemes and noise
    - L295 Instability/very different circulation with DST schemes (30, 33, 77) unless staggerTimeStep
    - L301 Negative salinity / overshoots (T=49 degC, S=57 at bottom) with centered advection — `GAD_ALLOW_TS_SOM_ADV`
    - L307 Grid-scale noise in non-hydrostatic wake/internal wave runs
    - L313 Non-hydrostatic + KPP warnings; 3-D Smagorinsky coefficients
  - L319 Sea ice dynamics that interact with ocean physics (brief)
    - L321 Stripes/melt at tile (CPU) boundaries and immobile ice at 20 m grids with LSR
    - L327 EVP solver unstable in near ice-free cells (denomU/V ~ 0) — `SEAICE_evpAreaReg`
  - L333 Build / run gotchas touching physics
    - L335 Segfault in adjoint or big runs: "Caught signal 11 (Segmentation fault: address not mapped...)"
    - L341 Output "zero every 4th level" / reading bathy gives zeros / MON stats all zero except time and vorticity

## `code_util/claude_code/mitgcm-ecco/references/troubleshooting/restart.md`
- L1 Troubleshooting: pickups, restarts, nIter0, checkpoints
  - L5 Choosing / naming the pickup (nIter0, pickupSuff, ckptA/B)
    - L6 Restart fails: "MDS_READ_FIELD: File does not exist" / reads wrong pickup / only pickup.ckptA exists — `timeStepNumber`, `nIter0`, `pChkptFreq`, `startTime`
    - L12 Restarted run deviates from the continuous run right after restart (forcing/phase wrong, SST jump) — `baseTime`, `startTime`
    - L18 Restart with a different deltaT (pickup named for old iteration, time axis wrong, EXF reads wrong forcing)
    - L24 Pickup iteration after switching calendar (360-day "model" vs 365/Gregorian) or calendarDumps — `repeatPeriod`
  - L30 Exact (bit-for-bit) restart; results differ after restart
    - L31 Restarted run differs from non-restart run (T,S,U diverge after the first step) — `GLOBAL_SUM_TILE`
    - L37 Results change when the number/size of tiles changes between pickups
  - L43 NetCDF (MNC) pickups
    - L44 WARNING MNC_READPARMS: incomplete MNC pickup files implementation / MNC_CW_RL_R: cannot open either a per-face or a per-tile file pickup...t001.nc / variable 'RC' already defined in file 'grid.t001.nc'
  - L50 Changing decomposition, resolution or initial fields
    - L51 Restart on a different number of MPI processes (or different tile size)
    - L57 Start a new/finer run from a coarser run's state, or replace theta/salt/eta/u/v at restart — `hydrogThetaFile`, `hydrogSaltFile`, `uVelInitFile`, `vVelInitFile`, `pSurfInitFile`, `HeffFile`, `AreaFile`, `HsnowFile`, `HsaltFile`, `UiceFile`, `ViceFile`
    - L63 Restart with changed physics: seaice initial values ignored; nonhydrostatic from hydrostatic pickup; EOS change — `phi_nh`
  - L69 Old/mismatched pickups, .meta handling
    - L70 Reading an old pickup with new code (garbled fields, CALC_R_STAR at first step, salt looks like temperature)
    - L76 Pickup read/write fixed upstream: missing-field list index typo, thsice MNC pickup, cheapaml passive tracer — `useCheapTracer`
    - L82 New field added to the pickup (e.g. QH stagger tendency, NH) breaks old restart — `ALLOW_QHYD_STAGGER_TS`
  - L88 Reading, editing and writing pickups by hand
    - L89 Edit a pickup in Matlab/Python (which record is which field; "floating invalid" after writing) — `fldList`, `readBinaryPrec`
    - L95 Nonlinear free surface restart: do not shift eta; Eta fields consistent
  - L101 Forcing / boundary conditions at restart
    - L102 Restart and forcing files: do I need to change OBCS/forcing file names or start records? — `externForcingPeriod`
    - L108 Starting tracers / packages from scratch on top of a physics pickup (ptracers, DIC, floats) — `PTRACERS_initialFile`, `PTRACERS_EvPrRn`, `pickupSuff`, `diag_pickup_read`
    - L114 Orlanski open boundary: W or boundary phase speed zero after restart
  - L120 Adjoint / control restarts
    - L121 Adjoint (or forward) restart with time-varying gentim2d controls: "Non-existing record number" in xx_*.effective.*.data / "attempting to access non-existing record in xx_precip"
    - L127 Divided adjoint (restart of adjoint) and pkg/profiles segfault
  - L133 Platform / scripting
    - L134 Segmentation fault when writing pickup (Darwin/DIC, big runs)
    - L140 Segfault on restart because of NaN in bathymetry/other input (not the pickup)
    - L146 Automating restarts under PBS/SLURM (find most recent pickup, update nIter0, resubmit) — `nIter0`, `most_recent_pickup.sh`
    - L152 Vector-computer / aggressive optimisation: crash in seaice_growth after unrelated change (stale restart behaviour)
    - L158 Misc setup errors seen alongside restarts

## `code_util/claude_code/mitgcm-ecco/references/troubleshooting/seaice.md`
- L1 Troubleshooting: pkg/seaice and thsice
  - L4 Blow-ups, NaNs, runaway ice thickness
    - L6 "STOP in CALC_R_STAR" / r* error / surface cell dries, with HEFF of 10-50 m in a few cells — `SEAICE_drag_south`, `SEAICE_drag`, `SEAICE_CAP_ICELOAD`
    - L13 Blow-up (MOM_IMPLICIT_R / extreme potential temperature) only with EVP; LSR runs fine — `SEAICEuseEVPrev`, `SEAICEuseEVP`, `SEAICEnEVPstarSteps`, `SEAICEaEVPcoeff`, `monitorFreq`, `SEAICE_evpAreaReg`, `lab_sea.hb87`
    - L19 NaNs / blow-up at high resolution (O(1 km)) in summer melt with LSR; smaller deltaT helps only temporarily — `LSR_ERROR`, `SEAICEnonLinIterMax`
    - L25 Stripes / "scars" in ice area, shear, divergence along MPI tile edges — `SEAICEnonLinIterMax`
    - L31 NaN in first steps with SEAICEuseDYNAMICS=.TRUE. and no ice present (old code) — `SEAICEdWatMin`, `SEAICE_CHECK`
    - L37 Sudden crash in seaice_growth ("Invalid operation PROG=seaice_growth") only with ALLOW_SITRACER on a vector machine
    - L43 Ice thickness jumps 50 m in one step next to a freshwater (EmP) patch
    - L49 Super-cooled water from ice-shelf cavities makes 40 m of ice / "extreme potential temperature" at cavity mouth — `SEAICE_gamma_t_frz`, `SEAICE_mcPheeTaper`, `SEAICE_frazilFrac`, `SEAICE_mcPheePiston`
    - L55 "SEAICE_CHECK: SEAICE_mcPheePiston is out of bounds ... must lie within 0. and drF(1)/SEAICE_deltaTtherm" — `SEAICE_mcPheePiston`, `SEAICE_availHeatFrac`
    - L61 Ocean+ice model blows up at high wind (40 m/s) at the seafloor next to the ice edge
    - L67 Water below the freezing point at depth in coarse ice-ocean runs (-2.4 C) — `useOldFreezing`
  - L73 Dynamics, solver and defaults
    - L75 Ice does not move at all / starts late at very high resolution (dx of metres to 100 m), LSR, free drift — `LSR_minIniGuess`, `LSR_mixIniGuess`
    - L81 Floating-point exception (div by zero) in seaice_lsr.F with SEAICEscaleSurfStress=T and LSR_mixIniGuess=1
    - L87 Results changed after updating: ice thicker / different drag / snow mass — `SEAICEscaleSurfStress`, `SEAICEaddSnowMass`, `SEAICE_useMultDimSnow`, `SEAICE_drag`, `SEAICE_waterDrag`, `global_ocean.cs32x15.icedyn`
    - L93 pkg/thsice + SEAICE: snow mass missing in ice load (snowHeight not mapped to HSNOW) — `snowHeight`
    - L99 Metric terms in sea-ice stress divergence changed results near the pole / on curvilinear grids — `SEAICEuseMetricTerms`, `SEAICEselectMetricTerms`, `SEAICEuseExtraMetricTerms`
    - L105 Compile error "Expecting END DO statement" in seaice_solve4temp.F
    - L111 Compile error in seaice_evp.F (multiple definition of nEVPstepMax) when ALLOW_AUTODIFF_TAMC undefined; STOP with SEAICE_deltaTdyn > SEAICE_deltaTtherm + autodiff
    - L117 "SEAICE_CHECK"/readparms error: deltaTdyn vs deltaTtherm consistency — `SEAICE_deltaTdyn`
    - L123 Ice seems to ignore SEAICE_waterDrag / drag has no effect in idealised runs — `SEAICEdWatMin`
    - L129 Ice velocity locked to a constant in a doubly periodic / uniform-forcing domain — `SEAICE_strength`, `SEAICE_waterDrag`
    - L135 Sea ice on cubed-sphere corners / Cartesian grid / periodic channel blows up
  - L141 Package interaction (OBCS, shelfice, thsice, exf)
    - L143 Ice piles up or forms spurious velocity at open boundaries (SEAICE + OBCS); how to relax ice at boundaries — `OBCS_SEAICE_AVOID_CONVERGENCE`, `OBCS_SEAICE_COMPUTE_UVICE`, `ALLOW_OBCS_SEAICE_SPONGE`, `seaice_obcs.sponge`
    - L149 Periodic seaice but OBCS for the ocean
    - L155 Thin ice-shelf (hFacC(k=1)>0) with seaice or thsice fails / shortwave warms ocean under shelf
    - L161 Large Eta / sea-ice loss next to a shelfice cavity (EVP/LSR regional, SHELFICE + SEAICE) — `SHELFICEloadAnomalyFile`
    - L167 seaice + thsice together: initial ice thickness set in data.seaice has no effect — `thSIceThick_InitFile`, `THSICE_REGULARIZE_CALC_THICKN`, `THICKNESS_AUTODIFF_REGULARIZE`
    - L173 pkg/seaice cannot be used in coupled atmosphere (cpl_aim+ocn); thsice instead — `ALLOW_ATM_TEMP`, `ALLOW_DOWNWARD_RADIATION`, `ALLOW_BULKFORMULAE`, `EXF_READ_EVAP`
    - L179 NLFS=4 (nonlinear free surface) with plain z-coordinates + seaice: allowed? (surface cell can dry)
    - L185 Imposing observed sea-ice (concentration from file) instead of computing it — `exf_iceFraction`, `SEAICE_ALLOW_AREA_RELAXATION`
    - L191 Ice drift under a stationary ocean / stand-alone sea ice
  - L197 Initialisation, restarts, parameter handling
    - L199 Initial ice area from HeffFile ignored / restart ignores SEAICE_initialHEFF — `HeffFile`, `AreaFile`
    - L205 Typo in data.seaice silently ignores everything after it
    - L211 Old SEAICE_OPTIONS.h / darwin code fails to compile after update ("COMMON block data object must not be an automatic", SEAICE.h)
    - L217 Interpolating pickup_seaice to a coarser grid
  - L222 Diagnostics and budgets
    - L224 Which diagnostic is stress on ocean under ice; heat/FW budgets that do not close — `oceQnet`
    - L230 SEAICE_USE_GROWTH_ADX: heat/mass budget not closed, +25 C spike when snow floods
    - L236 Salt/FW from sea ice: SEAICE_saltFrac and seaice_salt0 do nothing — `SEAICE_saltFrac`, `SEAICE_VARIABLE_SALINITY`, `SEAICE_salt0`, `ALLOW_SITRACER`
    - L241 Track or zero the freshwater from ice melt

## `code_util/claude_code/mitgcm-ecco/references/troubleshooting/tracers.md`
- L1 Troubleshooting: tracers, BGC, mixing and viscosity
  - L4 ptracers set-up, restart, namelists
    - L6 "invalid reference to variable in NAMELIST input" / "namelist not terminated with / or &end" / "End of file" in ptracers_readparms or any data.* read
    - L12 PTRACERS_Iter0 > nIter0 stops run ("wrong setting of PTRACERS_Iter0 ... nIter0 < PTRACERS_Iter0 < nEndIter not supported"), or tracers are all zero
    - L18 Restart with pickupSuff='ckptA' or after enabling bgc: DIC/ptracers NaN or zeros at domain edges
    - L24 All passive tracers behave like salinity / tracer concentrations track S
    - L30 How to add a source/sink, surface flux, runoff tracer, melt-water tracer — `ptracers_forcing_surf.F`, `ptracers_apply_forcing.F`, `gchem_calc_tendency.F`
    - L36 Turn off advection (or change vertical scheme) for a single passive tracer
    - L42 RBCS relaxation of tracers: "Non-existing record number", namelist error, tracer never reaches target
    - L48 Open boundaries and ptracers: tracer explodes (1e38, negative), accumulates at boundary, or zeros leak in
    - L54 Negative tracer concentrations after advection (undershoot)
    - L60 SIGFPE / blow-up in gad_fluxlimit_* (Rj ~ 1e-307, "Rjp,Rj" printed)
    - L66 Tracer content drifts / surface-cell budget does not close (linear free surface)
    - L72 Closing a tracer budget from diagnostics (ADV*, DF*, KPPg_*, TOTTTEND) does not balance
    - L78 Which flux diagnostics hold which tracer-mixing contribution (DFrE vs DFrI, GM, KPP)
    - L84 Slow runs with many ptracers (BLOCKING_EXCHANGES dominates), 51 tracers 13x slowdown
    - L90 Segmentation fault at start/pickup write/first step with more levels or Darwin
    - L96 Offline tracer runs: set-up and gotchas — `offlineOffsetIter`, `deltaToffline`
    - L101 deltaTtracer / deltaTmom / deltaTClock: raising deltaTtracer alone gives no speed-up, ice does not form, tendencies tiny
  - L107 gchem / DIC / BLING / Darwin
    - L109 GCHEM_CHECK: "Number of GCHEM tracers gchem_Tracer_num = N exceeds number of pTr: PTRACERS_numInUse" — `ALLOW_O2`
    - L115 "useDIC and useCFC cannot both be .TRUE." / running two bgc packages together
    - L121 BLING with `#define USE_SIBLING`: SIGSEGV in bling_bio_nitrogen_ (G_SI)
    - L127 BGC forcing files (wind, ice, atmospheric pCO2, silica, iron, PAR) not read / wrong units / "Cannot match namelist object name apco2startdate1" / "attempt to access non-existent record"
    - L133 DIC / BLING NaN in pH solver, DIC -> NaN after restart or at river runoff points — `DIC_NO_NEG`, `CARBONCHEM_SOLVESAPHE`, `selectPHsolver`
    - L139 Carbonate system bugs seen in old checkouts (check your checkpoint) — `DIC_CALCITE_SAT`
    - L144 pkg/dic forcing/initial-field loading bugs
    - L150 Closing a DIC budget: DICBIOA/DICCARB don't match tendency
    - L156 Fresh-water dilution of DIC/ALK/PO4 and PTRACERS_ref confusion — `ALLOW_OLD_VIRTUALFLUX`
    - L162 Darwin3 (ecco_darwin) compile errors "DARWIN_partScav / DARWIN_minFeLoss / ironSedFlux has no type"; radtrans
    - L168 Delaying biogeochemistry until physics is spun up
  - L174 GMRedi
    - L176 Only Redi (no GM), or only GM; per-tracer GM
    - L182 3-D (or 1-D x 2-D) GM/Redi diffusivity from files — `GM_READ_K3D_REDI`, `GM_READ_K3D_GM`, `GM_isopycK3dFile`, `GM_background_K3dFile`
    - L188 GMRedi under ice shelves / dry cells at top (shelfice): strong spurious mixing, tracer non-conservation, blow-up with advective GM
    - L194 GMREDI_WITH_STABLE_ADJOINT / GM_taper_scheme='stableGmAdjTap' (ECCO v4r5): wrong TAF keys in gmredi_slope_limit.F
    - L200 Results changed in 2026: GMRedi Bates/GEOMETRIC/QGLeith fixes (PRs #1024, #1025 from issue #1012) — `ALLOW_LEITH_QG`
  - L206 KPP
    - L208 "*** ERROR *** Some form of convection has been enabled" (KPP_CHECK)
    - L214 KPP stops/needs SHORTWAVE_HEATING; "KPP cannot work in the dark"
    - L220 KPP checkerboard / grid-scale noise in KPPviscAz/KPPdiffKzT, vertical banding at fine dz — `KPP_SMOOTH_SHSQ`, `KPP_SMOOTH_DBLOC`, `KPP_SMOOTH_DVSQ`, `KPP_SMOOTH_DENS`, `ALLOW_KPP_VERTICALLY_SMOOTH`, `KPP_DO_NOT_MATCH_DIFFUSIVITIES`, `KPP_DO_NOT_MATCH_DERIVATIVES`
    - L226 KPP_ESTIMATE_UREF does not compile; surface-velocity reference
    - L232 KPP diffusivity maxima at the bottom (shear mixing), KPPhbl stuck at minKPPhbl or 5 m, KPP in unforced runs
    - L238 KPP and ice shelves / non-hydrostatic / inverse stratification
  - L244 GGL90, other vertical mixing, convection
    - L246 GGL90: results changed / bugs fixed (check your checkpoint) — `GGL90_MISSING_HFAC_BUG`, `GGL90_IDEMIX_CVMIX_VERSION`, `ALLOW_GGL90_LANGMUIR`, `ALLOW_GGL90_HORIZDIFF`
    - L251 GGL90 too much deep mixing (Antarctic) / grid-scale velocity patterns on llc90; mixing-efficiency default
    - L257 Weak/no effect of RiLimit (pp81) / RiMax (my82); MY82 stops; MY82 units
    - L263 Convective adjustment: cAdjFreq vs ivdc_kappa; AD bug in convective_adjustment.F
    - L269 KL10 spurious huge viscosity with non-linear EOS
    - L275 "CONFIG_CHECK: diffKrFile is set but never used" / NaN with a 3-D diffKr file / changing diffKrT does nothing — `ALLOW_3D_DIFFKR`, `diffKrFile`
  - L281 Horizontal viscosity, noise, stability
    - L283 Choosing viscAh / viscA4 and viscAhGrid; "SOLUTION IS HEADING OUT OF BOUNDS" with viscosity CFL — `viscAhGrid`, `viscA4Grid`
    - L289 Striping/noise along basin boundaries on spherical grids (tutorial_baroclinic_gyre noisy)
    - L295 Biharmonic viscosity: ISOTROPIC_COS_SCALING/cosPower, L4rdt change, noise at 1/4-1/6 deg
    - L301 Leith / Smagorinsky usage pitfalls — `ALLOW_LEITH_QG`, `viscC2LeithQG`, `useSmag3D`, `smag3D_coeff`, `smag3D_diffCoeff`
    - L307 Dissipation diagnostics: Um_Diss etc. don't balance, units m^4/s^2
    - L313 Time-dependent or state-dependent viscosity/diffusivity (custom code)
    - L319 Spurious mixing / variance destruction with flux-limited schemes; choosing advection for tracers
    - L325 Mixed layer / vertical-mixing background choices (what to set under KPP)
- `code_util/transport/readme.txt` (1 lines) — ECCO-Darwin transport
- `idealized/1D_darwin/sea_ice_column/readme.txt` (27 lines) — v06 1D_ocean_ice_column_darwin idealized experiment — `darwin3-darwin`

## `idealized/1D_darwin/time_series_sites/README.md`
- L1 1-D ECCO-Darwin v05 columns at ocean time-series sites — `HydroS`
  - L31 Quick start: run the optimized set-ups — `opt_manifest.txt`
  - L55 Set-up — `delR`, `darwin_ckpt68g`, `xx_diffkr`, `make_runs.py`, `PTRACERS_EvPrRn`, `endTime`
    - L79 Grazer + nitrogen-fixer closure (optimized set-ups only)
  - L94 Directory contents — `code_darwin`, `code_v5r1`, `SIZE.h`, `RBCS_SIZE.h`
    - L106 Site input files
  - L124 1. Get code
  - L135 2. Build executable
  - L150 3. Make v05 baseline run directories and run
  - L164 4. (Re)extract site inputs from v5r1 (NASA Pleiades only)
  - L173 5. Green's-function optimization — `KZ_BG`, `diffkr_1x1x50`, `PCMAX_Syn`, `PCMAX_diazo`
  - L229 6. Results — `make_optimized_runs.py`, `code_rbcs_mask21`
  - L278 7. Observations
  - L286 Contact

## `idealized/1D_darwin/time_series_sites/obs/README_obs.md`
- L1 Observations for 1-D ECCO-Darwin columns (README_obs)
  - L11 Layout — `std_hot.py`, `std_bios.py`, `std_papa.py`, `std_pap.py`, `std_synth.py`, `std_tm.py`, `std_geotraces.py`, `std_sots.py`, `std_k2.py`, `std_cariaco.py`, `std_dyfamed.py`, `std_munida.py`, `std_owsm.py`, `fetch_synth.sh`, `subset_glodap.py`, `summarize_core.py`, `obs_summary.csv`, `obs_summary_core.csv`
  - L31 Tidy file columns — `depth_m`, `SiO2`, `Chl_HPLC`, `Chl_sat`, `pCO2`, `fCO2`, `pCO2_air`, `POC_flux`, `PON_flux`, `POP_flux`, `PIC_flux`, `bSi_flux`, `mass_flux`, `lith_flux`, `BactAbund`, `ProAbund`, `O2_CTD`, `NO3_LL`, `PO4_LL`, `SRP_LL`, `d13C_DIC`, `PP_dark`
  - L46 Flag convention (summary; details per source below)
  - L56 Units and density conversion
  - L64 General caveats
  - L76 Summary (site x variable, all sources combined; full breakdown in obs_summary.csv)
  - L277 Missing data and sources needing login or manual fetch
- L298 Per-source notes
  - L302 HOT / Station ALOHA (22.75N, 158.0W) - sources — `std_tm.py`
    - L306 Access problems (read first; re-tried 2026-10-08) — `bcodmo_dataset_3773_v3`
    - L320 1. raw/bcodmo_bottle/ - HOT Niskin bottle + CTD-at-bottle data — `hot_bottle_bcodmo3773_stn2.csv`, `hot_bottle_bcodmo3773_stn2_part2_YYYY.csv`
    - L342 2. raw/bcodmo_pp/ - 14C primary production — `hot_pp_bcodmo737163.csv`
    - L356 3. raw/bcodmo_flux/ - sediment traps (floating, 150 m standard; also 70-500 m early) — `hot_flux_bcodmo737393.csv`
    - L367 4. raw/hot_co2/ - HOT surface CO2 system product (Dore)
    - L377 5. raw/ncei_whots_pco2/ - WHOTS mooring MAPCO2 (PMEL), NCEI OCADS accession 0100080
    - L392 6. raw/scripps_co2/ - Scripps (Keeling) surface DIC/ALK at ALOHA
    - L401 7. raw/bcodmo_pp_v5/ - 14C primary production, BCO-DMO 737163 v5 (HOT-1 to HOT-355), added 2026-10-08 — `landing_page_737163.html`, `bcodmo_pp`, `bcodmo_pp_v5`, `Start_ISO_DateTime_UTC`, `Flag_Light`, `Flag_Bottle`, `Flag_Dark`
    - L418 8. raw/bcodmo_flux_v4/ - sediment traps, BCO-DMO 737393 v4 (HOT-2 to HOT-355), added 2026-10-08 — `bcodmo_flux`
    - L430 9. raw/bcodmo_spots/ - SPOTS synthesis, ALOHA bottle subset (fills HOT bottle chemistry to Dec 2019), added 2026-10-08 — `query_url.txt`, `bcodmo_dataset_896862_v2.das`, `spots_v2_ALOHA_2016on.csv`, `bcodmo_bottle`
    - L456 10. raw/bcodmo_hot_metals/ - dissolved trace metals at ALOHA, BCO-DMO 962986 v2 (Dec 2020 to Nov 2023), added 2026-10-08 — `bcodmo_dataset_962986_v2.csv`, `std_tm.py`
    - L470 Not obtained
  - L487 BATS (Bermuda Atlantic Time-series Study) - sources
    - L502 1. bcodmo_bottle - BATS discrete bottle file (v10)
    - L520 2. bcodmo_pp - BATS 14C primary production (v8)
    - L529 3. bcodmo_pigments - HPLC + fluorometric pigments (v10)
    - L539 4. bcodmo_flux - BATS PITS sediment traps (v8)
    - L550 5. ncei_ofp - Oceanic Flux Program deep moored traps (v2, 1978-2018)
    - L567 6. ncei_btm_mooring - Bermuda Testbed Mooring MAPCO2 (2005-2007)
    - L577 7. scripps_co2 - Scripps CO2 program surface DIC/ALK at BATS (BATS.csv)
    - L586 8. bcodmo_bait_fe - BAIT 2019 dissolved Fe + d56Fe (Conway lab), BCO-DMO 936824 v1, added 2026-10-08 — `std_tm.py`
    - L596 9. bcodmo_bait_tm_bottle - BAIT 2019 trace-metal rosette bottle data (Sedwick lab), BCO-DMO 937302 v1, added 2026-10-08
    - L608 Not fetched / notes
  - L617 Hydrostation S (Bermuda) - sources
    - L622 1. bcodmo_bottle - Hydrostation S discrete bottle file (v8)
    - L636 2. scripps_co2 - Scripps CO2 program surface DIC/ALK at Hydrostation S (BERM.csv)
    - L648 Not fetched
  - L654 Ocean Station Papa / Line P P26 (model column 50.1N, 144.9W): source notes
    - L658 1. DFO-IOS bottle profiles (Line P / Station P), `raw/ios_bot_erddap/` — `IOS_BOT_Profiles`, `IOS_BOT_Profiles_P26.csv`, `IOS_BOT_Profiles_metadata.csv`
    - L679 2. Line P carbonate system (NCEI OCADS), `raw/ncei_linep_carbon/`
    - L698 3. NOAA PMEL Papa CO2/OA surface mooring, `raw/pmel_co2_mooring/` — `papa_co2_mooring.csv`, `papa_co2_mooring_metadata.csv`
    - L712 4. NOAA PMEL OCS Papa daily gridded T/S (OceanSITES), `raw/oceansites_pmel_daily/` — `OS_PAPA_200706_M_TSVM_50N145W_dy.nc`
    - L725 5. OOI Global Station Papa flanking moorings A/B (riser BGC), `raw/ooi_papa_flanking/`
    - L740 6. Ocean Station P sediment traps (Wong et al. 1999), `raw/pangaea_osp_trap/` — `PANGAEA_92552.tab`
    - L752 Caveats
  - L763 PAP (Porcupine Abyssal Plain Sustained Observatory) – sources — `standardize_obs.py`
    - L769 1. OceanSITES PAP-1 / PAP-2 / PAP-3 mooring NetCDF (primary) — `depth_m`
    - L808 2. NCEI OCADS PAP-SO mooring pCO2 (accession 0312034) — `sbo_37_ox`
    - L821 3. NCEI OCADS RRS James Cook JC278 underway pCO2 (accession 0312066, EXPOCODE 740H20250530)
    - L830 4. PANGAEA – PAP sediment traps (Lampitt et al.)
    - L850 5. PANGAEA – JC087 bottle nutrients/O2 at PAP (Jun 2013)
    - L858 6. BODC NODB – PAP cruise bottle data, 2013-2019 (source tag `bodc_nodb_<cruise>`; added 2026-10-08)
    - L874 Not obtained / needs login
    - L906 Caveats for model comparison
  - L914 Synthesis-product subsets (all sites) — `standardize_obs.py`
    - L927 GLODAPv2.2023 (source tag `glodap`) — `GLODAPv2.2023_Atlantic_Ocean.csv.zip`, `GLODAPv2.2023_Pacific_Ocean.csv.zip`
    - L945 SOCAT v2026 (source tag `socat`)
    - L962 BGC-Argo synthetic profiles (source tag `bgcargo`)
    - L981 Satellite chlorophyll: ESA OC-CCI v6.0 monthly (source tag `satchl`) — `query_url.txt`
    - L994 GEOTRACES dissolved Fe: IDP2025 NOT downloaded (re-tried 2026-10-08); open BCO-DMO substitutes used — `std_tm.py`, `bcodmo_hot_metals`, `bcodmo_bait_fe`, `bcodmo_bait_tm_bottle`
    - L1015 Tidy-row counts from these products (after QC), 2026-10-07
    - L1027 Update 2026-10-09: four new sites (SOTS, K2, CARIACO, DYFAMED) — `fetch_synth.sh`, `subset_glodap.py`, `std_synth.py`, `std_geotraces.py`
    - L1041 Tidy-row counts from these products (after QC), 2026-10-09
    - L1051 Update 2026-10-09 (b): two more sites (MUNIDA, OWSM) — `fetch_synth.sh`, `subset_glodap.py`, `std_synth.py`, `std_geotraces.py`, `standardize_obs.py`, `std_munida`, `std_owsm`, `std_synth`, `std_geotraces`, `ncei_ocads_owsm`, `nmdc_noremso_mooring`, `ncei_ocads_0100218`, `ncei_ocads_0112883`
    - L1078 Tidy-row counts from these products (after QC), 2026-10-09 (b)
  - L1086 GEOTRACES IDP2025 dissolved Fe (added 2026-10-08)
  - L1101 G-POM suspended POC/PON/POP (added 2026-10-18) — `bin_obs`

## `idealized/1D_darwin/time_series_sites/obs/SOURCES_synthesis.md`
- L1 Synthesis-product subsets (all sites) — `standardize_obs.py`
  - L14 GLODAPv2.2023 (source tag `glodap`) — `GLODAPv2.2023_Atlantic_Ocean.csv.zip`, `GLODAPv2.2023_Pacific_Ocean.csv.zip`
  - L32 SOCAT v2026 (source tag `socat`)
  - L49 BGC-Argo synthetic profiles (source tag `bgcargo`)
  - L68 Satellite chlorophyll: ESA OC-CCI v6.0 monthly (source tag `satchl`) — `query_url.txt`
  - L81 GEOTRACES dissolved Fe: NOT downloaded
  - L92 Tidy-row counts from these products (after QC), 2026-10-07
  - L103 Update 2026-10-08
  - L105 GEOTRACES dissolved Fe: IDP2025 NOT downloaded (re-tried 2026-10-08); open BCO-DMO substitutes used — `std_tm.py`, `bcodmo_hot_metals`, `bcodmo_bait_fe`, `bcodmo_bait_tm_bottle`
  - L126 Update 2026-10-09: four new sites (SOTS, K2, CARIACO, DYFAMED) — `fetch_synth.sh`, `subset_glodap.py`, `std_synth.py`, `std_geotraces.py`
  - L140 Tidy-row counts from these products (after QC), 2026-10-09
  - L150 Update 2026-10-09 (b): two more sites (MUNIDA, OWSM) — `fetch_synth.sh`, `subset_glodap.py`, `std_synth.py`, `std_geotraces.py`, `standardize_obs.py`, `std_munida`, `std_owsm`, `std_synth`, `std_geotraces`, `ncei_ocads_owsm`, `nmdc_noremso_mooring`, `ncei_ocads_0100218`, `ncei_ocads_0112883`
  - L177 Tidy-row counts from these products (after QC), 2026-10-09 (b)
- `offline/V4r6_darwin_offline/readme.txt` (394 lines) — Offline ECCO-Darwin: LLC90 biogeochemistry driven by pre-computed ECCOv4r6 physics — `checkpoint68g`, `darwin3`, `darwin3/pkg/darwin`, `darwin3/tools/genmake2`, `darwin3/tools/darwin/cogapp/cogapp.py`, `linux_amd64_ifort+mpi_ice_nas.electra_skylake_20251013`, `linux_amd64_ifort+mpi_ice_nas`, `V4r6`, `LLC90`, `LLC270`, `llc270`, `git clone https://github.com/MITgcm/MITgcm.git -b`, `git clone https://github.com/MITgcm-contrib/ecco_darwin.git cd`, `git clone https://github.com/ECCO-GROUP/ECCO-v4-Configurations.git mv`

## `pkg/bbl/README.md`
- L1 pkg/bbl: an embedded bottom boundary layer for MITgcm and ECCO-Darwin
  - L17 Contents
  - L39 What the hooks do — `BBL_TENDENCY_APPLY_PTR`, `BBL_INIT_VARIA`, `PTRACERS_INIT_VARIA`, `CTRL_INIT_VARIABLES`, `BBL.h`, `BBL_PTR.h`, `the_main_loop.F`, `forward_step.F`
  - L50 Install in darwin3 (ECCO-Darwin)
  - L66 Install in MITgcm — `global_oce_latlon`
  - L81 Using it — `pickup_bbl`, `bbl_doPtracers`, `BBLuTrp`, `BBLvTrp`, `BBLwEnt`, `BBL_OPTIONS.h`
  - L106 Tests — `global_oce_latlon`
  - L123 Notes — `maskInW`, `maskInS`, `BBL_CHECK_BUDGET`

## `pkg/bbl/src/README_pkg_bbl.md`

## `pkg/wad/README.md`
- L1 pkg/wad — wetting and drying for MITgcm — `addMass`
  - L15 What changed in this update (2026-10-06) — `useRealFreshWaterFlux`, `EXTERNAL_FORCING_SURF`, `WAD_DRY_FORCING`, `PmEpR`, `wad_iceground`, `SALT_PLUME_VOLUME`, `wad_iceground_sp`
  - L38 Contents
  - L56 Install — `WAD_PARM01`
  - L85 Required model settings — `wad_check`, `SALT_PLUME_VOLUME`, `SEAICE_ITD`, `SEAICE_USE_GROWTH_ADX`, `PTRACERS_addSrelax2EmP`, `wadDryForcing`, `vectorInvariantMomentum`, `upwindShear`
  - L115 Parameters (`data.wad`, namelist `WAD_PARM01`) — `wadMinDepth`, `wadCritDepth`, `wadAdvDepth`, `wadGMDepth`, `wadKPPCapDepth`, `wadKPPDepth`, `wadDryForcing`, `wadIceDepth`, `wadIceGround`, `wadIceTopMelt`, `wadLoadDepth`, `wadLoadTau`, `wadIceSeal`, `wadIceMinHeff`, `wadIceOBpack`, `wadUpwindFace`, `wadCarryVel`, `momImplVertAdv`, `wadDragDepth`, `wadManningN`, `ALLOW_BOTTOMDRAG_ROUGHNESS`, `wadDragMax`, `wadManningFile`, `wadMaxFroude`, `wadMaxSpeed`, `wadConserveVol`, `wadMonFreq`, `wadSmoothWidth`
  - L146 Diagnostics and monitor — `WADdryC`, `pickup_wad`, `pickup_wadface`
  - L155 Verification experiments — `wad_thacker_1d`, `wad_balzano`, `wad_flat_xz`, `wad_estuary_3d`, `wad_mudflat`, `wad_iceground`, `wad_iceground_sp`, `darwin_arm64_gfortran_conda`
  - L184 Cost and time step
  - L191 Limitations — `wadSmoothWidth`
  - L202 References

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

## `pkg/wad/verification/wad_iceground/README.md`
- L1 wad_iceground: Grounded sea ice that melts and floats
  - L9 Set-up
  - L18 Variants
  - L24 Inputs — `_part`, `bathy_deep`
  - L30 What to check

## `pkg/wad/verification/wad_iceground_sp/README.md`
- L1 wad_iceground_sp: Sea-ice growth and brine rejection on a drying flat — `wad_iceground`
  - L8 Set-up — `wad_iceground`, `data.salt_plume`, `salt_plume`, `dumpFreq`
  - L18 What to check — `WAD_DRY_FORCING`

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

## `regions/GoM/GoM_1km/GoM_1km_setup_changelog.md`
- L1 GoM_1km setup changelog — `theta_max`
  - L12 Surface forcing — ERA5 — `ERA5_u10m_filled_1992`, `useExfYearlyFields`, `theta_max`
  - L37 River forcing — GloFAS runoff and temperature — `exf_runoff_max`, `crop_glofas_for_exf.py`, `glofas_to_mitgcm.py`, `salt_min`, `make_glofas_river_temp_forcing.py`
    - L51 Why reducing runoff did not fix it — `empmr_min`, `exf_runoff_max`
  - L72 Bathymetry and vertical grid — `GoM_1km_bathymetry_8_maskshallow_lagoons_deepenNE.bin`, `diffKr`, `ggl90_calc_diff.F`, `empmr_min`
  - L104 Open boundary conditions
    - L109 The merge bug — `dv_diag_all`, `ETAN_south.bin`, `UVEL_south.bin`
    - L122 The cadence, established twice — `ETAN_south`, `_merge_dv_chunks`
    - L134 Changes — `ETAN_east.bin`, `read_dv_diags`, `UTILS_DIR`, `gen_obcs_fast.py`, `gen_obcs.py`
    - L144 Two notes — `obcs_cadence.py`, `check_obcs_records.py`, `gen_obcs_patch.py`
  - L162 Initial conditions and pickups — `hydrogThetaFile`, `hydrogSaltFile`, `uVelInitFile`, `vVelInitFile`, `pSurfInitFile`, `PTRACERS_initialFile`, `clip_ptracers.py`, `theta_min`, `_zinterp_fixed`
  - L183 Input file manifest — `check_run_inputs.py`, `GoM_1km_bathymetry_8_maskshallow_lagoons_deepenNE.bin`, `GoM_1km_it42_diffkr_r8.data`, `delYFile`, `pickup_THETA`, `pickup_SALT`, `pickup_U`, `pickup_V`, `pickup_ETAN`, `pickup_ggl90`, `_v10m_filled`, `_spfh2m_filled`, `_tmp2m_degC_filled`, `_dlw`, `_dsw`, `_rain`, `GloFas_GoM`, `GloFas_temp_GoM`, `iron_monthly_clim_Hamilton_kgFem2s_GoM_1km`, `_YYYY`, `apCO2`
    - L210 Where this run actually stops — `GloFas_GoM_1994`, `GloFas_temp_GoM_1994`
    - L223 What finishing the forcing costs
  - L247 data file settings — `bathyFile`, `diffKrFile`, `GoM_1km_bathymetry_8_maskshallow_lagoons_deepenNE.bin`, `GoM_1km_it42_diffkr_r8.data`, `diffKrT`, `diffKrS`, `startDate_1`, `startDate_2`, `nIter0`, `deltaT`, `endTime`, `nonlinFreeSurf`, `select_rStar`, `implicitFreeSurface`, `implicitDiffusion`, `implicitViscosity`, `eosType`, `rhoConstFresh`, `tempAdvScheme`, `saltAdvScheme`, `PTRACERS_advScheme`, `multiDimAdvection`, `staggerTimeStep`, `vectorInvariantMomentum`, `viscAr`, `viscC4Leith`, `viscC4Leithd`, `viscA4GridMax`, `useGGL90`, `useRealFreshWaterFlux` (+14)
    - L283 Known-bad, not yet fixed — `DIAGNOSTIC_SIZE.h`, `vel_3D_mon_mean`
  - L299 Watch list for the job now running — `sss_del2`, `uvel_max`, `theta_max`, `vvel_max`, `salt_mean`, `advcfl_W_hf`

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

## `regions/PointDume/README.md`
- L1 Point Dume downscaled regional set-up
  - L4 Preliminary information
  - L9 1. Get code
  - L21 2. Build executable
  - L33 3. Instructions for running simulation
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
