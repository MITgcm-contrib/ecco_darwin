# Southern California Bight 2-km nest + Malibu ~460-m kelp nest (Jan–Jun 2025)

This folder holds two nested ECCO-Darwin regional configurations for the January–June 2025 period:

1. `socal/`: a 1/48° (~2 km) Southern California Bight nest. Its initial conditions and hourly open boundaries
   come from the global ECCO-Darwin v05 LLC270 solution on ECCO V5r1. It adds barotropic tides at the open
   boundaries.
2. `malibu/`: a 0.005° (~460 m) Malibu / Point Dume nest with the G-MACMODS giant kelp model, pkg/radtrans
   spectral light and pkg/sediment. Its physical initial conditions, daily boundaries and tides come from
   `socal/`.

The SoCal nest is there because the Malibu nest at first took its boundaries directly from the 1/3° global
parent. That was a ~64× jump in grid spacing, with monthly-mean boundaries and no tides.

| Level | Grid | Levels | Time step | Ranks | darwin3 |
|---|---|---|---|---|---|
| ECCO-Darwin v05 LLC270 (V5r1) | 1/3° | 50 | | | |
| `socal/` | 168 x 101, lon -121.0 to -117.5, lat 32.5 to 34.6 | 60 (20 x 2 m to 40 m, then stretched to 3,961 m) | 120 s | 24 (24 x 1 tiles of 7 x 101) | `840fc4fe4` |
| `malibu/` | 128 x 72, lon -119.10 to -118.46, lat 33.80 to 34.10 | 50 (20 x 2 m to 40 m, then stretched) | 90 s | 32 (8 x 4 tiles of 16 x 18) | `24885b71` |

The parent global solution is described in Carroll et al. 2022 (GBC), https://doi.org/10.1029/2021GB007162.
The kelp model is G-MACMODS (Arzeno-Soltero et al. 2023, Commun. Earth Environ.,
https://doi.org/10.1038/s43247-023-00833-2), as set up for the CCS in `../v05_kelp`.

## Contents

| Path | What it is |
|---|---|
| `socal/code/` | Code changes for the SoCal build (`-mods`). `obcs_add_tides.F` ramps the boundary tides up over the first 2 days. |
| `socal/input/` | Namelists of the SoCal reference run (`run_v4_prod`). |
| `socal/scripts/` | Build, set-up and run scripts, and the Python scripts that made the SoCal grid, temperature correction and tides. |
| `malibu/code/` | Code changes for the Malibu build, used with `../v05_kelp/code` and `../v05/code`. Every change is tagged `C malibu_v01:`. |
| `malibu/pkg_sediment/` | The pkg/sediment version the Malibu reference run used (see note in section 4). |
| `malibu/input/` | Namelists and PBS script of the Malibu reference run (`run_2025h1_v4w_smg155`). |
| `malibu/scripts/` | The Python scripts that made the Malibu grid, initial and boundary conditions, rivers, waves, kelp initial condition and namelists. |

The scripts are the versions that made the reference inputs. They keep the paths of the development project,
so read each script header before you run it. The binary inputs of both runs are on Box and Pleiades (section 6), so
you do not need the scripts to repeat the reference runs.

## 1. Get the code

```
git clone https://github.com/MITgcm-contrib/ecco_darwin.git
git clone https://github.com/darwinproject/darwin3 darwin3_socal
git clone https://github.com/darwinproject/darwin3 darwin3_malibu
cd darwin3_socal  && git checkout 840fc4fe4 && cd ..
cd darwin3_malibu && git checkout 24885b71 && cd ..
cp -r ecco_darwin/regions/CCS/socal_malibu/malibu/pkg_sediment darwin3_malibu/pkg/sediment
```

The two nests use different darwin3 versions: `socal/` follows the LLC270 V5r1 parent (as `../../kerguelen/v05`),
and `malibu/` follows `../v05_kelp`.

## 2. Build `socal/` on Pleiades

```
module purge
module load comp-intel/2023.2.1 mpi-hpe/mpt.2.30 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12
cd darwin3_socal && mkdir build && cd build
../tools/genmake2 -rootdir .. -mpi -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas \
  -mo ../../ecco_darwin/regions/CCS/socal_malibu/socal/code
sed -i '/^LINK = /{/-mcmodel=medium/!s/$/ -mcmodel=medium -Wl,--no-relax/}' Makefile
sed -i '/^FFLAGS = /{/-mcmodel=medium/!s/$/ -mcmodel=medium/}' Makefile
sed -i '/^F90FLAGS = /{/-mcmodel=medium/!s/$/ -mcmodel=medium/}' Makefile
sed -i 's/ -axCORE-AVX2 -xSSE4\.2//' Makefile
make depend
make -j 16
```

`socal/scripts/job_build_tides.sh` is the PBS version of these steps. Run `make CLEAN` before you build again
after a change to a header file (`*.h`).

## 3. Run `socal/`

1. Get `socal_run_inputs_2025h1.tar` from Box (<https://nasa-ext.box.com/s/u7vmrhb5wggjuy4cer9icsjbaedytsrm>, 0.41 GB) and untar it. It gives `socal_run/` with the
   bathymetry, the initial conditions (`pickup_*.0000867888.data`), the boundary tides, the sea-filled ERA5
   atmosphere, the atmospheric CO2 and the iron dust files.
2. Link the hourly boundary files (on Pleiades, see section 6), copy the namelists and link the executable:

```
cd socal_run
ln -s /nobackup/dcarrol2/pub/regions/CCS/socal_malibu/socal_obcs_2025h1/*.bin .
cp ../ecco_darwin/regions/CCS/socal_malibu/socal/input/* .
ln -s ../darwin3_socal/build/mitgcmuv .
mkdir -p diags/daily diags/monthly diags/hourly
mpiexec -np 24 ./mitgcmuv
```

`socal/scripts/job_run.sh` runs it on one Skylake node. The run without tides took 28 h 19 min for the
181 days. Stop the run at or before `endtime = 15631200.` (4,342 h): the hourly boundary files end at
23:00 on 30 June 2025.

## 4. Build `malibu/` on Pleiades

```
module purge
module load comp-intel mpi-hpe/mpt hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12
cd darwin3_malibu && mkdir build && cd build
E=../../ecco_darwin/regions/CCS
../tools/genmake2 -mpi -of $E/v05_kelp/code/linux_amd64_ifort+mpi_ice_nas \
  -mo "$E/socal_malibu/malibu/code $E/v05_kelp/code $E/v05/code"
make depend
make -j 16
```

Notes:

1. darwin3 `24885b71` runs `tools/darwin/cog` with `python`. If `genmake2` stops at "pkg/darwin from
   templates", put a `python -> python3` link first on your `PATH`.
2. `malibu/pkg_sediment/` is the pkg/sediment of the reference run. It is older than `ecco_darwin/pkg/sediment`:
   three files are different (`SEDIMENT_OPTIONS.h` comment, and the pickup read and write routines).
3. The core hooks for pkg/sediment on this darwin3 version (`forward_step.F`, `packages_*.F`, `PARAMS.h`) are
   in `malibu/code/`.
4. `malibu/code/` has the seasonal kelp mortality (`mp_spp_death_mon`), the kelp drag on currents (off in the
   reference run), the kelp optics in pkg/radtrans, the ramped boundary tides and extra CDOM diagnostics.

## 5. Run `malibu/`

1. Get `malibu_run_inputs_2025h1.tar` from Box (<https://nasa-ext.box.com/s/4w7549u4767mmmmi0w4xzx113pzehq3h>, 2.2 GB) and untar it. It gives `malibu_run/` with the
   grid (`BATHY_128x72_malibu`, `delYFile_72_malibu`), `init/`, `obcs/`, `tides/`, `rivers/`, `waves/`,
   `sediment/`, `oasim/`, `era5_seafill/`, `v5r1/` (iron dust), `OPTICS_COEFF2/` and the `apCO2_*` files.
2. Copy the namelists and link the executable:

```
cd malibu_run
cp ../ecco_darwin/regions/CCS/socal_malibu/malibu/input/* .
ln -s ../darwin3_malibu/build/mitgcmuv .
mkdir -p diags/daily
qsub job_pace.pbs
```

Change `group_list` in `job_pace.pbs` to your own group. The 180-day reference run took 9 h 36 min on one
Skylake node (32 ranks).

## 6. Inputs on Box

| Where | What it holds |
|---|---|
| Box: <https://nasa-ext.box.com/s/4w7549u4767mmmmi0w4xzx113pzehq3h> (`malibu_run_inputs_2025h1.tar`, 2.2 GB) | All binary inputs of the Malibu reference run. |
| Box: <https://nasa-ext.box.com/s/u7vmrhb5wggjuy4cer9icsjbaedytsrm> (`socal_run_inputs_2025h1.tar`, 0.41 GB) | SoCal bathymetry, initial conditions, boundary tides, ERA5 atmosphere, atmospheric CO2, iron dust. |
| Pleiades `/nobackup/dcarrol2/pub/regions/CCS/socal_malibu/socal_obcs_2025h1/` | SoCal hourly boundary files: T (temperature-corrected), S, U, V and 31 Darwin tracers on 4 boundaries (140 files, 19.6 GB). |

`md5sums.txt` in each Box tar lists a checksum for every file. `sizes.txt` in the Pleiades folder lists the size of every file.

The Malibu biogeochemical initial conditions and boundaries use CalCOFI PRELIMINARY bottle data (cruises 2411SR,
2501RL, 2502RL). CalCOFI marks these data for non-publication use: oxygen, nitrate and chlorophyll can change
after the post-cruise calibrations.

How the inputs were made:

| Input | Script | Source |
|---|---|---|
| SoCal grid and bathymetry | `socal/scripts/make_socal_grid.py`, `make_socal_ncgrid_full.py` | NOAA Coastal Relief Model vol. 6 |
| SoCal initial conditions and hourly boundaries | `regions/downscaling` STEP1–STEP3 (`gen_obcs_fast.py`, `gen_pickups.py`) | LLC270 V5r1 rerun with pkg/diagnostics_vec, 1 Jan 2025 |
| SoCal temperature correction (initial conditions and boundaries, 0–75 m) | `socal/scripts/socal_ctd_bias.py`, `make_socal_ctdfix.py` | CalCOFI CTD casts, cruises 2501RL and 2504SH |
| SoCal boundary tides | `socal/scripts/socal_tides.py`, `calibrate_socal_tides.py` | FES2014a currents, calibrated against FES2022b sea level |
| Atmosphere (both nests) | `socal/scripts/era5_seafill_box.py` | ERA5 hourly, land cells filled from the nearest ocean cell |
| Malibu grid and bathymetry | `malibu/scripts/make_grid_bathy.py` | NOAA Coastal Relief Model vol. 6 |
| Malibu physics initial conditions and daily boundaries | `malibu/scripts/extract_socal_box.py`, `make_ic_obcs.py --mode socal` | `socal/` daily output |
| Malibu tides | `malibu/scripts/socal_tides_malibu.py` | `socal/` boundary tides |
| Malibu BGC initial conditions and boundaries | `malibu/scripts/calcofi_bgc.py`, `make_pace_obcs_phyto.py`, `make_obcs_community.py` | CalCOFI preliminary bottle data (cruises 2411SR, 2501RL, 2502RL) for upper-ocean NO3, PO4, SiO2 and plankton; PACE Chl scaling Feb–Apr; 70% diatom share Feb–Apr |
| Malibu TSS boundaries | `malibu/scripts/make_tss_obcs.py` | Zero (the parents carry no TSS) |
| Kelp initial condition (`init/Bmag`, `Qmag`, `Cmag`) | `malibu/scripts/make_kelp_ic.py` | KelpWatch Landsat canopy 2025 Q1 (EDI knb-lter-sbc.74), SBC LTER canopy fraction (knb-lter-sbc.112) |
| Rivers (water, sediment, nutrients) | `malibu/scripts/make_rivers.py` | LA County Public Works gauge F130-R (Malibu Creek) and USGS daily discharge; sediment rating curves |
| Waves | `malibu/scripts/make_waves.py --source waverys` | Copernicus WAVERYS (`cmems_mod_glo_wav_my_0.2deg_PT3H-i`) |
| Seabed sediment | `malibu/scripts/make_sed_bed.py` | usSEABED mud % |
| Light (OASIM) | `malibu/scripts/oasim_window.py` | OASIM C61 monthly fields below the surface, 2025 |
| Iron dust | `malibu/scripts/extract_v5r1_box.py` | LLC270 V5r1 (Mahowald et al. 2009) |
| Malibu namelists | `malibu/scripts/make_namelists.py`, `make_inputs_2025h1.py` | |

## 7. Kelp settings of the reference run

| | Value | Source |
|---|---|---|
| Maximum growth (`mp_spp_Gmax_cap`) | 0.155 d⁻¹ | Lowered from 0.173 d⁻¹ so that the Q1 → Q2 2025 growth is near KelpWatch |
| Background mortality (`mp_spp_death_mon`) | 0.017–0.042 d⁻¹ by month, linear between mid-month values | SBC LTER monthly median loss (knb-lter-sbc.112.15) minus the model wave loss of 0.0035 d⁻¹ |
| Kelp drag on currents | off | `data.kelpdrag` |
| Initial kelp | 106 t dry weight on 1 January 2025 | KelpWatch 2025 Q1 |

## 8. Evaluation (1 January to 30 June 2025)

| Observation | Model |
|---|---|
| NOAA tide gauges Santa Monica and Los Angeles (SoCal, 15-day calibrated tide test) | hourly sea level r = 0.996, M2 amplitude x1.03–1.04 of observed, phase +1° |
| KelpWatch Malibu canopy: 106 t (Q1), 299 t (Q2), growth x2.82 | 107 t (15 Feb), 322 t (16 May), growth x3.00 |
| SBC LTER kelp loss Apr–Jun: 0.032, 0.023, 0.024 d⁻¹ | 0.031, 0.024, 0.023 d⁻¹ |
| Malibu CTD casts (n = 8), median model − CTD | 2 m −0.06 °C, 10 m +0.10 °C, 20 m +0.45 °C, 30 m +0.99 °C (measured on run v4, before the BGC-boundary, wave and kelp changes; not checked again) |
| PACE surface Chl, shelf, 8–23 March | mean log10(model/PACE) = −0.34 |

## 9. Known limitations

1. Spring shelf warm bias at 20–30 m, about +1 °C at the late-March CalCOFI casts. The SoCal nest has it too.
   Lower GGL90 mixing and tides off did not remove it.
2. The March inner-shelf bloom is too small (shelf Chl about half of PACE).
3. Model TSS is mineral sediment only and is far below OLCI TSM near the coast.
4. Kelp loss in January–March is lower than SBC LTER (0.027–0.034 vs 0.033–0.046 d⁻¹).
5. The kelp initial condition depends on the canopy fraction of the SBC LTER sites, which are shallower
   (4.5–10.7 m) than much of the Malibu kelp. This fraction is uncertain at Malibu.
6. SoCal sea level has no remote signal: there is no sea-level boundary condition, and the volume balance holds
   the domain mean fixed.
