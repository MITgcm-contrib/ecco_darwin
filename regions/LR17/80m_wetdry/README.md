# LR17 80 m with wetting and drying (baseline v4.1)

LR17 is an ~80 m MITgcm configuration of the Pertuis Charentais (Bay of Biscay, France). This folder holds
the v4.1 baseline: physics with wetting and drying (pkg/wetdry) on Adeola's LR17 80 m grid
(`../80m`), with real March 2022 forcing. pkg/sediment is compiled in and is off by default.

The reference run is 43 days from 23 February 2022 00 UTC, which covers the observation window
of the 2022 calibration package (SOMLIT Antioche, CocoriCO2, AWAC, SHOM tide gauges, IFREMER SST).

## Contents

| Path | What it is |
|---|---|
| `code/` | Code changes for the build (`-mods`). `SIZE.h` is for 84 MPI ranks (12 x 7 tiles); `SIZE.h.10ranks` is for 10 ranks (2 x 5 tiles). |
| `code/lr17_flather.F` | Flather radiation correction of the prescribed normal velocity on the west and south boundaries. |
| `code/lr17_rivers.F`, `LR17_RIVERS.h` | Daily river discharge (`data.rivers`) as volume sources (`addMass`) at the tip of each river channel. |
| `code/obcs_calc.F`, `code/load_fields_driver.F` | Stock MITgcm `d861cd501` files with one added call each (to `LR17_FLATHER` and `LR17_RIVERS`). |
| `input/` | Namelists of the v4.1 reference run. |
| `scripts/` | Python scripts that make every binary input from public data (see "Make the inputs"). |

## 1. Get the code

The build uses stock MITgcm at commit `d861cd501` with pkg/wetdry and pkg/sediment from ecco_darwin
(tested with ecco_darwin `09d448a`; pkg/wetdry as of `ef048f0`). LR17 has no sea ice and does not use the
landfast-ice options of pkg/wetdry.

```
git clone https://github.com/MITgcm/MITgcm.git
git clone https://github.com/MITgcm-contrib/ecco_darwin.git
cd MITgcm
git checkout d861cd501
cp -r ../ecco_darwin/pkg/wetdry/src pkg/wetdry
git apply ../ecco_darwin/pkg/wetdry/mitgcm_changes/wetdry_hooks.patch
cp -r ../ecco_darwin/pkg/sediment/src pkg/sediment
git apply ../ecco_darwin/pkg/sediment/mitgcm_changes/sediment_hooks_with_wetdry.patch
```

If a patch does not apply, see `ecco_darwin/pkg/wetdry/README.md` and `ecco_darwin/pkg/sediment/README.md`.

## 2. Build on Pleiades (84 ranks)

```
module purge
module load comp-intel mpi-hpe hdf4 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
mkdir build && cd build
../tools/genmake2 -mods=../../ecco_darwin/regions/LR17/80m_wetdry/code \
  -of=../tools/build_options/linux_amd64_ifort+mpi_ice_nas -mpi
make depend
make -j 40
```

Run `make CLEAN` before you build again after a change to a header file (`*.h`).

## 3. Run

1. Get the run template from Box: <https://nasa-ext.box.com/s/211e45aj2ve5j7raqbvf8yty8vkex602>. It holds `80m_wetdry_run_template.tar` (binary inputs, 0.53 GB;
   `tar -xf` gives `run_template/`).
2. Copy the namelists and link the executable:

```
cd run_template
cp ../ecco_darwin/regions/LR17/80m_wetdry/input/* .
ln -s ../MITgcm/build/mitgcmuv .
mpiexec -np 84 ./mitgcmuv
```

On Pleiades the 43-day run takes about 33 hours on 3 Skylake nodes (84 ranks). Use the `long` queue.
LR17 needs at least 10 ranks: a 6-rank build stops at start-up on macOS.

## 4. Make the inputs

All inputs come from public data. The scripts expect the data in `data/LR17/` next to `scripts/`
(or set `LR17_ROOT` to the folder that holds `data/LR17/`). The files that no script downloads are in the
Box data tarball `80m_wetdry_data.tar` (<https://nasa-ext.box.com/s/zd2oinz1k106ixuwq9a34zr3lfrvqqin>; `tar -xf` gives `data/LR17/`): Adeola's `bathy_LR17_936x875`, the
OpenStreetMap coast, dike and river files, the Hub'Eau discharge used for the reference run, the v4 bed
`LR17_zb_v4.bin` and the EMODnet bed maps for sediment runs.

| Step | Script | Data | Account |
|---|---|---|---|
| IBI ocean reanalysis (sea level, currents, T, S, waves) | `LR17_ibi.py START END` | CMEMS `IBI_MULTIYEAR_PHY_005_002`, `IBI_MULTIYEAR_WAV_005_006` | Copernicus Marine |
| Atmosphere | `LR17_era5.py START END` | ERA5 hourly single levels | Copernicus CDS (`~/.cdsapirc`) |
| Rivers | `LR17_rivers.py START END` | Hub'Eau daily discharge | none |
| FES tide atlases | `fes_crop.py FILE ...` | FES2022b elevations, FES2014 currents | AVISO |
| Bed (v4 river channels) | `LR17_bathy.py --zmouth=-4 --ztip=-3 --out=LR17_zb_v4.bin` | Adeola's `bathy_LR17_936x875`, OpenStreetMap coast, dikes and rivers | none |
| Run folder | `LR17_gendata_somlit.py` (calls `gendata.py`) | all of the above | |

The command for the v4.1 reference run:

```
LR17_S_OFFSET_IC=1.3 python3 scripts/LR17_gendata_somlit.py 2022-02-23 43 run_v41 2022-02-20 2022-04-10 \
  --fes-boundary=eta --band-download=2022-01-25_2022-05-05 --zb=data/LR17/work/LR17_zb_v4.bin
cp input/data* run_v41/
```

`--band-download` needs IBI sea level and currents for 25 January to 5 May 2022, so that the band-pass
filter does not leak into the run window. Then copy the namelists from `input/`: they add the hourly
level-1 U, V output (`uv1`) that `gendata.py` does not write. We checked that these scripts make 30 of
the 31 input files of the reference run byte for byte; the other one is `data.diagnostics`.

## Sediment (optional)

pkg/sediment is compiled in. To make a sediment run, add `--sed` to the `LR17_gendata_somlit.py` command:
it writes `data.sediment` and `data.ptracers`, the initial bed from the EMODnet seabed substrate map
(`LR17_bed_mud.bin`, `LR17_bed_sand.bin`), the mud carried by each river and the IBI waves, and it turns on
rain and evaporation. A first 2-day test ran on an earlier physics version; sediment is not yet
calibrated.

## 5. What is in v4.1

| | Adeola's `../80m` | v4.1 (this folder) |
|---|---|---|
| Period | 48 h from 16 Jan 1993 | 43 days from 23 Feb 2022 |
| Vertical grid | 9 levels, 60 m | 14 levels (8 x 2.5, 2 x 5, 3 x 10, 8.5 m), rest level 5.5 m above MSL |
| Free surface | surface level, no wetting and drying | r*, pkg/wetdry (5 cm film, 10 cm critical depth) |
| Bathymetry | `bathy_LR17_936x875` | the same, plus intertidal flats and marsh up to +5 m MSL and river channels from -4 m MSL (mouth) to -3 m MSL (tip) |
| Initial T, S | IBI | IBI daily mean of the start day, S + 1.3 (IBI is about 1.75 too fresh at SOMLIT Antioche) |
| Boundary sea level | IBI | IBI, with its semidiurnal and quarter-diurnal tide replaced by FES2022b |
| Boundary currents | IBI | IBI depth-mean current and daily shear, with a Flather radiation condition |
| Atmosphere | ERA-Interim, 6-hourly | ERA5 hourly |
| Rivers | none | Charente (with the Boutonne), Sevre Niortaise, Lay, Seudre (Hub'Eau, daily) |
| Salt vertical advection | | explicit, flux-limited third order (scheme 33): no negative salinity at the river sources |

The full description, with the tests that led to each choice, is in
`docs/LR17/LR17_boundary_initial_conditions.md` of the development project.

## 6. Evaluation of the reference run (23 February to 7 April 2022)

| Observation | v4.1 |
|---|---|
| SHOM tide gauges La Pallice / Ile d'Aix | sea-level RMSE 0.14 / 0.15 m for a ~6 m tidal range (IBI: 0.21 / 0.22 m) |
| Tidal constituents (41 days) | M2 0.07 m low (~4%), phase within 1.5 degrees; M4 within 0.04 m and 2.5 degrees; MS4 about half of observed |
| AWAC 2 current speed (hourly, level 1) | bias +0.03 m/s, rms 0.16 m/s, r = 0.84 (observed mean 0.35 m/s) |
| CocoriCO2 MA21 / ROCHL (daily) | temperature bias -0.06 / +0.10 C; salinity bias -0.89 / -1.07 |
| IFREMER L4 SST | bias +0.15 C (domain mean) |
| SOMLIT Antioche casts (4) | mean salinity error 0.44 (IBI: 1.31); temperature within -0.33 to +0.08 C |
| Negative salinity | no surface cell below S = 0 in the 6-hourly output (v4 with implicit vertical advection: 178 cells) |

## 7. Known limitations

1. M2 at the SHOM gauges is about 4% low; most of this comes from the IBI and FES forcing.
2. Spring low waters are about 0.1 m too high, and MS4 has about half of the observed amplitude.
   Lower and zoned bottom friction, higher flats and a thinner wet/dry critical depth did not change this.
3. Salinity is 0.4 to 1.1 too fresh at the observation sites, mostly from IBI.
4. The river sources are single cells at the tip of each channel.
