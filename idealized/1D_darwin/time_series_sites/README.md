# 1-D ECCO-Darwin v05 columns at ocean time-series sites

One-dimensional ECCO-Darwin v05 water columns, extracted from the ECCO-v5r1 (llc270) solution, at five ocean
time-series sites. The columns run 1992–2025 in about an hour on a single laptop core, which makes them a fast,
data-constrained testbed for the Darwin ecosystem: they come with an observation set for each site and a
Green's-function workflow for optimizing ecosystem parameters against those observations.

| Site | Station | Model point |
|---|---|---|
| `HOT` | Station ALOHA, Hawaii Ocean Time-series | 22.75°N, 158.00°W |
| `BATS` | Bermuda Atlantic Time-series Study | 31.67°N, 64.17°W |
| `HydroS` | Hydrostation "S", Bermuda | 32.17°N, 64.50°W |
| `PAPA` | Ocean Weather Station Papa | 50.10°N, 144.90°W |
| `PAP` | Porcupine Abyssal Plain Sustained Observatory | 49.00°N, 16.50°W |

Each column sits at the nearest wet llc270 point to the station.

**Parent solution:** ECCO-Darwin v05 llc270 v5r1 (`v05/llc270/readme_v5r1.txt`), documented in
Carroll et al. (2022), *GBC*, 36, e2021GB007162, <https://doi.org/10.1029/2021GB007162>.
**Ecosystem:** five phytoplankton and two zooplankton types, as in Carroll et al. (2020), *JAMES*, 12,
e2019MS001888, <https://doi.org/10.1029/2019MS001888>.

## Set-up

- 1×1 spherical-polar cell at the llc270 point, 50 levels (v5r1 `delR`), f-plane, 1200-s time step.
- darwin3 `darwin_ckpt68g`, with the same Darwin code, namelists and traits as v05/llc270 (31 tracers).
- GGL90 mixing (v5r1 `data.ggl90`) plus the v5r1 `xx_diffkr` background diffusivity at the site.
- v5r1 iter70 6-hourly atmospheric state (air temperature, humidity, precipitation, winds, downward short- and
  longwave) sampled at the site; bulk formulae and zenith-angle shortwave as in v5r1; `pkg/seaice` thermodynamics.
- Darwin forcing: Mahowald (2009) soluble dust Fe and NOAA MBL atmospheric pCO2 sampled at the site.
- Initial T, S and all 31 ptracers from v5r1 pickup `0000000002` (1992-01-01).
- **T/S relaxation** (`pkg/rbcs`, 30-day timescale, full column) toward v5r1 monthly-mean profiles at the site, as a
  stand-in for 3-D advection.
- **Deep BGC supply:** NO3, PO4, FeT, SiO2, DIC, ALK and O2 are relaxed (1-year timescale) toward v5r1 monthly
  profiles below the euphotic zone only (mask 0 above 150 m, linear ramp, 1 below 250 m). The upper 150 m and all
  plankton, DOM and POM tracers evolve freely.
- **rbcs timing:** with `rbcsForcingCycle=0`, record *k* is centered at
  `rbcsForcingOffset + (k-0.5)*rbcsForcingPeriod` (`eesupp/src/get_periodic_interval.F`), so
  `rbcsForcingOffset = -2629746` (one period) centers record 2 on mid-January 1992.
- **Freshwater:** virtual salt flux and linear free surface (`useRealFreshWaterFlux=.FALSE.`), so local E–P does
  not drain or fill the column over 34 years. `PTRACERS_EvPrRn` is left unset, so tracers see no local E–P. (With
  `EvPrRn=0`, surface DIC and ALK drifted by roughly +90 mmol m⁻³ at HOT and −95 mmol m⁻³ at Papa in 16 years.)
- 1992-01-01 to 2025-12-31 (`endTime` as v5r1), with daily-mean NetCDF diagnostics (`pkg/mnc`).

## Directory contents

| Path | Contents |
|---|---|
| `code/` | MITgcm/darwin3 overrides: v05/llc270 `code_darwin` plus the column-relevant v05/llc270 `code_v5r1` files, a 1×1×50 `SIZE.h` and `packages.conf` (adds `rbcs`, `mnc`) |
| `code_rbcs_mask21/` | Optional `RBCS_SIZE.h` (`maskLEN = 21`) for per-tracer relaxation masks (`relaxMaskFile(2+iTr)`). Add it to `code/` and `make CLEAN`. Bit-identical to the standard build when all masks are equal (10-day test). |
| `input/` | Namelist templates (`data` has `@SITE@`, `@F0@`, `@XG@`, `@YG@` placeholders) |
| `scripts/` | Set-up, Green's-function optimization and experiment scripts (sections 3–6) |
| `pfe/` | PBS scripts used on NASA Pleiades (edit the `/nobackup` paths) |
| `obs/` | Readers that build the tidy per-site observation files (section 7) |
| `results/` | Optimized parameter sets and the forward-run misfit table (section 6) |

### Site input files

The site input files are **not in the repo**. They are on NASA Pleiades, readable by all users:

```
/nobackup/dcarrol2/pub/1-D/sites/<SITE>/     # HOT, BATS, HydroS, PAPA, PAP; ~12 MB each
/nobackup/dcarrol2/pub/1-D/README.txt        # what each file is
```

Each site folder holds the EXF forcing, apCO2, dust Fe, background diffusivity, bathymetry, initial T/S/ptracers,
rbcs relaxation profiles and v5r1 validation series. Copy them with, e.g.,

```
scp -r pfe:/nobackup/dcarrol2/pub/1-D/sites .
```

or regenerate them from v5r1 (section 4).

## 1. Get code

```
git clone https://github.com/MITgcm-contrib/ecco_darwin.git
git clone https://github.com/darwinproject/darwin3
cd darwin3
git checkout darwin_ckpt68g
mkdir build_1D
cd build_1D
```

## 2. Build executable

The build is serial; a 1992–2025 column takes about an hour on one laptop core.

```
MOD="../../ecco_darwin/idealized/1D_darwin/time_series_sites"
../tools/genmake2 -rootdir=.. -mods=$MOD/code -of=<your optfile>
make depend
make -j 4
```

> **macOS:** darwin3's `cog` step needs a `python` older than 3.12 (it imports `imp`) on `PATH`, e.g. a small
> wrapper that runs `/usr/bin/python3`. With a conda gfortran, the environment's `bin/` must also be on `PATH` so
> the assembler finds `clang`.

## 3. Make run directories and run

```
cd ..
SITES=/nobackup/dcarrol2/pub/1-D/sites       # on Pleiades; elsewhere, your copy of it
python3 $MOD/scripts/make_runs.py $SITES $MOD/input runs_1D build_1D/mitgcmuv
cd runs_1D/HOT
./mitgcmuv > output.txt
```

Repeat for `BATS`, `HydroS`, `PAPA` and `PAP`. Daily-mean diagnostics are written to `mnc_out/`.
`scripts/run_queue.sh` runs several columns at once (`MAXJOBS`, default 4).

## 4. (Re)extract site inputs from v5r1 (NASA Pleiades only)

Edit `SITES` in `scripts/extract_v5r1_sites.py` to add a site, then:

```
qsub scripts/job_extract
python3 scripts/make_relax.py sites/*        # rbcs T/S + BGC relaxation profiles
```

## 5. Green's-function optimization

There are 20 controls (`scripts/gf_controls.py`):

- 19 ecosystem parameters, each a multiplicative perturbation of `data.traits` or `data.darwin`: maximum growth,
  nutrient half-saturation, maximum Chl:C, phytoplankton and zooplankton mortality, grazing, PIC:POC, POM sinking,
  POM and DOM remineralization, PIC dissolution, Fe scavenging, dust-Fe solubility, light attenuation, quantum
  yield and picoplankton palatability;
- `KZ_BG`, the background diffusivity at 75–250 m (edits `diffkr_1x1x50`), a stand-in for the eddy and
  internal-wave nutrient supply that a 1-D column cannot resolve.

**a) Perturbation runs:** one control run plus 19 perturbations per site, i.e. 100 runs (see `pfe/job_gf2`), plus
5 `KZ_BG` runs on any machine, differenced against the matching baseline run.

```
python3 scripts/gf_make_runs.py runs/ctrl runs/gf
python3 scripts/gf_kz_runs.py
```

**b) Bin the observations** (site × variable × depth layer × month) and compute model equivalents:

```
python3 scripts/gf_bins.py export <obs_root> gf_results/bins
python3 scripts/gf_bins.py equiv gf_results/bins runs/ctrl runs/gf gf_results/cache    # pfe/job_gf_equiv
python3 scripts/gf_bins.py cache gf_results/bins gf_results/cache
python3 scripts/gf_kz_cache.py gf_results/cache gf_results/cache_kz gf_results/bins
```

**c) Solve:** global, hold-one-site-out, single site, or a comma-separated list of sites (regional):

```
python3 scripts/gf_solve.py runs/ctrl none <obs_root> out/global --cache gf_results/cache_kz --kz --min_layer_bins 3
python3 scripts/gf_solve.py ... --only HOT,BATS,HydroS      # gyre set
python3 scripts/gf_solve.py ... --only PAPA,PAP             # subpolar set
```

The cost is the mean normalized misfit (d − m)²/r for each site–variable pair, with every pair weighted
equally. The error variance r is the larger of the bin standard deviation and 10% of the mean, squared, per site,
variable and depth layer. Pairs with fewer than 20 bins and layers with fewer than 3 bins are dropped, and POC flux
uses 100–500 m sediment traps only.

**d) Confirm** with forward runs using the optimized set, scored on the same cost:

```
python3 scripts/gf_make_opt.py out/global/gf_parameters.csv runs/ctrl runs/opt global HOT BATS HydroS PAPA PAP --fix KPOM
python3 scripts/score_exp.py out/global gf_results/bins gf_results/cache_kz runs/opt/global_HOT:HOT ...
python3 scripts/gf_importance.py <results_dir> gf_results/cache      # which controls matter at each site
```

## 6. Results

Forward 1992–2025 runs, scored on the cost of 2026-10-08. Mean normalized misfit, averaged over the five sites:

| Parameter set | Misfit | Notes |
|---|---|---|
| Baseline (control parameters) | 2.93 | |
| Global set, 19 controls | 2.42 | diatom max growth ×1.17 |
| Global set + background mixing, 20 controls | 2.23 | mixing ×2.5 |
| Regional sets (gyre: HOT, BATS, HydroS; subpolar: PAPA, PAP) | **2.05** | |
| Site-only sets | 2.07 | 1.98 with the PAP set at HydroS |

- **Gyre vs subpolar:** the gyre fit wants phytoplankton mortality ×0.59 and mixing ×2.4; the subpolar fit wants
  mortality ×1.33, mixing ×1.5 and diatom maximum growth ×1.23.
- **Files:** `results/params_*.csv` (column `factor` multiplies the control value) and
  `results/forward_run_misfit_by_variable.csv` (misfit for every run and variable).
- **Known structural biases:** surface primary production (~5× low) and export (~2.5× low) in the subtropical
  gyres, and surface SiO2 at OWS Papa, where the 1-D column lacks Ekman upwelling of the nutricline.
- **Follow-up experiments:** `scripts/make_exp.py` (mixing, picophytoplankton growth, Fe half-saturation, an
  upwelling proxy) and `scripts/make_papa_si.py` (SiO2-only relaxation and diatom Si:C; needs `code_rbcs_mask21`).

## 7. Observations

`obs/scripts/standardize_obs.py` builds `obs/<SITE>/<SITE>_obs.csv` (time, depth, variable, value, units, source,
QC flag) from HOT-DOGS and BCO-DMO, BATS/BIOS, Line P and the PMEL Papa mooring, OceanSITES PAP, BODC PAP cruises,
GLODAPv2.2023, SOCATv2026, BGC-Argo, OC-CCI v6.0 and GEOTRACES IDP2025. Sources, licenses and download steps are
in `obs/README_obs.md`; the raw data are not in the repo.

## Contact

Dustin Carroll (dustin.carroll@sjsu.edu)
