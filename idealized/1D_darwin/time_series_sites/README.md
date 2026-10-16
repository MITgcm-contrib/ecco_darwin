# 1-D ECCO-Darwin v05 columns at ocean time-series sites

One-dimensional ECCO-Darwin v05 water columns, extracted from the ECCO-v5r1 (llc270) solution, at 11 ocean
time-series sites, from the subtropical gyres to the Southern Ocean and the Nordic Seas. The columns run 1992–2025 in
about an hour on a single laptop core, which makes them a fast, data-constrained testbed for the Darwin ecosystem:
each site comes with an observation set and a Green's-function workflow for optimizing ecosystem parameters against
those observations. We provide both the v05 baseline and our current optimized set-ups, which lower the mean
model-data misfit at nine sites from 2.94 to 1.44 (section 6).

| Site | Station | Station position | Region (optimized set) |
|---|---|---|---|
| `HOT` | Station ALOHA, Hawaii Ocean Time-series | 22.75°N, 158.00°W | gyre |
| `BATS` | Bermuda Atlantic Time-series Study | 31.67°N, 64.17°W | gyre |
| `HydroS` | Hydrostation "S", Bermuda | 32.17°N, 64.50°W | gyre |
| `PAPA` | Ocean Weather Station Papa | 50.10°N, 144.90°W | high latitude |
| `PAP` | Porcupine Abyssal Plain Sustained Observatory | 49.00°N, 16.50°W | high latitude |
| `SOTS` | Southern Ocean Time Series | 46.94°S, 142.28°E | high latitude |
| `K2` | JAMSTEC K2, western subarctic Pacific | 47.00°N, 160.00°E | high latitude |
| `MUNIDA` | Munida offshore subantarctic station, New Zealand | 45.83°S, 171.53°E | high latitude |
| `OWSM` | Ocean Weather Station M, Norwegian Sea | 66.00°N, 2.00°E | high latitude |
| `CARIACO` | Cariaco Basin, Venezuela | 10.50°N, 64.67°W | none (v05 + closure) |
| `DYFAMED` | DYFAMED, Ligurian Sea | 43.42°N, 7.87°E | none (v05 + closure) |

Each column sits at the nearest wet llc270 point to the station (7–20 km away).

**Parent solution:** ECCO-Darwin v05 llc270 v5r1 (`v05/llc270/readme_v5r1.txt`), documented in
Carroll et al. (2022), *GBC*, 36, e2021GB007162, <https://doi.org/10.1029/2021GB007162>.
**Ecosystem:** five phytoplankton and two zooplankton types, as in Carroll et al. (2020), *JAMES*, 12,
e2019MS001888, <https://doi.org/10.1029/2019MS001888>.

## Quick start: run the optimized set-ups

1. Get the site inputs (not in the repo; see "Site input files" below) and build `mitgcmuv` (sections 1–2).
2. Build all 11 run directories with the current best parameters in one command:

   ```
   MOD=ecco_darwin/idealized/1D_darwin/time_series_sites
   python3 $MOD/scripts/make_optimized_runs.py <sites_dir> $MOD/input runs_opt build_1D/mitgcmuv
   ```

   Add site names at the end (e.g. `HOT K2`) to build only those. Needs `numpy` and `pandas`.
3. Run each column:

   ```
   cd runs_opt/HOT
   ./mitgcmuv > output.txt          # ~1 hr; daily-mean diagnostics in mnc_out/
   ```

   or several at once: `MAXJOBS=4 $MOD/scripts/run_queue.sh runs_opt HOT BATS HydroS ...`.

Each run directory holds an `opt_manifest.txt` that lists every namelist change, and `experiment.txt`. The final
`data.traits` and `data.darwin` for each site are also in `results/optimized/namelists/<SITE>/`, so you can see
or copy the parameters without building anything. For the v05 baseline instead, use `scripts/make_runs.py` (section 3).

## Set-up

- 1×1 spherical-polar cell at the llc270 point, 50 levels (v5r1 `delR`), f-plane, 1200-s time step.
- darwin3 `darwin_ckpt68g`, with the same Darwin code, namelists and traits as v05/llc270 (31 tracers).
- GGL90 mixing (v5r1 `data.ggl90`) plus the v5r1 `xx_diffkr` background diffusivity at the site.
- v5r1 iter70 6-hourly atmospheric state (air temperature, humidity, precipitation, winds, downward short- and
  longwave) sampled at the site; bulk formulae and zenith-angle shortwave as in v5r1; `pkg/seaice` thermodynamics.
- Darwin forcing: Mahowald (2009) soluble dust Fe and NOAA MBL atmospheric pCO2 sampled at the site.
- Initial T, S and all 31 ptracers from v5r1 pickup `0000000002` (1992-01-01). The v5r1 initial conditions and
  relaxation profiles contain small negative tracer values (71–268 per site); `make_runs.py` clips them to zero,
  because a plankton type that starts negative never recovers and stays extinct for all 34 years.
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

### Grazer + nitrogen-fixer closure (optimized set-ups only)

With v05 parameters, the gyre zooplankton sit below their break-even prey concentration (~1.2 mmol C m⁻³) and starve,
in the 1-D columns and in 3-D v05 alike. The optimized set-ups therefore add a closure (`scripts/make_closure_runs.py`):

1. **Quadratic zooplankton mortality:** linear mortality `MORT(6:7) = 0` and `MORT2(6:7) = MORT/0.1 mmol C m⁻³`, which
   keeps both grazers alive at all 11 sites;
2. **Nitrogen fixer:** plankton type 2 becomes a non-quota diazotroph (`DIAZO = 1`, `USENH4/NO2/NO3 = 0`; integer
   traits) with a slow maximum growth rate of 0.41 d⁻¹ (`PCMAX(2) = 4.70e-6 s⁻¹`, as in v06); its initial
   carbon is raised to at least 1e-3 mmol C m⁻³ (with matching Chl), since it starts at zero at most sites.

The nitrogen fixer persists only at the three subtropical gyre sites (9–12 mmol C m⁻² over 0–150 m in 2021–2025),
the expected biogeography. With darwin3's non-quota fixer, nitrogen fixation does not respond to ambient DIN; the
faster v05 type 2 (2.4 d⁻¹) drives surface NO3 to 5–50 mmol m⁻³ at HOT, which is why the slow rate matters.

## Directory contents

| Path | Contents |
|---|---|
| `code/` | MITgcm/darwin3 overrides: v05/llc270 `code_darwin` plus the column-relevant v05/llc270 `code_v5r1` files, a 1×1×50 `SIZE.h` and `packages.conf` (adds `rbcs`, `mnc`) |
| `code_rbcs_mask21/` | Optional `RBCS_SIZE.h` (`maskLEN = 21`) for per-tracer relaxation masks (`relaxMaskFile(2+iTr)`). Add it to `code/` and `make CLEAN`. Bit-identical to the standard build when all masks are equal (10-day test). |
| `input/` | Namelist templates (`data` has `@SITE@`, `@F0@`, `@XG@`, `@YG@` placeholders) |
| `scripts/` | Set-up, optimized set-ups, Green's-function optimization, scoring and experiment scripts (sections 3–6) |
| `pfe/` | PBS scripts used on NASA Pleiades (edit the `/nobackup` paths) |
| `obs/` | Readers that build the tidy per-site observation files (section 7) |
| `results/` | Optimized parameter sets, final namelists and misfit tables (section 6) |

### Site input files

The site input files are **not in the repo**. They are on NASA Pleiades, readable by all users:

```
/nobackup/dcarrol2/pub/1-D/sites/<SITE>/     # all 11 sites above; ~12 MB each
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

## 3. Make v05 baseline run directories and run

```
cd ..
SITES=/nobackup/dcarrol2/pub/1-D/sites       # on Pleiades; elsewhere, your copy of it
python3 $MOD/scripts/make_runs.py $SITES $MOD/input runs_1D build_1D/mitgcmuv
cd runs_1D/HOT
./mitgcmuv > output.txt
```

Repeat for the other sites. Daily-mean diagnostics are written to `mnc_out/`.
`scripts/run_queue.sh` runs several columns at once (`MAXJOBS`, default 4). A finished run prints `NORMAL END`
(`STOP NORMAL END` with some compilers).

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

With the closure in the control runs, set `GF_CONTROLS=closure`: zooplankton mortality is then perturbed through
`MORT2` (the linear `MORT` is zero), and the shared growth control of types 2–3 splits into separate
*Synechococcus* (`PCMAX_Syn`) and nitrogen-fixer (`PCMAX_diazo`) controls. `pfe/job_gf_cl2` + `pfe/setup_cl2.csh`
are the PBS scripts for that ensemble (231 runs at 11 sites).

**a) Perturbation runs:** one control run plus one run per control at each site (see `pfe/job_gf2`,
`pfe/job_gf_cl2`), plus the `KZ_BG` runs on any machine, differenced against the matching control run.

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
python3 scripts/gf_solve.py ... --only HOT,BATS,HydroS                       # gyre set
python3 scripts/gf_solve.py ... --only PAPA,PAP,SOTS,K2,MUNIDA,OWSM          # high-latitude set
```

The cost is the mean normalized misfit (d − m)²/r for each site–variable pair, with every pair weighted
equally. The error variance r is the larger of the bin standard deviation and 10% of the mean, squared, per site,
variable and depth layer. Pairs with fewer than 20 bins and layers with fewer than 3 bins are dropped, and POC flux
uses 100–500 m sediment traps only. Observations deeper than the column's sea floor have no model equivalent and
are skipped (since 2026-10-15; before, they were compared with the zero-filled dry levels, which inflated the misfit
at Hydrostation S by ~0.9 and at Munida by ~0.4 but changed the solved parameters by less than 0.001).

**d) Confirm** with forward runs using the optimized set, then score any set of runs on the same cost:

```
python3 scripts/gf_make_opt.py out/global/gf_parameters.csv runs/ctrl runs/opt global HOT BATS HydroS PAPA PAP --fix KPOM
python3 scripts/score_sites.py scores.csv runs/opt/global_HOT:HOT runs/opt/global_BATS:BATS ...
python3 scripts/gf_importance.py <results_dir> gf_results/cache      # which controls matter at each site
```

## 6. Results

Forward 1992–2025 runs, scored on the cost of 2026-10-15. Mean normalized misfit over the nine sites with a regional
set (`results/misfit_by_site_20261015.csv`; per variable in `results/misfit_by_site_variable_20261015.csv`):

| Parameter set | Misfit (9 sites) | Notes |
|---|---|---|
| v05 baseline | 2.94 | clipped initial conditions |
| Regional sets (gyre; high latitude) | 1.69 | step 1: one Green's-function fit per region |
| Regional sets + closure | 1.54 | grazers survive at all 11 sites |
| **Optimized** (regional sets + closure + second iteration) | **1.44** | step 2: a small Green's-function step around the line above |

The optimized set-ups improve every site relative to the v05 baseline, e.g. HOT 4.08 to 1.90, K2 4.45 to 1.49
and PAP 3.18 to 1.01. CARIACO (2.95) and DYFAMED (11.3) get the closure only; their misfit is structural (anoxic
basin, Mediterranean physics) and parameters do not fix it.

- **What the optimized sets change** (step 1 × step 2, relative to v05): the gyre set uses ~3.4× background mixing
  at 75–250 m, 0.63× phytoplankton mortality, 1.9× faster POM remineralization, 0.71× POM sinking and 1.2× diatom
  growth; the high-latitude set uses ~2.3× mixing, 2.8× faster POM remineralization, 1.6× phytoplankton mortality
  and 0.57× POM sinking. Every factor is in `results/optimized/{gyre,highlat}_step{1,2}.csv` (column `factor`
  multiplies the control value).
- **Known structural biases:** surface primary production is still ~4–8× low in the subtropical gyres, which
  points to missing new nitrogen supply rather than mixing or recycling. At HOT, chlorophyll at 100–200 m runs
  ~60% high (0.22 vs 0.13 mg m⁻³), and diatoms living in the deep chlorophyll maximum draw surface SiO2 down to
  0.4 vs 1.2 mmol m⁻³ observed. Hydrostation S has only five cruises per depth layer and no surviving diatoms,
  so its surface SiO2 stays too high.
- **Earlier results:** `results/params_*.csv` and `results/forward_run_misfit_by_variable.csv` are the five-site
  sets of 2026-10-08 (old initial conditions and cost), kept for reference.
- **Follow-up experiments:** `scripts/make_exp.py` (mixing, picophytoplankton growth, Fe half-saturation, an
  upwelling proxy) and `scripts/make_papa_si.py` (SiO2-only relaxation and diatom Si:C; needs `code_rbcs_mask21`).

## 7. Observations

`obs/scripts/standardize_obs.py` builds `obs/<SITE>/<SITE>_obs.csv` (time, depth, variable, value, units, source,
QC flag) from the program data at each site (HOT-DOGS and BCO-DMO, BATS/BIOS, Line P and the PMEL Papa mooring,
OceanSITES PAP, BODC PAP cruises, SOTS, JAMSTEC K2, Munida/SPOTS, OWS M, CARIACO and DYFAMED/MOOSE) plus
GLODAPv2.2023, SOCATv2026, BGC-Argo, OC-CCI v6.0 and GEOTRACES IDP2025. Sources, licenses and download steps are in
`obs/README_obs.md`; the raw data are not in the repo.

## Contact

Dustin Carroll (dustin.carroll@sjsu.edu)
