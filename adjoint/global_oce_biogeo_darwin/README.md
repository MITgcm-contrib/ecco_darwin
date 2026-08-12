# global_oce_biogeo_darwin

Forward-only Darwin ecosystem model, configured to run on **exactly** the
same grid, physics, and physical forcing as the sibling
`global_oce_biogeo_bling_SOCAT` experiment -- same bathymetry, wind stress,
SST/SSS restoring, sea-ice climatology, and circulation -- with only the
biogeochemistry package swapped (Darwin's trait-based, multi-species
ecosystem model instead of BLING's simpler implicit ecosystem). No
adjoint, no optimization -- this experiment exists purely as the "other
model" half of the pilot study in `../darwin_bling_comparison/`, which
compares Darwin's forward finite-difference sensitivities against BLING's
adjoint sensitivities at matched grid points.

See `../darwin_bling_comparison/README.md` for the full comparison
methodology, results, and figures.

## Get the code

```sh
git clone https://github.com/MITgcm-contrib/ecco_darwin.git
git clone https://github.com/MITgcm/MITgcm.git
cd MITgcm
git clone https://github.com/darwinproject/darwin3.git   # Darwin ecosystem model
```

## Directory layout

```
code_ad/        MITgcm source-code customizations (headers, package list)
input_ad/       Runtime namelists (binary/NetCDF input files are NOT included here -- see below)
scripts/        IC setup: regrid_darwin_ics.py, match_ics_to_bling.py
```

**This repo intentionally excludes large binary/NetCDF files** to keep
the repo small -- only the text namelists (`data*`) are checked in.

## The ecosystem configuration: 6+4+0

Darwin here uses a 10-plankton-type configuration: 6 phototrophic
(phytoplankton) types + 4 non-phototrophic (zooplankton) types + 0
explicit bacteria (DOM remineralization handled implicitly), i.e.
`nPhoto=6`, `nplank=10`, `nGroup=7` functional groups (`code_ad/DARWIN_SIZE.h`).
This is a much richer community-ecology model than BLING's 2-3 implicit
biomass pools with Liebig (single-limiting-nutrient) growth -- see
`../darwin_bling_comparison/README.md` for why that structural difference
matters for interpreting the comparison results.

`input_ad/data.traits` (the trait/palatability parameters for all 10
types) and the matching `code_ad/DARWIN_SIZE.h` dimensions originate from
a reference "6+4+0" Darwin configuration; obtain a compatible
`data.traits` from a Darwin reference run with the same `DARWIN_SIZE.h`
dimensions (`nplank=10, nGroup=7, nopt=12, nPhoto=6`).

## Reproducing the initial conditions

Darwin's biogeochemical tracer ICs were originally regridded from a
separate, independent source (a prior Darwin LLC270 ecosystem run,
`scripts/regrid_darwin_ics.py`) -- this is **not** cross-validated against
BLING's own ICs, and for the 6 tracers Darwin and BLING have in common
(DIC, ALK, O2, NO3, PO4, Fe/FeT) that mismatch was found to explain most
of a large adjoint-vs-forward-FD discrepancy in the comparison study (see
`../darwin_bling_comparison/README.md`). **Run `scripts/match_ics_to_bling.py`
after obtaining `global_oce_biogeo_bling_SOCAT/input_ad`'s tracer IC
files** to replace Darwin's independently-sourced ICs for those 6 tracers
with BLING's own (unit-converted mol/m3 -> mmol/m3) before running any
comparison. Darwin's other ~30 tracers (NO2, NH4, SiO2, DOC/DON/DOP/DOFe,
POC/PON/POP/POFe/POSi, PIC, CDOM, per-functional-type biomass and Chl
pools) have no BLING equivalent and keep their LLC270-regridded values.

Both scripts read/write via environment variables rather than hardcoded
paths -- see each script's docstring.

## Physics/forcing setup notes

Getting Darwin to run identically to BLING's physics required matching
several settings that are easy to get subtly wrong:

- `useEXF=.TRUE.` is required by `darwin_check.F` even though BLING's own
  `EXTERNAL_FIELDS_LOAD` still does the actual SST/SSS forcing
  afterward and overwrites whatever EXF computed -- confirmed via source
  inspection that this doesn't change the physics, it's just a
  required-package check.
- `tauThetaClimRelax`/`tauSaltClimRelax` must live in `data.exf` as
  `climsstTauRelax`/`climsssTauRelax` once `useEXF=.TRUE.` (not in
  `data`), with valid `climsststartdate1/2` (EXF's `CAL_FULLDATE` will
  abort on the default `0,0` if `useCAL=.TRUE.`).
- **Sea-ice forcing period**: BLING's `fice.bin` is a 12-record monthly
  climatology, read with `iceperiod=-12.` (the MITgcm convention for "N
  repeating monthly records"). Darwin's `data.darwin` defaults to
  `iceperiod=86400.0` (a literal daily period) -- silently wrong for a
  monthly file, and MITgcm won't catch this until it runs off the end of
  the 12 records. Set `iceperiod=-12.` to match.
- `mah_flux_smooth.bin` is referenced by `data.darwin` but isn't part of
  Darwin's own standard input set -- copy it from
  `global_oce_biogeo_bling_SOCAT/input_ad/`.
- `data.diagnostics`: remove any adjoint-only diagnostic fields (e.g.
  `ADJptr01..08`) if reusing a BLING-derived `data.diagnostics` as a
  starting point -- they don't exist in a forward-only (non-TAF) build
  and `DIAGNOSTICS_SET_POINTERS` will abort.

## Known limitation: pkg/profiles' PCO/PH/CHL/POC support

If you want Darwin's *own* `pkg/cost`/`pkg/profiles` machinery to compute
a SOCAT-pCO2 cost function directly (rather than replicating the formula
externally in Python, as `../darwin_bling_comparison/scripts/compute_cost_fd.py`
does), note that the generic `pkg/profiles/profiles_interp.F` has no case
for the `'PCO'`/`'PH'`/`'CHL'`/`'POC'` variable names -- it silently
returns 0 for all of them. BLING's own experiment gets working PCO
support from an experiment-specific override,
`darwin3/verification/global_oce_biogeo_bling/code_ad/profiles_interp.F`
(adds `pCO2`/`pH`/`CHL`/`POC` array lookups), which was never copied into
this Darwin experiment's `code_ad/`. This didn't block the pilot study
(which computes the cost function independently in Python from Darwin's
raw `pCO2` diagnostic), but would need fixing -- copy that override file
into `code_ad/` -- for Darwin's own internal cost/profiles output to be
meaningful.
