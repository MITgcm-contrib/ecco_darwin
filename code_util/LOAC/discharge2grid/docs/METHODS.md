# Methods

discharge2grid turns river discharge (m³ s⁻¹ at points) into model runoff (m s⁻¹ on wet
cells), and optionally river water temperature into a flux-weighted temperature on the same
cells. Volume and heat are conserved exactly. The method grew out of the river forcing of two
MITgcm/DARWIN regional set-ups downscaled from ECCO-Darwin:
- **pan_AO**: pan-Arctic, curvilinear, 1–4 km cells.
- **WestAO**: Bering Sea, Gulf of Alaska and Beaufort; 860 × 680 cells of ~4 km, with wetting-drying.

Both use GloFAS v4 daily discharge for 1990–2026. This document explains each step, the
problems we hit with GloFAS, and how the package was checked.

## 1. Overview

```
dataset ──reader──> SourceSet ──Mapper.snap──> nearest wet cell ──Mapper.template──> footprint
 (GloFAS,           n points,                  (KD-tree on the       (Gaussian in km,
  tables,           Q(t) m3/s                   sphere; drop too      sum t·area = 1)
  fields)                                       far / outside)
                                   W (n_wet × n_source, sparse, m⁻²)
                                   runoff(t) = W · Q(t)   [m s⁻¹]
                                   temperature(t) = W · (Q T)(t) / W · Q(t)   [°C]
                                   └── writers: MITgcm binary (dense) | netCDF (receiving cells only)
```

The grid enters only through cell-centre coordinates, areas and the wet mask. The dataset
enters only through the SourceSet. W is built once per grid and source list, and reused
for every record and every year.

## 2. Sources: GloFAS, other products, temperature

### 2.1 What GloFAS provides

GloFAS v4 (`avg_dis`, m³ s⁻¹) is a daily field on a 0.05° grid from 89.975 N to 59.975 S.
It has a value on every land cell and NaN over the sea.
- The value of a cell is the flow leaving it, so along a river every cell carries roughly
  the river's discharge.
- Summing land cells therefore counts each river once per cell of its course. **Only
  outlets may be used**: the cells where a channel meets the sea.

### 2.2 Outlets

Two rules are provided. Both produce a **fixed** list of outlet cells that is then read for
every year (§2.3).

**(a) The drainage map (`outlets_from_ldd`, preferred).** GloFAS routes water with a local
drain direction map (LDD; PCRaster codes 1–9, 5 = pit) published with its static maps.
- An outlet is a land cell that drains to the sea: a pit, or a cell whose downstream
  neighbour is not land.
- Each river has exactly one, so there are no duplicates.
- Pits within `max_sea_km` of the sea are rivers that stop short of the coast (§2.5). They
  are kept and delivered at the pit; farther pits are endorheic and dropped.
- **We did not have the LDD file, so this path is tested on synthetic maps only**
  (`tests/test_ldd.py`, north-up and south-up).

**(b) The coast rule (`coastal_outlets`).** Uses the discharge field alone.
- An outlet is a land cell with at least one sea neighbour (8-connectivity), where sea means
  NaN on every scanned day. The longitude axis wraps on a global grid.
- Small coastal catchments are included.
- The rule must be paired with duplicate removal (§2.4) and the dead-end list (§2.5).

**Comparison with our production input (V2: 2015, 40–90 N).** That input was made outside
this package, by a script described as keeping "coastal cells only, duplicate coastal
outlets removed". We do not have the script; we compared its output with the raw files.
- **Cells.** The input keeps 35,005 cells carrying 243,739 m³ s⁻¹. The coast rule finds
  106,080 cells carrying 335,965 m³ s⁻¹ (+38 %).
  - 26,229 cells are in both.
  - 8,776 input cells are not coastal in the raw sea mask. In deltas and lagoons (Neva,
    Rhône, Po) the input puts the outlet 1–3 cells upstream, at the cell of largest
    discharge, where the coast rule finds 2–3 coastal cells that all carry the river (§2.4).
- **Values.** At the shared cells the input is 0.7 % below raw `avg_dis` (median ratio
  0.9928; 1st–99th percentile 0.991–0.999). This is a cell-area convention, not a
  routing difference (V6a). The input maps were made by dividing discharge by GloFAS's own
  pixel areas (`GLOFAS_pixarea_Global_03min.nc`); `glofas_to_mitgcm.py` converts back with
  a spherical formula.
  - The formula is 0.1 % smaller than the GloFAS pixel areas at 42.5 N, 0.6 % at 60 N and
    0.9 % at 85 N.
  - So **our production freshwater is 0.6–0.8 % below raw GloFAS**. For WestAO 2015 the
    runoff-map part is 36,148 m³ s⁻¹ with the formula against 36,372 with the GloFAS areas
    (V6d).
  - The raw path reads m³ s⁻¹ directly and has no area step.
  - `runoff_map_sources(..., area_file=)` uses the GloFAS areas.
- **Size classes** (coast rule, before duplicate removal):

  | coastal cells | count | share of the water |
  |---|---|---|
  | ≥ 1,000 m³ s⁻¹ | 37 | 44.5 % |
  | ≥ 100 m³ s⁻¹ | 324 | 67.9 % |
  | ≥ 10 m³ s⁻¹ | 2,260 | 84.1 % |
  | ≥ 1 m³ s⁻¹ | 15,175 | 94.6 % |
  | < 1 m³ s⁻¹ | 90,905 | 5.4 % |

- **After duplicate removal** (V6b), the coast rule is **+4.7 %** over the input: 256,549
  against 245,054 m³ s⁻¹, with the input at GloFAS pixel areas. The difference is not
  uniform:
  - **Small cells.** Cells < 10 m³ s⁻¹ carry 53,481 m³ s⁻¹ in 103,820 cells under the coast
    rule, against 27,891 in 33,323 cells in the input. The input dropped most small coastal
    cells (its duplicate removal seems to apply at every size). With the discharge field
    alone it cannot be told whether adjacent small coastal cells are separate catchments
    (the coast rule is right) or drain into each other (the input is right).
    - Small cells are 11 % of the input's water but 21 % of the coast rule's.
    - **Treat the difference, about 10 % of the total, as the method uncertainty** until the
      drainage map settles it.
  - **Large cells.** Cells ≥ 10 m³ s⁻¹ carry 203,067 m³ s⁻¹ under the coast rule against
    217,163 in the input. Several large rivers end one cell short of the raw sea mask, where
    the input's coastline still reaches them: Indigirka 1,900 m³ s⁻¹ (6.6 km from the sea),
    a James Bay river 1,244 (11 km), Alsek 827 (10 km), Kemi 574 (9 km). Under the coast
    rule these are dead ends. All four are in the candidate list (§2.5) and must be taken
    from it.
  - **Regions.** Differences by 10° × 5° box are in `validation/v6b_boxes.csv`. The
    largest are James Bay (−1,900 m³ s⁻¹) and the Indigirka box (−1,880), both dead ends,
    then the SE Alaska / BC fjord coast (+1,360), which is small cells.

### 2.3 A fixed outlet list, and why: dropouts

Our production input (the yearly "runoff maps") lost some river mouths in some years,
and the loss was intermittent:
- **Ölfusá, Iceland** (~260 m³ s⁻¹): absent in 11 of 36 years (1991–93, 2000, 2007–09,
  2012, 2021, 2024–25).
- **Pur, West Siberia** (~1,000 m³ s⁻¹): absent in 1999.
- **Thjórsá, Iceland**: absent in 2004 and 2006.
- 37 pan_AO rivers ≥ 20 m³ s⁻¹ were affected (~300 m³ s⁻¹ on average, ~0.2 % of the
  domain's water), and 21 WestAO rivers (~70 m³ s⁻¹ on average).

The water was not moved elsewhere; it was missing. We found it from a per-river atlas,
not from domain totals, which hide it.

**V5: mechanism.** We examined four cases: Ölfusá 1990/1991, Pur 1998/1999, a Chukotka
river 1991/1992 and a Gulf of Alaska river 2007/2008.
- In each, the outlet cell has an adjacent cell carrying 97–100 % of its discharge, e.g.
  Pur 1999 has 1,101.4 and 1,101 m³ s⁻¹ side by side. This is the configuration a
  duplicate-removal step has to resolve, and the same configuration exists in the years
  the river was kept.
- In every dropout year the raw discharge at the outlet cell is present: Ölfusá 292, Pur
  1,101, Chukotka 12 and Gulf of Alaska 51 m³ s⁻¹.
- So a fixed list read from the raw files keeps all four. The year-by-year conversion lost
  them.

The cure built into this package:
- Find the outlets **once**, from several years, and read the discharge at that **fixed**
  list for every year.
- A river whose last cell moves between years (the Mackenzie alternates between
  -136.075/68.925 and -136.125/68.875) then has all its cells in the list. In any year only
  one of them carries the river, so nothing is counted twice.
- `python -m discharge2grid dropouts` scans per-year source tables for rivers that vanish
  in some years. It tells a real loss (`lost_frac` ≈ 1) from a move to a nearby cell
  (≈ 0; e.g. the Gaula, Norway, which must not be filled).

For the production forcing we patched the existing input instead, with a "fill-if-missing"
step: raw discharge added at 36 listed pan_AO and 21 WestAO outlets in the years their
water was missing. A fixed list makes that patch unnecessary, and it is not part of this
package.

### 2.4 Along-coast duplicates

A channel that runs along the coast makes every cell of it coastal, and each of those
cells carries the river's discharge. `flag_coastal_duplicates` marks a coastal cell whose
8-neighbour coastal cell has a mean discharge within 10 % (both ≥ 10 m³ s⁻¹).
- The smaller twin points to the larger, which is downstream because discharge grows
  downstream. Ties are broken deterministically, so exactly one cell of a chain survives.
- **Flagged cells are removed by default** (`use=0`; `--keep-duplicates` keeps them).

V2 (2015, 40–90 N): 731 duplicate cells carrying 79,417 m³ s⁻¹; only 3 of them are in our
production input. Examples:
- the Mackenzie, two adjacent coastal cells of ~10,000 m³ s⁻¹ (-136.075/68.925 and
  -136.125/68.875);
- the Rhône, three coastal cells of 1,662 m³ s⁻¹;
- the Po, three cells of 1,488 m³ s⁻¹.

Kept, they made the WestAO total 43 % larger than production; removed, the excess is
10.6 % (V4). The 10 % / 10 m³ s⁻¹ thresholds are a judgement. Smaller chains of tiny coastal cells are not caught, which is one reason the
drainage map (§2.2a) is preferred.

### 2.5 Dead ends

Some GloFAS channels stop inland: the last cell of the river has no sea neighbour, so the
coast rule never sees it.
- The **Ob** ends 9 km from the Gulf of Ob; others include the Anadyr, Kuskokwim, Yana, Thelon
  and Mezen, and a dozen Gulf of Alaska / BC rivers (Copper, Stikine, Skeena, Taku …).
- Without them the forcing loses large rivers: the Ob alone is ~12,000 m³ s⁻¹.

`dead_end_candidates` uses a channel test:
1. Seed on local discharge maxima ≥ `min_q`.
2. The channel is the connected set of land cells carrying ≥ 80 % of the seed's discharge.
3. If no cell of that set touches a sea-adjacent cell, the river ends inland. It is
   reported with the distance to the nearest sea, if within 250 km.

The test is deliberately permissive, and the list needs hand triage:
- Discharge in GloFAS can dip slightly downstream (evaporation, floodplains), so points
  inside deltas and along big rivers are flagged too. The Yukon at -159.975/62.275 is not a
  dead end.
- Our triaged list is `examples/glofas_extra_rivers_arctic_alaska.csv` (27 rivers).
  - Each river's discharge is read at its terminal cell and delivered there, or at a given
    point: the Thelon is read inland and delivered at Chesterfield Inlet.
  - Homathko was left out because it drains to the Strait of Georgia, outside WestAO.

**V3 (2015, 40–90 N).** 136 candidates ≥ 100 m³ s⁻¹ within 250 km of the sea. **All 27
hand-triaged rivers are among them**: 24 on the same cell, the others within 2.8–13.6 km
(Arnaud, Thelon, Yana). The other 109 include:
- rivers outside our two domains;
- false alarms rejected in triage;
- rivers that are dead ends only under the raw coast rule (Indigirka, Alsek, a James Bay
  river, Kemi; §2.2), which our production input already had.

The list is therefore a complete starting point here, but much of it is not dead ends, and
**it has to be triaged in full for the raw path**.

### 2.6 Reading GloFAS fast

The raw yearly files are deflate-compressed in chunks of 29 days × 231 × 554 cells.
- Indexing `v[:, rows, cols]` at scattered points decompresses a chunk per point and per
  day: 903 s per year for 50 outlets.
- `readers.glofas.chunk_points` groups the points by chunk and reads one hyperslab per
  chunk: 9.6 s per year, with identical values.
- `scan_raw` reads whole chunks (all longitudes, one row-chunk and one time-chunk at a
  time) to build the mean field and sea mask.

### 2.7 The legacy runoff map

`readers.glofas.runoff_map_sources` reads our production input: one big-endian float32
file per year, (days, 3000, 7200), m s⁻¹ per 0.05° cell, north-to-south, non-zero only at
outlet cells.
- Every GloFAS cell in the grid's latitude band (± 1°) is a candidate source, as in the
  production code.
- The cells that can reach the grid and ever carry water become sources, in the same order.
- It exists so the production forcing can be reproduced (V1). New work should use §2.2–2.5.

### 2.8 Other products

`readers.points` reads river/station tables (CSV, or station netCDF such as Dai &
Trenberth). `readers.gridded` reads coastal runoff fields: every cell that ever carries
water becomes a source, and its discharge is the value times the cell area for flux units.

Three cautions:
- **Units.** On our system the JRA55-do river files were m s⁻¹, not the kg m⁻² s⁻¹ we expected.
- **Gauges are inland.** Move them to the mouth.
- **Routing fields.** Never feed a whole routing-model field (all land cells) to the
  gridded reader; see §2.1.

### 2.9 River temperature

Discharge products rarely carry water temperature (GloFAS has none). Temperature therefore
comes from a separate gridded product and is **sampled at the sources**
(`readers/temperature.py`):
- `TemperatureClimatology`: a daily calendar climatology, the method of our production
  forcing (below).
- `TemperatureSeries`: a dated, year-specific product, matched to the sources' dates.
- Tables: `points.from_csv(..., temp_csv=)`.
- The legacy temperature map: `runoff_map_sources(..., temp_path=)`, for exact reproduction.

**What our production forcing used.** A river-temperature product on the JRA55-do 0.25° grid
(`wattmp<year>.bin`, daily, 1991–2021, plus two climatology files).
- **Calendar slots.** Every record of every file goes to a calendar slot (Jan 1 = 0 …
  Dec 31 = 365, Feb 29 = 59 and Mar 1 = 60 in every year), so 365- and 366-day files line up.
- **Quality filter.** A value counts if it is finite, not a fill value and 0–45 °C after the
  Kelvin conversion.
- **Usable cells.** A cell is usable with ≥ 30 good values and a non-zero mean. Slots it
  lacks are filled with its own mean.
- **Lookup.** Each river outlet takes the value of the **nearest usable cell** on each day's
  slot. So the river temperature is a climatology: the same seasonal cycle every year.
- **Storage.** It was then stored on the GloFAS grid (one file per year) by giving each
  outlet cell and its 8 neighbours that value, averaging where neighbours overlap. That is
  the temperature map `glofas_to_mitgcm.py` read.

**Format traps** (each gave plausible-looking numbers when wrong):
- little-endian float32;
- rows run north to south;
- longitudes −179.875 … 179.875;
- Kelvin;
- 1e20 over the sea.

`readers.temperature.WATTMP` holds this layout. Check any new product the same way: map its
valid cells and compare known rivers.

**Unknown temperatures.** Some sources may have no temperature (rivers added by hand, sources
beyond `--temp-max-km` of a usable cell). `fill_missing` gives each one the daily temperature
of the **largest** source with known temperature within 150 km. This is the production rule
for the hand-added rivers, which took the largest real outlet within 1.5°.

## 3. Mapping (`mapping.py`)

### 3.1 Snap

All geometry uses 3-D unit vectors on the sphere, so the dateline, the poles and any grid
layout (lat-lon, curvilinear, LLC tiles) need no special handling.
- A KD-tree over the **wet** cells gives each source its nearest wet cell and the
  great-circle distance to it.
- Sources farther than `max_snap_km` are dropped. The report lists the water dropped
  within 3 × `max_snap_km`, so a too-small value is visible.

**Outside-domain test (`exclude_outside`).** A source is dropped when its nearest grid
cell, *wet or dry*, is more than `domain_tol_km` away (default 1.5 × mean cell size).
- Land cells count, so a river inside the domain whose estuary the wet mask does not
  resolve (a filled estuary, a fjord head, a smoothed delta) still counts as inside and can
  be snapped far to water.
- A river beyond the grid's edge is dropped however close it is to model water.
- This is what allows `max_snap_km` = 95 km on regional grids without importing
  neighbouring seas.

**Wet mask.** Use the mask the model treats as wet at the surface. On WestAO (wetting-drying)
the mask without floodable land keeps rivers off mudflats.

**Narrow fjords.** Where fjords are narrower than the cells, two rivers can cross land to the
same delivery cell (Kitlope and Dean, BC, at 4 km). Both are kept; the footprint shares them.

### 3.2 Footprint and conservation

Around the destination wet cell c, every wet cell k within `cutoff_sigma × sigma_km`
(3.5 σ) gets

    g_k = exp(-d_ck² / 2σ²),      t_k = g_k / Σ_j g_j A_j      so that  Σ_k t_k A_k = 1,

where d_ck is the great-circle distance (km) and A the cell area. A source of discharge Q
adds R_k = Q t_k (m s⁻¹) to cell k, and Σ_k R_k A_k = Q **exactly**.
- This holds whatever the geometry, including footprints clipped by land or the domain edge.
- Clipping does not lose water; it concentrates it.
- The participation-ratio area 1 / Σ t_k² A_k measures this. It equals 4πσ² for an
  unclipped Gaussian, and the report prints its median and minimum.
- σ is in **kilometres**, not cells, so footprints stay isotropic on curvilinear and polar
  grids.

All sources together form the sparse operator W (n_wet × n_source, m⁻²) with
`runoff = W · Q`. The rows of W that are not empty are the only cells that ever receive
water; they are the "points" of the sparse output.

### 3.3 Choosing σ

The surface cell of depth h is renewed at rate R/h, so the residence time is h/R. If it is
shorter than the advective time scale, the receiving cells simply become river water.
- Size σ on **peak** discharge, not the annual mean. Arctic rivers have peak/mean ≈ 10.
- A rough number of cells needed is N ≈ Q_peak Δt / (h A f) with f ~ 0.01 (1 % of the
  column per step).
- Too small: salinity collapses at the mouth. Too large: freshwater (and, in the full
  pipeline, nutrients and carbon) is spread onto the shelf where it does not belong.
- Our 1–4 km grids use σ = 18 km.
- On ≈ 1° (LLC90) or ⅓° (LLC270) grids, σ below the cell size reduces to "nearest wet
  cell". σ of about one cell spreads the largest rivers over a few cells.

Salinity going **negative** is a namelist problem, not a forcing one. With a fixed
`convertFW2Salt` = 35 and `useRealFreshWaterFlux` off, the salt sink does not depend on local
salinity, and a strong river mouth goes negative by construction. Use
`useRealFreshWaterFlux = .TRUE.` (with a nonlinear free surface) or `convertFW2Salt = -1`.

### 3.4 Report

`Delivery.report()` prints:
- the sources kept, too far and outside, with the water each group carries;
- the conservation ratio (must be 1 to 12 digits);
- footprint areas;
- the largest long snaps (check each lands in the right basin);
- the largest dropped sources.

`residence_days()` gives the smallest h/R on a record.

### 3.5 River temperature on the grid

Temperature is **intensive**: it is not spread like water. Each receiving cell c gets the
flux-weighted temperature of the water that reaches it:

    T_c = Σ_i q_i t_ic T_i / Σ_i q_i t_ic      so that  Σ_c R_c T_c A_c = Σ_i q_i T_i,

so the heat carried by the rivers is conserved exactly. Spreading T itself through the
Gaussian would pull it toward 0 °C at the plume edge.
- **Unknown temperature.** Water whose temperature is unknown (NaN) is left out of both
  sums. A cell whose water is all of unknown temperature gets the fill value. The report
  counts such cells and prints a heat out/in ratio, which is 1 to 12 digits only when every
  source has a temperature.
- **Fill values.** Where a cell gets no water, the MITgcm file holds 0 (as our production
  files do) and the sparse netCDF holds NaN.
- **Exactness.** The products are formed as (q t) T and summed in source order, exactly as in
  `glofas_to_mitgcm.py`, so the legacy path reproduces its temperature files bit for bit (V7).

**MITgcm.** The temperature file is the companion of the runoff file: EXF `runoftempfile`,
with `runoftemp_interpMethod = 0` and the same period as the runoff. It needs
`ALLOW_RUNOFTEMP` in `EXF_OPTIONS.h`; check how your MITgcm version uses it in the surface
heat flux.

## 4. Outputs (`writers.py`)

**MITgcm binary.** Big-endian float32, one record per time step, each record of the grid's
array shape, m s⁻¹, one file per year. With temperature, a second file
(`write_mitgcm_temperature`, °C, 0 where no water) has the same layout.
- EXF settings: `runoff_interpMethod = 0`, `useExfYearlyFields = .TRUE.`, `runoffperiod` set to
  the record length.
- With yearly files EXF reads the last record of the previous year and the first of the next
  around each new year, so write one extra year on each side of the run.
- **LLC**: build the grid from the model's own MDS output so records are in the model's
  compact layout.

**Sparse netCDF.** No zeros outside the river footprints.

| variable | dims | meaning |
|---|---|---|
| `runoff` | (time, point) | m s⁻¹, the forcing |
| `point_index` | (point) | flat C-order index into the grid array (`grid_shape` attribute) |
| `<dim>` e.g. `tile`, `j`, `i` | (point) | 0-based index along each grid dimension |
| `lon`, `lat`, `area` | (point) | cell centre and area (m²) |
| `source_lon/lat/name/kind/status/snap_km/dest_point` | (source) | what was read and what happened to it |
| `discharge` | (time, source) | m³ s⁻¹ as read (optional) |
| `link_point`, `link_source`, `link_weight` | (link) | the operator: runoff[p] = Σ w · discharge[s] |
| `temperature` | (time, point) | °C, flux-weighted; NaN where no water (only if the sources have temperature) |
| `source_temperature` | (time, source) | °C as assigned to each source (NaN = unknown) |

Global attributes record the mapping parameters. `read_sparse_netcdf(p).to_dense(k)`
rebuilds a dense record. The layout is a starting point; change it in `writers.py`.

## 5. Validation

**V0. Synthetic tests** (`tests/`):
- volume conserved to 1e-12;
- snapping, too-far and outside-domain decisions;
- unclipped footprint area within 3 % of 4πσ²;
- an LLC-like (tile, j, i) grid;
- the sparse netCDF reproduces the dense field, and the stored operator reproduces the runoff;
- on a synthetic GloFAS file: coastal outlets, along-coast duplicates, a dead end, and the
  chunk reader equal to direct slicing;
- table and gridded readers, unit conversion, and the dropout scan (real loss vs move);
- the LDD outlet finder on synthetic drainage maps, north-up and south-up: pit and non-pit
  mouths, a near-sea pit and an inland pit;
- temperature (`tests/test_temperature.py`):
  - flux-weighted mixing between two rivers, with heat conserved to 1e-12;
  - a cell reached by one river has its temperature exactly;
  - unknown temperatures, and filling from the largest neighbour within the distance;
  - the climatology builder on synthetic little-endian, north-first Kelvin files of 365 and
    366 records (calendar slots, Feb 29 from leap years only, usable cells, sampling and
    distance limit);
  - temperature tables;
  - the sparse and dense temperature outputs.

**V1. Bitwise reproduction of the production code.**
- **Synthetic input.** On the real WestAO grid with a synthetic 3-day runoff map, the package
  and `glofas_to_mitgcm.py` (σ 18 km, snap 95 km, outside-domain on) give **bitwise
  identical** files: 0 of 1,754,400 values differ.
- **Real input.** On the real input, the package's legacy reader on the
  2015 GloFAS runoff map gives **bitwise identical** runoff to `glofas_to_mitgcm.py` (same
  settings, no extra rivers): 0 of 213,452,000 values differ over 365 days. The package took
  17 s and the production script 68 s. So everything downstream of the SourceSet (snap,
  footprint, operator, output) is the method that made our forcing.

**V2. Outlet rule vs our production input** (2015, 40–90 N): §2.2 and §2.4.

**V3. Dead-end candidates vs the hand-triaged list:** §2.5. All 27 found.

**V4. Raw path on WestAO, 2015** (coast rule + the WestAO part of the hand list,
σ 18 km, snap 95 km).
- Production forcing (runoff map + 18 extra rivers; no fills that year): 49,200 m³ s⁻¹.
- Coast rule, duplicates **kept**: 70,369 m³ s⁻¹ (**+43 %**). The largest excess is the
  Mackenzie delta (+9,900 m³ s⁻¹ in three 1° boxes).
- Coast rule, duplicates **removed**: 54,393 m³ s⁻¹ (**+10.6 %**). Daily correlation with
  production 0.9988.
- Output size: 81,971 receiving cells out of 584,800 (14 %). The sparse netCDF is 137 MB
  (including the source discharge and the operator), against 854 MB for the dense binary.

**V5. Dropout mechanism:** §2.3. The raw discharge is present at the outlet cell in each
dropout year, so the fixed-list path keeps those rivers.

**V6. Where the remaining difference sits.**
- At 40–90 N: small coastal cells and short dead ends (§2.2).
- On WestAO the +10.6 % is mostly on the SE Alaska / BC fjord coast. The 1° boxes there sum
  to +6,600 m³ s⁻¹, against −1,400 elsewhere.
- The largest negative boxes are around the Alsek. The production input had it, but the
  raw run used only our 27-river hand list, which was triaged against that input and so
  does not contain it. Adding it from the candidate list gives about +12 %.
- The area convention explains −0.6 % of the production values (V6a/V6d).

**V7. River temperature** (WestAO, 2015; `validation/validate_v7.log`).
- **Legacy path: bitwise.** Reading the production runoff and temperature maps, the package
  writes temperature files identical to `glofas_to_mitgcm.py`'s: 0 of 213,452,000 values
  differ (synthetic input: 0 of 1,754,400). Heat out/in is 1.000000000000.
- **Climatology rebuilt from the raw files.** Built from the same 33 `wattmp` files as
  production (12,054 daily records), with 130,039 usable cells at 35–90 N.
  - Sampled at the 4,656 runoff-map outlets, it equals the production temperature map on
    79 % of the 1.45 million outlet-days with water. Median difference 0, 99th percentile
    2.1 °C, maximum 11.9 °C.
  - **V7d.** The cause is the production map's 3×3 averaging: each outlet passes its value
    to its 8 neighbours, and each cell takes the mean of what it receives. Adding that step
    to the package's climatology reproduces the production map on every outlet-day to
    within 0.00001 °C. 95.8 % match bitwise; the rest differ in the last float32 bit
    because the package stores the climatology in float32.
  - So **the climatology is the production one**. The package deliberately leaves out the
    averaging, which mixes the temperatures of neighbouring rivers (up to 12 °C apart at a
    few small Alaska rivers); each source keeps its own nearest-cell value.
- **Raw path.** Coast rule with duplicates removed, plus the hand list. Each outlet samples
  the climatology at the nearest usable cell (median 11 km away, at most 362 km).
  - WestAO's flux-weighted river temperature over 2015 is 7.24 °C, against 7.26 °C in
    production. The daily difference is 0.48 °C on average (at most 1.15); correlation 0.9989.
  - The raw path is warmer outside spring: Jan–Mar 1.23 vs 0.70 °C, Jul–Sep 10.87 vs 10.54,
    Oct–Dec 2.22 vs 1.79. Apr–Jun is 8.12 vs 8.07.
  - Two things differ between the paths: the extra small coastal cells (+10.6 % water, V4),
    and the dead-end rivers, which here take the temperature at their own location instead
    of borrowing a neighbour's.
  - Heat flux Σ Q T is +10.3 %, following the water.
  - A lookup of 362 km is far. `--temp-max-km` rejects such lookups, and those sources are
    then filled from a neighbouring river.

**What is and is not verified:**
- **Verified.** The mapping is exact against production. The raw reading, the fixed list
  and duplicate removal behave as described on real 2015 data. The dead-end list contains
  every river we know of.
- **Not verified.** Which small-cell treatment is right (about 10 % of the water at
  40–90 N). The LDD outlet finder on real data.

## 6. Not included, and known limits

- **Biogeochemistry.** Our production pipeline also writes 15 BGC tracer fluxes on the same
  footprints. A tracer flux is Q × C mixed on the same W, the way temperature is (§3.5), and
  could be added the same way.
- **No vertical spreading.** Runoff enters the surface cell.
- **GloFAS limits.**
  - Glacier melt is probably missing in the Gulf of Alaska.
  - Small coastal catchments are only as good as the land model.
  - Dead ends need hand triage (§2.5).
- **Delivery point.** Water enters at the nearest wet cell to the mouth. Estuaries narrower
  than the grid are not represented; the snap distance in the report shows where that matters.

## 7. Provenance

The mapping core is the water part of `glofas_to_mitgcm.py` (versions used for the pan_AO
and WestAO forcing, Oct–Nov 2026), which in turn generalised a Gulf of Mexico script from
lat-lon boxes to any grid:
- the KD-tree replaced index lookups;
- σ in km replaced σ in cells;
- a latitude band replaced lon/lat boxes, removing the dateline problem.

The river-temperature climatology repackages `make_glofas_river_temp_forcing.py` (the
script that made our temperature maps): same calendar slots, filters and nearest-cell
lookup, now sampled directly at the sources.

New in this package: the SourceSet interface, the sparse operator, the sparse netCDF
writer, N-D grids, the raw-GloFAS fixed-outlet path, the generic readers, and the
temperature sampling and filling.
