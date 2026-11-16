# Validation on real data

What was run on real GloFAS files, with logs and tables. Results are summarised in
docs/METHODS.md §5 (V1–V7). Paths in these copies are replaced by `$DATA`.

- `validate_pfe.py`, `validate_d2g.log`:
  - V1: package vs production script, WestAO grid, 2015 (bitwise).
  - V2: coast rule vs the production runoff input (40–90 N, 2015).
  - V3: dead-end candidates vs the hand-triaged list.
  - V4: raw path on WestAO vs the production forcing.
  - V5: dropout cases.
- `validate_v6.py`, `validate_v6.log`:
  - V6a: GloFAS pixel areas vs the spherical formula.
  - V6b: where the coast rule (duplicates removed) differs from the input, by region and
    cell size.
  - V6c: WestAO by 1° box.
  - V6d: production map total with each area convention.
- `validate_v7_temp.py`, `validate_v7.log`:
  - V7a: river temperature, package vs production (bitwise).
  - V7b: the climatology rebuilt from the raw river-temperature files vs the production
    temperature map.
  - V7c: raw-path temperature vs the production forcing.
- `validate_v7d_avg.py`, `validate_v7d.log`: V7d, the production map's 3×3 averaging
  emulated on the package climatology.
- `v2_outlets_2015_40N.csv.gz`: every coastal cell at 40–90 N in 2015, with mean
  discharge, sea neighbours, duplicate flag (`dup_of`) and whether the production input
  had it (`in_map`).
- `v2_only-map.csv`, `v2_only-raw.csv`: largest cells found by only one of the two.
- `v3_dead_end_candidates_40N.csv`: the 136 dead-end candidates (2015, ≥ 100 m3 s-1,
  within 250 km of the sea), before triage.
- `v4_box_differences.csv`, `v6b_boxes.csv`: box differences listed in the logs.

The scripts ran in batch jobs on our cluster. They expect our file layout (raw GloFAS
yearly files, the production runoff maps and forcing, the WestAO grid) and are included as
a record, not as tools.
