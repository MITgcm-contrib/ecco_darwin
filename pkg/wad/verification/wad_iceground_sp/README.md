# wad_iceground_sp: Sea-ice growth and brine rejection on a drying flat

`wad_iceground` under a cold atmosphere, so the ice grows and rejects brine,
with pkg/salt_plume on. Tests that the plume flux follows the WAD
dry-column ramp: brine removed from the surface and re-injected at depth
must vanish together on drying columns, so the closed basin conserves salt.

## Set-up

As `wad_iceground` (same code, grid and bathymetry), plus:

- `data.exf`: air temperature 253.15 K, no short wave, 180 W/m² long wave.
- `data.pkg`: `useSALT_PLUME=.TRUE.`; `data.salt_plume`:
  `SaltPlumeCriterion = 0.4`, `SPsalFRAC = 1.0`.
- `code/packages.conf` adds `salt_plume`. Snapshots of S and η at the start
  and end (`dumpFreq` = 10 days) for the salt budget.

## What to check

- Total salt change over 10 days: +0.006 % (reference run, arm64 gfortran).
  Without the salt_plume scaling in `WAD_DRY_FORCING` the plume adds about
  1 % of the basin's salt in 10 days (brine re-injected at depth while its
  surface removal is ramped off).
- Mean ice thickness grows: `%MON seaice_heff_mean` 0.80 m at the start, 0.840 m at day 10 (iter 14400).
- `%WAD_MON budgetErr/vol0` ~1e-15.
