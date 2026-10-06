# pkg/gridalt

Alternative vertical grid support for atmospheric physics (fizhi).

**runtime switch:** `useGRIDALT`-style flag in `data.pkg` (check exact name in packages_boot.F)
**manual:** `doc/phys_pkgs/gridalt.rst`, `doc/outp_pkgs/outp_pkgs.rst`

## Headers
- `GRIDALT_OPTIONS.h` — Package-specific options go here
- `gridalt_mapping.h` — Alternate grid Mapping Common

## Routines (6)
`dyn2phys.F`, `gridalt_diagnostics_init.F`, `gridalt_initialise.F`, `gridalt_update.F`, `make_phys_grid.F`, `phys2dyn.F`

## Called from outside the package
- `GRIDALT_UPDATE` ← `model/src/forward_step.F:1117`
- `GRIDALT_UPDATE` ← `model/src/initialise_varia.F:367`
- `GRIDALT_INITIALISE` ← `model/src/packages_init_fixed.F:591`
- `DYN2PHYS` ← `pkg/fizhi/fizhi_init_vars.F:159,170,179,188`
- `PHYS2DYN` ← `pkg/fizhi/fizhi_wrapper.F:291,300,309,318`
- `DYN2PHYS` ← `pkg/fizhi/step_fizhi_corr.F:213,223,232,241`
- `PHYS2DYN` ← `pkg/fizhi/step_fizhi_corr.F:143,152,161,170`

## Verification experiments compiling it (3)
`fizhi-cs-32x32x40` `fizhi-cs-aqualev20` `fizhi-gridalt-hs`
