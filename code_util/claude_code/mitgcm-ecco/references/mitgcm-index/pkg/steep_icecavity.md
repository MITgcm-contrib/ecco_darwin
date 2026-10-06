# pkg/steep_icecavity

Ice cavities with steep ice draft (alternative shelfice treatment).

**runtime switch:** `useSTEEP_ICECAVITY`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.stic`
**adjoint support files:** stic_ad_diff.list

## Namelist parameters
### STIC_PARM01
- `STIClengthFile` — name of icefront length file (m/m^2) 2D file containing the ratio of the horizontal length of the ice front in each model grid cell divided by the grid cell area
- `STICdepthFile` — name of icefront depth file (m) 2D file containing depth of the ice front at each model grid cell

## CPP options (defaults as shipped)
- `ALLOW_SHITRANSCOEFF_3D` (define, STIC_OPTIONS.h) — use 3D version of transfer coefficients (needed for variable transfer coefficients)

## Headers
- `STIC.h` — BOP
- `STIC_OPTIONS.h` — BOP

## Routines (7)
`stic_check.F`, `stic_init_depths.F`, `stic_init_fixed.F`, `stic_init_varia.F`, `stic_readparms.F`, `stic_solve4fluxes.F`, `stic_thermodynamics.F`

## Called from outside the package
- `STIC_THERMODYNAMICS` ← `model/src/do_oceanic_phys.F:510`
- `STIC_INIT_DEPTHS` ← `model/src/ini_masks_etc.F:62`
- `STIC_CHECK` ← `model/src/packages_check.F:353`
- `STIC_INIT_FIXED` ← `model/src/packages_init_fixed.F:478`
- `STIC_INIT_VARIA` ← `model/src/packages_init_variables.F:402`
- `STIC_READPARMS` ← `model/src/packages_readparms.F:286`

## Verification experiments compiling it (1)
`isomip`
