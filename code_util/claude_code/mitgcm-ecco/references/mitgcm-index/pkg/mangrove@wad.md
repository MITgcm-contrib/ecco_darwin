# pkg/mangrove  (wad: ~/Documents/research/ECCO/wetting_drying/MITgcm)

---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+---- BOP

**vs its upstream base (merge-base, see README):** new (not in its upstream base)

**runtime switch:** `useMANGROVE`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.mangrove`

## Namelist parameters
### MANGROVE_PARM01
- `mgNsp` — number of species (<= mgMaxSp)
- `mgCd` — drag coefficient of each species [-]
- `mgTrunkN` — trunk density [stems/m^2]
- `mgTrunkD` — trunk diameter [m]
- `mgRootA0` — root frontal area per unit volume at the bed [1/m]
- `mgRootH` — root height [m]
- `mgRootP` — root profile exponent [-]
- `mgRootD` — root (or pneumatophore) diameter [m], for the stem-wake mixing
- `mgCoverFile` — cover fraction (0-1) of each species [2-D file]
- `mgNepfAlpha` — stem-wake eddy viscosity nu = alpha*(Cd a d)^(1/3) *|u|*d (Nepf 1999), added to the vertical viscosity and diffusivity ; 0: off
- `mgMaxKz` — upper limit of the stem-wake viscosity [m^2/s]
- `mgBedShelter` — bed shear stress (pkg/sediment) times 1 - mgBedShelter*total cover ; 0: off
- `mgWaveAtt` — wave height (pkg/sediment wave modes 2, 3) decays through the forest, H = H0/(1 + beta x) (Mendez and Losada 2004, emergent), x = mgDistFile
- `mgDistFile` — distance into the forest from its seaward edge [m]

## Headers
- `MANGROVE.h` — BOP
- `MANGROVE_OPTIONS.h` — BOP

## Routines (10)
`mangrove_calc_visc.F`, `mangrove_check.F`, `mangrove_diagnostics_init.F`, `mangrove_drag.F`, `mangrove_init_fixed.F`, `mangrove_init_varia.F`, `mangrove_readparms.F`, `mangrove_wave_att.F`

## Called from outside the package
- `MANGROVE_CALC_DIFF` ← `model/src/calc_3d_diffusivity.F:249`
- `MANGROVE_CALC_VISC` ← `model/src/calc_viscosity.F:124`
- `MANGROVE_DRAG` ← `model/src/forward_step.F:943`
- `MANGROVE_CHECK` ← `model/src/packages_check.F:521`
- `MANGROVE_INIT_FIXED` ← `model/src/packages_init_fixed.F:664`
- `MANGROVE_INIT_VARIA` ← `model/src/packages_init_variables.F:562`
- `MANGROVE_READPARMS` ← `model/src/packages_readparms.F:433`
- `MANGROVE_WAVE_ATT` ← `pkg/sediment/sediment_waves.F:181`

## Verification experiments compiling it (1)
`wad_mangrove@wad`
