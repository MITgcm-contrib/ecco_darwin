# pkg/atm_phys

Grey-radiation idealized atmospheric physics (Frierson/O'Gorman-type moist physics).

**runtime switch:** `useATM_PHYS`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.atm_gray`, `data.atm_phys`
**manual:** `doc/examples/examples.rst`

## Namelist parameters
### ATMOSPHERE_NML
- `turb`
- `ldry_convection`
- `lwet_convection`
- `do_virtual`
- `two_stream`
- `mixed_layer_bc`
- `roughness_heat`
- `roughness_moist`
- `roughness_mom`
### ATM_PHYS_PARM01
- `atmPhys_addTendT` — apply ATM_PHYS tendency to temperature
- `atmPhys_addTendS` — apply ATM_PHYS tendency to Specific Humid
- `atmPhys_addTendU` — apply ATM_PHYS tendency to U-component wind
- `atmPhys_addTendV` — apply ATM_PHYS tendency to V-component wind
- `atmPhys_tauDampUV` — damping time-scale (s)
- `atmPhys_dampUVfac` — damping coefficient for each level
- `atmPhys_stepSST` — step forward SST
- `atmPhys_sstFile` — name of initial SST [in K] file
- `atmPhys_qFlxFile` — name of Q-flux file
- `atmPhys_mxldFile` — name of Mixed-Layer Depth file
- `atmPhys_albedoFile` — name of Albedo file
- `atmPhys_ozoneFile` — name of annual mean ozone concentration file (units: mol/mol i.e. volume mixing ratio)
### MY25_TURB_NML
### SHALLOW_CONV_NML
### VERT_TURB_DRIVER_NML
- `do_shallow_conv`
- `do_mellor_yamada`

## Headers
- `ATM_PHYS_OPTIONS.h` — Package-specific options go here
- `ATM_PHYS_PARAMS.h` — --   ATM_PHYS parameters atmPhys_addTendT :: apply ATM_PHYS tendency to temperature atmPhys_addTendS :: apply ATM_PHYS tendency to Specific Humid atmP
- `ATM_PHYS_VARS.h` — -    AtmPhys 2-dim. fields

## Routines (36)
`atm_phys_check.F`, `atm_phys_diagnostics_init.F`, `atm_phys_driver.F`, `atm_phys_dyn2phys.F`, `atm_phys_init_fixed.F`, `atm_phys_init_varia.F`, `atm_phys_read_pickup.F`, `atm_phys_readparms.F`, `atm_phys_tendency_apply.F`, `atm_phys_write_pickup.F`, `dargan_bettsmiller_mod.F90`, `lscale_cond_mod.F90`, `my25_turb_mod.F90`, `shallow_conv_mod.F90`, `simple_sat_vapor_pres_mod.F90`

## Called from outside the package
- `ATM_PHYS_TENDENCY_APPLY_S` ← `model/src/apply_forcing.F:861`
- `ATM_PHYS_TENDENCY_APPLY_T` ← `model/src/apply_forcing.F:491`
- `ATM_PHYS_TENDENCY_APPLY_U` ← `model/src/apply_forcing.F:113`
- `ATM_PHYS_TENDENCY_APPLY_V` ← `model/src/apply_forcing.F:303`
- `ATM_PHYS_DRIVER` ← `model/src/do_atmospheric_phys.F:134`
- `ATM_PHYS_TENDENCY_APPLY_S` ← `model/src/external_forcing.F:691`
- `ATM_PHYS_TENDENCY_APPLY_T` ← `model/src/external_forcing.F:369`
- `ATM_PHYS_TENDENCY_APPLY_U` ← `model/src/external_forcing.F:76`
- `ATM_PHYS_TENDENCY_APPLY_V` ← `model/src/external_forcing.F:216`
- `ATM_PHYS_CHECK` ← `model/src/packages_check.F:400`
- `ATM_PHYS_INIT_FIXED` ← `model/src/packages_init_fixed.F:571`
- `ATM_PHYS_INIT_VARIA` ← `model/src/packages_init_variables.F:481`
- `ATM_PHYS_READPARMS` ← `model/src/packages_readparms.F:363`
- `ATM_PHYS_WRITE_PICKUP` ← `model/src/packages_write_pickup.F:207`

## Verification experiments compiling it (1)
`atm_gray`
