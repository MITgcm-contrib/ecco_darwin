# pkg/overflood  (wad: ~/Documents/research/ECCO/wetting_drying/MITgcm)

pkg/overflood: a thin sheet of water on top of floating sea ice (river overflood at breakup), coupled to the ocean beneath and to the open (ice-free or bottom-fast-ice) cells around it. The sheet lives on the cells of ovfIceFile (floating i

**vs its upstream base (merge-base, see README):** new (not in its upstream base)

**runtime switch:** `useOVERFLOOD`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.overflood`

## Namelist parameters
### OVERFLOOD_PARM01
- `ovfIceFile` — 1 on sheet (floating ice) cells, 0 elsewhere
- `ovfHoleFile` — drainage-hole area per cell [m^2]
- `ovfFreeboard` — ice top above the ocean surface [m]
- `ovfManningN` — Manning roughness of the ice surface [s m^-1/3]
- `ovfHoleCd` — discharge coefficient of the holes
- `ovfMinDepth` — no flow on a face shallower than this [m]
- `ovfMaxFroude` — cap on the face speed, as a Froude number
- `ovfCFL` — substeps keep dt < ovfCFL*dx/sqrt(g h_max)
- `ovfOceanKeep` — an ocean cell gives the sheet at most this fraction of its water above 2*wadMinDepth per step
- `ovfMonFreq` — monitor interval [s]
- `ovfIceCh` — heat transfer from the sheet to the ice top (with pkg/seaice): Q = rho*cp*ovfIceCh*|u|*(T - Tf) with |u| the sheet speed (at least ovfIceUmin) and Tf its freezing point; Q melts snow, then ice (HEFF), limited to the sheet heat above freezing, and the
- `ovfIceUmin` — minimum sheet speed in that flux [m/s]
- `ovfUseSeaice` — the sheet follows pkg/seaice: a cell carries it where AREA >= ovfAreaMin and HEFF >= ovfHeffMin (and ovfIceFile is 1, if given), its base is the ice top: the model surface + the floe thickness HEFF/AREA under loaded ice (real fresh-water flux:
- `ovfAreaMin` — ice cover [0-1] and thickness [m] that carry a sheet (ovfUseSeaice)
- `ovfHeffMin` — ice cover [0-1] and thickness [m] that carry a sheet (ovfUseSeaice)
- `ovfLeadFrac` — with ovfUseSeaice, drainage area through leads, ovfLeadFrac*(1-AREA)*rA, added to ovfHoleFile

## Headers
- `OVERFLOOD.h` — pkg/overflood: a thin sheet of water on top of floating sea ice (river overflood at breakup), coupled to the ocean beneath and to the open (ice-free o
- `OVERFLOOD_OPTIONS.h` — CPP options file for pkg/overflood

## Routines (12)
`overflood_check.F`, `overflood_diagnostics_init.F`, `overflood_init_fixed.F`, `overflood_init_varia.F`, `overflood_read_pickup.F`, `overflood_readparms.F`, `overflood_step.F`, `overflood_tendency_apply.F`, `overflood_write_pickup.F`

## Called from outside the package
- `OVERFLOOD_TENDENCY_APPLY_S` ← `model/src/apply_forcing.F:984`
- `OVERFLOOD_TENDENCY_APPLY_T` ← `model/src/apply_forcing.F:745`
- `OVERFLOOD_STEP` ← `model/src/load_fields_driver.F:132`
- `OVERFLOOD_CHECK` ← `model/src/packages_check.F:528`
- `OVERFLOOD_INIT_FIXED` ← `model/src/packages_init_fixed.F:673`
- `OVERFLOOD_INIT_VARIA` ← `model/src/packages_init_variables.F:568`
- `OVERFLOOD_READPARMS` ← `model/src/packages_readparms.F:438`
- `OVERFLOOD_WRITE_PICKUP` ← `model/src/packages_write_pickup.F:252`

## Verification experiments compiling it (1)
`wad_overflood@wad`
