# pkg/regrid  (bbl: ~/Documents/research/ECCO/BBL/MITgcm)

---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+---- BOP 0

**vs its upstream base (merge-base, see README):** added: regrid_scalar_out_RL.F, regrid_scalar_out_RS.F

**runtime switch:** `useREGRID`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.regrid`
**manual:** `doc/getting_started/getting_started.rst`

## Namelist parameters
### REGRID_PARM01
- `regrid_MNC`
- `regrid_MDSIO`
- `regrid_ngrids`
- `regrid_fbname_in`
- `regrid_nout`

## CPP options (defaults as shipped)
- `RL_IS_REAL4` (undef, REGRID_OPTIONS.h) — This CPP is not set (neither def nor undef) in CPP_EEMACROS.h so it is set (to undef) here

## Headers
- `REGRID.h` — Package flags
- `REGRID_OPTIONS.h` — BOP
- `REGRID_SIZE.h` — EH3 ;;; Local Variables: *** EH3 ;;; mode:fortran *** EH3 ;;; End: ***

## Routines (7)
`regrid_check.F`, `regrid_init_fixed.F`, `regrid_init_varia.F`, `regrid_mnc_init.F`, `regrid_readparms.F`, `regrid_scalar_out_RL.F`, `regrid_scalar_out_RS.F`

## Called from outside the package
- `REGRID_CHECK` ← `model/src/packages_check.F:464`
- `REGRID_INIT_FIXED` ← `model/src/packages_init_fixed.F:600`
- `REGRID_INIT_VARIA` ← `model/src/packages_init_variables.F:491`
- `REGRID_READPARMS` ← `model/src/packages_readparms.F:379`
