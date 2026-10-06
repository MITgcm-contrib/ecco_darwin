# pkg/wad  (wadcheckin: ~/Documents/research/ECCO/MITgcm_wad_checkin)

---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+---- BOP

**vs its upstream base (merge-base, see README):** new (not in its upstream base)

**runtime switch:** `useWAD`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.wad`
**adjoint support files:** wad_check.F, wad_check_speed.F

## Namelist parameters
### WAD_PARM01
- `wadMinDepth` — thickness of the water film kept in a dry column [m]
- `wadCritDepth` — a face is open when the highest surface on either side stands more than this above the face bed crest [m] ; must be > wadMinDepth
- `wadDragDepth` — below this face depth, velocity is relaxed linearly to zero at wadMinDepth [m] ; no relaxation if <= wadMinDepth (default 0)
- `wadAdvDepth` — below this face flow depth, momentum advection (pkg/mom_vecinv) is ramped linearly to zero at wadCritDepth [m] ; off if <= wadCritDepth
- `wadGMDepth` — below this column depth, the GM/Redi tensor (pkg/gmredi) is ramped linearly to zero at wadCritDepth [m] ; off if <= wadCritDepth
- `wadKPPDepth` — below this column depth, the KPP non-local transport (KPPghat) is ramped linearly to zero at wadCritDepth [m] ; off if <= wadCritDepth (default)
- `wadKPPCapDepth` — in columns shallower than this [m], the KPP non-local flux fraction Kz*KPPghat is capped at 1 (WAD_KPP_TAPER); deeper columns keep stock KPP. <= 0: no cap
- `wadSmoothWidth` — width of tanh face-mask ramp [m] (not implemented)
- `wadMonFreq` — frequency of WAD monitor output [s]
- `wadMaxSpeed` — stop, and report where, when a velocity on an open face exceeds this [m/s] ; off if <= 0
- `wadManningN` — Manning roughness [s m^-1/3] for a depth-dependent quadratic drag Cd = g n^2/D^(1/3), limited to bottomDragQuadratic .. wadDragMax (WAD_BOTDRAG, needs ALLOW_BOTTOMDRAG_ROUGHNESS) ; off if <= 0
- `wadDragMax` — upper limit of that drag coefficient
- `wadMaxFroude` — cap the speed on each open interior face at wadMaxFroude*sqrt(g D) (D the face flow depth, as for the Manning drag), applied to u(n+1) before the outflow limiter, so volume and tracers stay conserved ; for fast thin sheets pouring off
- `wadManningFile` — optional 2-D extra Manning roughness [s m^-1/3] at tracer points (e.g. river ice cover and ice jams), added to wadManningN on each face (the larger of the two neighbours), times a factor 1 until wadManningT1, falling linearly to 0 at
- `wadManningT1`
- `wadManningT2`
- `wadConserveVol` — remove the volume added by the film floor from wet neighbours (through open faces)
- `wadCarryVel` — a face that opens takes the depth-mean u* of the upstream face of its wet cell (the advancing front then converges with resolution) ; with Nr > 1 it needs momImplVertAdv=.TRUE.
- `wadUpwindFace` — thickness of an open face from the upwind surface height (default ; if .FALSE.: the lower one with select_rStar=0, the mean with r*)
- `wadDryForcing` — no surface heat/salt/fresh-water/short-wave or pTracer surface forcing on dry columns (ramp from wadCritDepth to 2*wadCritDepth)

## CPP options (defaults as shipped)
- `WAD_DEBUG` (undef, WAD_OPTIONS.h) — Print extra per-tile information on face flips and clipped cells
- `WAD_SMOOTH_MASK` (undef, WAD_OPTIONS.h) — Smooth (tanh) face mask for adjoint differentiability (Phase 3, not implemented yet: wad_check stops if wadSmoothWidth > 0)

## Headers
- `WAD.h` — BOP
- `WAD_OPTIONS.h` — BOP

## Routines (16)
`wad_apply_floor.F`, `wad_botdrag.F`, `wad_check.F`, `wad_check_speed.F`, `wad_diagnostics_init.F`, `wad_dry_forcing.F`, `wad_gmredi_taper.F`, `wad_init_fixed.F`, `wad_init_varia.F`, `wad_kpp_taper.F`, `wad_limit_froude.F`, `wad_limit_outflow.F`, `wad_read_pickup.F`, `wad_readparms.F`, `wad_update_masks.F`, `wad_write_pickup.F`

## Called from outside the package
- `WAD_DRY_FORCING` ← `model/src/external_forcing_surf.F:404`
- `WAD_APPLY_FLOOR` ← `model/src/forward_step.F:960`
- `WAD_CHECK_SPEED` ← `model/src/forward_step.F:944`
- `WAD_LIMIT_OUTFLOW` ← `model/src/forward_step.F:943`
- `WAD_UPDATE_MASKS` ← `model/src/forward_step.F:842,862`
- `WAD_UPDATE_MASKS` ← `model/src/initialise_varia.F:309,327`
- `WAD_CHECK` ← `model/src/packages_check.F:508`
- `WAD_INIT_FIXED` ← `model/src/packages_init_fixed.F:646`
- `WAD_INIT_VARIA` ← `model/src/packages_init_variables.F:550`
- `WAD_READPARMS` ← `model/src/packages_readparms.F:423`
- `WAD_WRITE_PICKUP` ← `model/src/packages_write_pickup.F:245`
- `WAD_GMREDI_TAPER` ← `pkg/gmredi/gmredi_calc_tensor.F:1122`
- `WAD_KPP_TAPER` ← `pkg/kpp/kpp_calc.F:592`

## Verification experiments compiling it (5)
`wad_balzano@wadcheckin` `wad_estuary_3d@wadcheckin` `wad_flat_xz@wadcheckin` `wad_mudflat@wadcheckin` `wad_thacker_1d@wadcheckin`
