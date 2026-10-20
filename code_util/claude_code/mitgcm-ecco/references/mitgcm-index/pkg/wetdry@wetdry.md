# pkg/wetdry  (wetdry: ~/Documents/research/ECCO/wetting_drying/MITgcm)

---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+---- BOP

**vs its upstream base (merge-base, see README):** new (not in its upstream base)

**runtime switch:** `useWETDRY`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.wetdry`

## Namelist parameters
### WETDRY_PARM01
- `wdMinDepth` — thickness of the water film kept in a dry column [m]
- `wdCritDepth` — a face is open when the highest surface on either side stands more than this above the face bed crest [m] ; must be > wdMinDepth
- `wdDragDepth` — below this face depth, velocity is relaxed linearly to zero at wdMinDepth [m] ; no relaxation if <= wdMinDepth (default 0)
- `wdAdvDepth` — below this face flow depth, momentum advection (pkg/mom_vecinv) is ramped linearly to zero at wdCritDepth [m] ; off if <= wdCritDepth
- `wdGMDepth` — below this column depth, the GM/Redi tensor (pkg/gmredi) is ramped linearly to zero at wdCritDepth [m] ; off if <= wdCritDepth
- `wdIceDepth` — below this face flow depth, the sea-ice surface tilt force (pkg/seaice) is ramped linearly to zero at wdCritDepth [m] ; off if <= wdCritDepth
- `wdLoadDepth` — with the sea-ice load on the surface (real fresh- water flux), the load a column passes to the ocean (and to the ice tilt) is ramped from 0 at wdCritDepth of water to the full load at this depth: the bed carries grounded ice. Loaded ice
- `wdLoadTau` — wdLoadFac relaxes to its depth ramp with this time scale [s] (<= 0: at once). The load reacts to the water depth explicitly with a gain of draft/(wdLoadDepth - wdCritDepth); without a lag this is unstable at ordinary time steps
- `wdIceMaxSlip` — pkg/seaice: after the ice solve, each ice velocity component is kept within wdIceMaxSlip of the surface current at its point. Floating broken ice cannot outrun the water carrying it; at the edge of grounded (masked) fast ice the solver shot partly
- `wdIceBondH` — with wdIceSeal: ice once sealed to the bed stays frozen to it (bonded) until its thickness (HEFF, m) falls below this, whatever the water beneath does. Without it the seal follows the film under the ice, which grows with the top meltwater
- `wdBondHead` — with wdIceBondH: a bonded column lets go when the water on an open wet neighbour stands more than this above its ice top (floe top, HEFF/AREA over the column surface). Flooded bottom-fast ice is always buoyant once water gets beneath it; only
- `wdBondTau` — with wdIceBondH: a released bond (wdBond) falls from 1 to 0 over this time [s]; the faces of the column open at once but the thin-water velocity factor (wdRelaxW/S) is scaled by 1 - wdBond, so the flow through them ramps up instead of a dam
- `wdIceMinHeff` — sea-ice faces whose two cells hold less ice than this (HEFF sum, m) are left out of the ice solver (default 0: only ice-free faces)
- `wdKPPDepth` — below this column depth, the KPP non-local transport (KPPghat) is ramped linearly to zero at wdCritDepth [m] ; off if <= wdCritDepth (default)
- `wdKPPCapDepth` — in columns shallower than this [m], the KPP non-local flux fraction Kz*KPPghat is capped at 1 (WETDRY_KPP_TAPER); deeper columns keep stock KPP. <= 0: no cap
- `wdSmoothWidth` — width of tanh face-mask ramp [m] (not implemented)
- `wdMonFreq` — frequency of wetdry monitor output [s]
- `wdMaxSpeed` — stop, and report where, when a velocity on an open face exceeds this [m/s] ; off if <= 0
- `wdManningN` — Manning roughness [s m^-1/3] for a depth-dependent quadratic drag Cd = g n^2/D^(1/3), limited to bottomDragQuadratic .. wdDragMax (WETDRY_BOTDRAG, needs ALLOW_BOTTOMDRAG_ROUGHNESS) ; off if <= 0
- `wdDragMax` — upper limit of that drag coefficient
- `wdMaxFroude` — cap the speed on each open interior face at wdMaxFroude*sqrt(g D) (D the face flow depth, as for the Manning drag), applied to u(n+1) before the outflow limiter, so volume and tracers stay conserved ; for fast thin sheets pouring off
- `wdManningFile` — optional 2-D extra Manning roughness [s m^-1/3] at tracer points (e.g. river ice cover and ice jams), added to wdManningN on each face (the larger of the two neighbours), times a factor 1 until wdManningT1, falling linearly to 0 at
- `wdManningT1`
- `wdManningT2`
- `wdConserveVol` — remove the volume added by the film floor from wet neighbours (through open faces)
- `wdCarryVel` — a face that opens takes the depth-mean u* of the upstream face of its wet cell (the advancing front then converges with resolution) ; with Nr > 1 it needs momImplVertAdv=.TRUE.
- `wdIceGround` — pkg/seaice: ice does not move across a face where its draft exceeds the flow depth (grounded ice)
- `wdUpwindFace` — thickness of an open face from the upwind surface height (default ; if .FALSE.: the lower one with select_rStar=0, the mean with r*)
- `wdDryForcing` — no surface heat/salt/fresh-water/short-wave or pTracer surface forcing on dry columns (ramp from wdCritDepth to 2*wdCritDepth)
- `wdIceTopMelt` — pkg/seaice: ice on a dry column (grounded ice) may still melt (atmospheric melt, meltwater into the film); growth there stays off (WETDRY_SEAICE_DRY)
- `wdIceSeal` — bottom-fast ice: with the load, ice with less than max(wdCritDepth, wdLoadDepth) of water under it is frozen to the bed and closes the column's faces (floods go over it, not under), until it melts
- `wdIceOBpack` — keep ice-free (or wdIceMinHeff) faces within 2 cells of an open boundary in the ice solver, so that an offshore pack can drift in (default .TRUE.). .FALSE. where no pack enters (ice-free leads along the boundaries): trace ice there ran away at 1 m/s

## CPP options (defaults as shipped)
- `WETDRY_DEBUG` (undef, WETDRY_OPTIONS.h) — Print extra per-tile information on face flips and clipped cells
- `WETDRY_SMOOTH_MASK` (undef, WETDRY_OPTIONS.h) — Smooth (tanh) face mask for adjoint differentiability (Phase 3, not implemented yet: wetdry_check stops if wdSmoothWidth > 0)

## Headers
- `WETDRY.h` — BOP
- `WETDRY_OPTIONS.h` — BOP

## Routines (18)
`wetdry_apply_floor.F`, `wetdry_botdrag.F`, `wetdry_check.F`, `wetdry_check_speed.F`, `wetdry_diagnostics_init.F`, `wetdry_dry_forcing.F`, `wetdry_gmredi_taper.F`, `wetdry_init_fixed.F`, `wetdry_init_varia.F`, `wetdry_kpp_taper.F`, `wetdry_limit_froude.F`, `wetdry_limit_outflow.F`, `wetdry_read_pickup.F`, `wetdry_readparms.F`, `wetdry_seaice_dry.F`, `wetdry_seaice_freeslip.F`, `wetdry_update_masks.F`, `wetdry_write_pickup.F`

## Called from outside the package
- `WETDRY_DRY_FORCING` ← `model/src/external_forcing_surf.F:444`
- `WETDRY_APPLY_FLOOR` ← `model/src/forward_step.F:969`
- `WETDRY_CHECK_SPEED` ← `model/src/forward_step.F:953`
- `WETDRY_LIMIT_OUTFLOW` ← `model/src/forward_step.F:952`
- `WETDRY_UPDATE_MASKS` ← `model/src/forward_step.F:842,862`
- `WETDRY_UPDATE_MASKS` ← `model/src/initialise_varia.F:309,327`
- `WETDRY_CHECK` ← `model/src/packages_check.F:508`
- `WETDRY_INIT_FIXED` ← `model/src/packages_init_fixed.F:646`
- `WETDRY_INIT_VARIA` ← `model/src/packages_init_variables.F:550`
- `WETDRY_READPARMS` ← `model/src/packages_readparms.F:423`
- `WETDRY_WRITE_PICKUP` ← `model/src/packages_write_pickup.F:245`
- `WETDRY_GMREDI_TAPER` ← `pkg/gmredi/gmredi_calc_tensor.F:1122`
- `WETDRY_KPP_TAPER` ← `pkg/kpp/kpp_calc.F:592`
- `WETDRY_SEAICE_FREESLIP` ← `pkg/seaice/seaice_calc_viscosities.F:482`
- `WETDRY_SEAICE_DRY` ← `pkg/seaice/seaice_model.F:279,283`

## Verification experiments compiling it (8)
`wetdry_balzano@wetdry` `wetdry_estuary_3d@wetdry` `wetdry_flat_xz@wetdry` `wetdry_iceground@wetdry` `wetdry_mangrove@wetdry` `wetdry_mudflat@wetdry` `wetdry_overflood@wetdry` `wetdry_thacker_1d@wetdry`
