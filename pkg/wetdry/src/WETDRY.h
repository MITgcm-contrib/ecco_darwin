CBOP
C     !ROUTINE: WETDRY.h
C     !INTERFACE:
C     #include "WETDRY.h"

C     !DESCRIPTION:
C     *================================================================*
C     | WETDRY.h
C     | o Header file for pkg "wetdry" (wetting and drying): thin-film
C     |   scheme with dynamic face masking, for r* non-linear
C     |   free surface.
C     *================================================================*
CEOP

C--   wetdry parameters (namelist WETDRY_PARM01 in data.wetdry)
C     wdMinDepth    :: thickness of the water film kept in a dry column [m]
C     wdCritDepth   :: a face is open when the highest surface on either
C                       side stands more than this above the face bed
C                       crest [m] ; must be > wdMinDepth
C     wdDragDepth   :: below this face depth, velocity is relaxed
C                       linearly to zero at wdMinDepth [m] ;
C                       no relaxation if <= wdMinDepth (default 0)
C     wdAdvDepth    :: below this face flow depth, momentum advection
C                       (pkg/mom_vecinv) is ramped linearly to zero at
C                       wdCritDepth [m] ; off if <= wdCritDepth
C     wdGMDepth     :: below this column depth, the GM/Redi tensor
C                       (pkg/gmredi) is ramped linearly to zero at
C                       wdCritDepth [m] ; off if <= wdCritDepth
C     wdKPPDepth    :: below this column depth, the KPP non-local
C                       transport (KPPghat) is ramped linearly to zero at
C                       wdCritDepth [m] ; off if <= wdCritDepth
C                       (default)
C     wdKPPCapDepth :: in columns shallower than this [m], the KPP
C                       non-local flux fraction Kz*KPPghat is capped at 1
C                       (WETDRY_KPP_TAPER); deeper columns keep stock KPP.
C                       <= 0: no cap
C     wdIceDepth    :: below this face flow depth, the sea-ice surface
C                       tilt force (pkg/seaice) is ramped linearly to
C                       zero at wdCritDepth [m] ; off if <= wdCritDepth
C     wdLoadDepth   :: with the sea-ice load on the surface (real fresh-
C                       water flux), the load a column passes to the
C                       ocean (and to the ice tilt) is ramped from 0 at
C                       wdCritDepth of water to the full load at this
C                       depth: the bed carries grounded ice. Loaded ice
C                       is grounded where less than this depth of water
C                       is under it. Off (full load, grounded at
C                       wdCritDepth) if <= wdCritDepth [m]
C     wdLoadTau     :: wdLoadFac relaxes to its depth ramp with this
C                       time scale [s] (<= 0: at once). The load reacts
C                       to the water depth explicitly with a gain of
C                       draft/(wdLoadDepth - wdCritDepth); without a lag
C                       this is unstable at ordinary time steps
C     wdIceMaxSlip  :: pkg/seaice: after the ice solve, each ice velocity
C                       component is kept within wdIceMaxSlip of the
C                       surface current at its point. Floating broken ice
C                       cannot outrun the water carrying it; at the edge of
C                       grounded (masked) fast ice the solver shot partly
C                       covered cells out at 2-10 m/s in a 2 m/s flood
C                       (Colville freshet), which blew up the free surface.
C                       Off if <= 0 [m/s]
C     wdIceSeal     :: bottom-fast ice: with the load, ice with less than
C                       max(wdCritDepth, wdLoadDepth) of water under it
C                       is frozen to the bed and closes the column's faces
C                       (floods go over it, not under), until it melts
C     wdIceBondH    :: with wdIceSeal: ice once sealed to the bed stays
C                       frozen to it (bonded) until its thickness (HEFF,
C                       m) falls below this, whatever the water beneath
C                       does. Without it the seal follows the film under
C                       the ice, which grows with the top meltwater
C                       (wdIceTopMelt), so 1.7 m bottom-fast ice let the
C                       flood under it and floated off (Colville freshet).
C                       Off if <= 0 [m]
C     wdBondHead    :: with wdIceBondH: a bonded column lets go when the
C                       water on an open wet neighbour stands more than
C                       this above its ice top (floe top, HEFF/AREA over
C                       the column surface). Flooded bottom-fast ice is
C                       always buoyant once water gets beneath it; only
C                       the bond holds it, and a thin remnant cannot hold
C                       back metres of water (Colville: 0.25 m of ice held
C                       a 2.5 m step until HEFF < wdIceBondH, then the
C                       dam broke at 5-7 m/s). Off if <= 0 [m]
C     wdBondLag     :: with wdBondHead: a bonded column lets go only after
C                       the head condition has held for this long in
C                       total [s] (wdLiftT). Bottom-fast ice under the
C                       overflood is held to the bed by suction and stays
C                       submerged for ~1-2 weeks; water seeps beneath over
C                       "days rather than hours" (Reimnitz 2002, Hearon et
C                       al. 2009). The HEFF < wdIceBondH release stays
C                       immediate. Off if <= 0 (default)
C     wdBondTau     :: with wdIceBondH: a released bond (wdBond) falls
C                       from 1 to 0 over this time [s]; the faces of the
C                       column open at once but the thin-water velocity
C                       factor (wdRelaxW/S) is scaled by 1 - wdBond, so
C                       the flow through them ramps up instead of a dam
C                       break. At once if <= 0 (default)
C     wdIceMinHeff  :: sea-ice faces whose two cells hold less ice than
C                       this (HEFF sum, m) are left out of the ice solver
C                       (default 0: only ice-free faces)
C     wdIceOBpack   :: keep ice-free (or wdIceMinHeff) faces within 2
C                       cells of an open boundary in the ice solver, so
C                       that an offshore pack can drift in (default .TRUE.).
C                       .FALSE. where no pack enters (ice-free leads along
C                       the boundaries): trace ice there ran away at 1 m/s
C                       in the Colville cycle run (2026-10-04)
C     wdSmoothWidth :: width of tanh face-mask ramp [m] (not implemented)
C     wdMonFreq     :: frequency of wetdry monitor output [s]
C     wdMaxSpeed    :: stop, and report where, when a velocity on an
C                       open face exceeds this [m/s] ; off if <= 0
C     wdConserveVol :: remove the volume added by the film floor from
C                       wet neighbours (through open faces)
C     wdUpwindFace  :: thickness of an open face from the
C                       upwind surface height (default ; if .FALSE.: the
C                       lower one with select_rStar=0, the mean with r*)
C     wdDryForcing  :: no surface heat/salt/fresh-water/short-wave or
C                       pTracer surface forcing on dry columns (ramp
C                       from wdCritDepth to 2*wdCritDepth)
C     wdIceTopMelt  :: pkg/seaice: ice on a dry column (grounded ice) may
C                       still melt (atmospheric melt, meltwater into the
C                       film); growth there stays off (WETDRY_SEAICE_DRY)
C     wdTopMeltSheet:: with wdIceTopMelt, wdIceBondH and pkg/overflood: the
C                       top melt of ice frozen to the bed (wdBond = 1)
C                       stays on the ice (joins the overflood sheet)
C                       instead of draining into the water beneath it.
C                       Bottom-fast ice is frozen to the bed, so its melt
C                       ponds on top; draining it beneath floated the
C                       near-shore ice before the Colville freshet (the
C                       water beneath rose 0.23 -> 0.48 m, days 0-9).
C                       Default .FALSE.
C     wdTopMeltFloat:: with wdTopMeltSheet: floating ice keeps its top melt
C                       too (melt ponds on all landfast ice, not only on
C                       the bonded band, whose edge otherwise shows as a
C                       line of meltwater). Default .FALSE.
C     wdPondMax     :: with wdTopMeltSheet: the kept melt fills the sheet on
C                       a bonded column only to this depth [m]; the rest
C                       drains beneath as without wdTopMeltSheet (ponds
C                       drain sideways to cracks). No cap if <= 0 (default)
C     wdBondFormEnd :: with wdIceBondH: no new bond (bottom-fast ice) forms
C                       after this model time [s]; bottom-fast ice forms
C                       in winter, while under the freshet sheet a loaded
C                       column pushed down to the bed sealed and bonded
C                       over 3-8 m beds (Colville test, runaway). Default
C                       1e30 (no limit)
C     wdIceGround   :: pkg/seaice: ice does not move across a face where
C                       its draft exceeds the flow depth (grounded ice)
C     wdManningN    :: Manning roughness [s m^-1/3] for a depth-dependent
C                       quadratic drag Cd = g n^2/D^(1/3), limited to
C                       bottomDragQuadratic .. wdDragMax (WETDRY_BOTDRAG,
C                       needs ALLOW_BOTTOMDRAG_ROUGHNESS) ; off if <= 0
C     wdDragMax     :: upper limit of that drag coefficient
C     wdManningFile :: optional 2-D extra Manning roughness [s m^-1/3] at
C                       tracer points (e.g. river ice cover and ice jams),
C                       added to wdManningN on each face (the larger of
C                       the two neighbours), times a factor 1 until
C                       wdManningT1, falling linearly to 0 at
C                       wdManningT2 [s, model time] (constant if
C                       wdManningT2 <= wdManningT1)
C     wdMaxFroude   :: cap the speed on each open interior face at
C                       wdMaxFroude*sqrt(g D) (D the face flow depth, as
C                       for the Manning drag), applied to u(n+1) before
C                       the outflow limiter, so volume and tracers stay
C                       conserved ; for fast thin sheets pouring off
C                       banks (breakup floods) ; off if <= 0
C     wdCarryVel    :: a face that opens takes the depth-mean u* of the
C                       upstream face of its wet cell (the advancing
C                       front then converges with resolution) ; with
C                       Nr > 1 it needs momImplVertAdv=.TRUE.
      _RL wdMinDepth
      _RL wdCritDepth
      _RL wdDragDepth
      _RL wdAdvDepth
      _RL wdGMDepth
      _RL wdIceDepth, wdLoadDepth, wdLoadTau, wdIceMinHeff
      _RL wdIceMaxSlip, wdIceBondH, wdBondHead, wdBondTau, wdBondLag
      _RL wdPondMax, wdBondFormEnd
      _RL wdKPPDepth, wdKPPCapDepth
      _RL wdSmoothWidth
      _RL wdMonFreq
      _RL wdMaxSpeed
      _RL wdManningN, wdDragMax, wdMaxFroude
      _RL wdManningT1, wdManningT2
      CHARACTER*(MAX_LEN_FNAM) wdManningFile
      LOGICAL wdConserveVol, wdCarryVel, wdIceGround, wdUpwindFace
      LOGICAL wdDryForcing, wdIceTopMelt, wdIceSeal, wdIceOBpack
      LOGICAL wdTopMeltSheet, wdTopMeltFloat
      COMMON /WETDRY_PARAMS_R/
     &     wdMinDepth, wdCritDepth, wdDragDepth, wdAdvDepth,
     &     wdGMDepth, wdIceDepth, wdLoadDepth, wdLoadTau,
     &     wdIceMinHeff, wdIceMaxSlip, wdIceBondH,
     &     wdBondHead, wdBondTau, wdBondLag, wdPondMax, wdBondFormEnd,
     &     wdKPPDepth, wdKPPCapDepth,
     &     wdSmoothWidth, wdMonFreq, wdMaxSpeed,
     &     wdManningN, wdDragMax, wdMaxFroude,
     &     wdManningT1, wdManningT2
      COMMON /WETDRY_PARAMS_C/ wdManningFile
      COMMON /WETDRY_PARAMS_L/
     &     wdConserveVol, wdCarryVel, wdIceGround, wdUpwindFace,
     &     wdDryForcing, wdIceTopMelt, wdIceSeal, wdIceOBpack,
     &     wdTopMeltSheet, wdTopMeltFloat

C--   wetdry fields
C     wdMaskW/S  :: 1 = face open, 0 = face closed by wetdry
C     wdMaskW0/S0:: static (initial) maskW/S, before wetdry masking
C     wdDryC     :: 1 where the column is at the film floor (dry)
C     wdVolAdj   :: cumulative volume added by the film floor, net of
C                    what was taken from wet neighbours [m^3]
C     wdVolAdjStep :: same, for the latest time step [m^3]
C     wdNflips   :: number of face open/close changes, latest step
C     wdRelaxW/S :: thin-water velocity relaxation factor at U/V points
C                    (1 = no relaxation, 0 = velocity set to zero)
C     wdAdvFacW/S:: momentum-advection factor at U/V points (0..1)
C     wdIceFacW/S:: sea-ice tilt-force factor at U/V points (0..1)
C     wdLoadFac  :: fraction of the sea-ice load a column passes to the
C                    ocean surface pressure and the ice tilt (0..1, see
C                    wdLoadDepth, wdLoadTau; 1 when off)
C     wdLoadRead :: wdLoadFac was read from a pickup
C     wdOutFac   :: outflow limiter factor of the latest step: outgoing
C                    velocities of a column are scaled by this so that
C                    it keeps at least wdMinDepth (1 = not limited)
C     wdNlimit   :: number of columns limited, latest step
C     wdVol0     :: initial volume of the interior (maskInC) domain [m^3]
C     wdCumIn    :: cumulative volume that entered the interior domain
C                    through its lateral boundary (e.g. OBCS) [m^3]
      _RS wdMaskW  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS wdMaskS  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS wdDryC   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS wdMaskW0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS wdMaskS0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
C     wdSIMaskU0/V0 :: static seaice masks SIMaskU/V (pkg/seaice), which
C                       are multiplied by wdMaskW/S when useSEAICE
C     wdSIsolU0/V0  :: same for seaiceMaskU/V, the masks of the ice
C                       momentum solvers (C-grid)
      _RS wdSIMaskU0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS wdSIMaskV0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS wdSIsolU0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS wdSIsolV0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     wdManEx    :: extra Manning roughness at tracer points
C                    (wdManningFile; 0 without it)
      _RS wdManEx  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /WETDRY_FIELDS_RS/
     &     wdMaskW, wdMaskS, wdDryC, wdMaskW0, wdMaskS0,
     &     wdSIMaskU0, wdSIMaskV0, wdSIsolU0, wdSIsolV0,
     &     wdManEx
      _RL wdRelaxW (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL wdRelaxS (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL wdAdvFacW(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL wdAdvFacS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL wdIceFacW(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL wdIceFacS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL wdOutFac (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL wdLoadFac(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL wdVolAdj, wdVolAdjStep, wdVol0, wdCumIn
C     wdEvLim    :: net evaporation (EmPmR > 0) removed by the dry-column
C                    ramp and the outflow limiter, latest step [kg/m^2/s]
      _RL wdEvLim  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /WETDRY_FIELDS_RL/
     &     wdRelaxW, wdRelaxS, wdAdvFacW, wdAdvFacS,
     &     wdIceFacW, wdIceFacS, wdOutFac, wdLoadFac,
     &     wdVolAdj, wdVolAdjStep, wdVol0, wdCumIn, wdEvLim
C     wdUpW/S    :: sign (+1, -1, 0) of the surface velocity that chose
C                    the upwind side of each face (wdUpwindFace); set in
C                    the in-step WETDRY_UPDATE_MASKS and saved in the pickup
C                    (pickup_wetdryface) so that a restart uses the same side
      _RL wdUpW(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL wdUpS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /WETDRY_UPWIND_RL/ wdUpW, wdUpS
C     wdQswIn/Out :: Qsw before / after the dry-column scaling of the
C                     last call of WETDRY_DRY_FORCING; Qsw is not reloaded
C                     every step, so an unchanged Qsw (= wdQswOut) is
C                     scaled again from wdQswIn, not from itself
      _RS wdQswIn  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS wdQswOut (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /WETDRY_QSW_RS/ wdQswIn, wdQswOut
C     wdSIheff/area/snow/EmPmR :: pkg/seaice state before SEAICE_GROWTH,
C                       for the dry-column ramp of the thermodynamic
C                       change (S/R WETDRY_SEAICE_DRY)
      _RL wdSIheff (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL wdSIarea (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL wdSIsnow (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL wdSIempr (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL wdSIsalt (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /WETDRY_SEAICE_RL/
     &     wdSIheff, wdSIarea, wdSIsnow, wdSIempr, wdSIsalt
C     wdBond     :: 1 where the ice is frozen to the bed (wdIceBondH);
C                    between 0 and 1 while a released bond ramps off
C                    (wdBondTau)
C     wdBondRead :: wdBond was read from a pickup
      _RL wdBond(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /WETDRY_BOND_RL/ wdBond
C     wdMeltSh   :: top melt of bonded ice kept on the ice (wdTopMeltSheet)
C                    [kg m^-2 s^-1], set in WETDRY_SEAICE_DRY and added to
C                    the overflood sheet in the next OVERFLOOD_STEP (not in
C                    the pickup: one step of melt is lost at a restart)
      _RL wdMeltSh(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /WETDRY_MELTSH_RL/ wdMeltSh
C     wdMeltAtm  :: ice + snow melted by the atmosphere in the latest
C                    SEAICE_GROWTH [kg m^-2] (wdTopMeltSheet keeps only it)
      _RL wdMeltAtm(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /WETDRY_MELTATM_RL/ wdMeltAtm
C     wdLiftT    :: time the head release condition has held for a bonded
C                    column (wdBondLag) [s]; 0 once the bond is gone
      _RL wdLiftT(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /WETDRY_LIFTT_RL/ wdLiftT
      LOGICAL wdLoadRead, wdBondRead
      COMMON /WETDRY_LOAD_L/ wdLoadRead, wdBondRead
      INTEGER wdNflips, wdNlimit
      COMMON /WETDRY_FIELDS_I/
     &     wdNflips, wdNlimit
