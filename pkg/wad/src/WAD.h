CBOP
C     !ROUTINE: WAD.h
C     !INTERFACE:
C     #include "WAD.h"

C     !DESCRIPTION:
C     *================================================================*
C     | WAD.h
C     | o Header file for pkg "wad" (wetting and drying): thin-film
C     |   scheme with dynamic face masking, for r* non-linear
C     |   free surface.
C     *================================================================*
CEOP

C--   WAD parameters (namelist WAD_PARM01 in data.wad)
C     wadMinDepth    :: thickness of the water film kept in a dry column [m]
C     wadCritDepth   :: a face is open when the highest surface on either
C                       side stands more than this above the face bed
C                       crest [m] ; must be > wadMinDepth
C     wadDragDepth   :: below this face depth, velocity is relaxed
C                       linearly to zero at wadMinDepth [m] ;
C                       no relaxation if <= wadMinDepth (default 0)
C     wadAdvDepth    :: below this face flow depth, momentum advection
C                       (pkg/mom_vecinv) is ramped linearly to zero at
C                       wadCritDepth [m] ; off if <= wadCritDepth
C     wadGMDepth     :: below this column depth, the GM/Redi tensor
C                       (pkg/gmredi) is ramped linearly to zero at
C                       wadCritDepth [m] ; off if <= wadCritDepth
C     wadKPPDepth    :: below this column depth, the KPP non-local
C                       transport (KPPghat) is ramped linearly to zero at
C                       wadCritDepth [m] ; off if <= wadCritDepth
C                       (default)
C     wadKPPCapDepth :: in columns shallower than this [m], the KPP
C                       non-local flux fraction Kz*KPPghat is capped at 1
C                       (WAD_KPP_TAPER); deeper columns keep stock KPP.
C                       <= 0: no cap
C     wadSmoothWidth :: width of tanh face-mask ramp [m] (not implemented)
C     wadMonFreq     :: frequency of WAD monitor output [s]
C     wadMaxSpeed    :: stop, and report where, when a velocity on an
C                       open face exceeds this [m/s] ; off if <= 0
C     wadConserveVol :: remove the volume added by the film floor from
C                       wet neighbours (through open faces)
C     wadUpwindFace  :: thickness of an open face from the
C                       upwind surface height (default ; if .FALSE.: the
C                       lower one with select_rStar=0, the mean with r*)
C     wadDryForcing  :: no surface heat/salt/fresh-water/short-wave or
C                       pTracer surface forcing on dry columns (ramp
C                       from wadCritDepth to 2*wadCritDepth)
C     wadManningN    :: Manning roughness [s m^-1/3] for a depth-dependent
C                       quadratic drag Cd = g n^2/D^(1/3), limited to
C                       bottomDragQuadratic .. wadDragMax (WAD_BOTDRAG,
C                       needs ALLOW_BOTTOMDRAG_ROUGHNESS) ; off if <= 0
C     wadDragMax     :: upper limit of that drag coefficient
C     wadManningFile :: optional 2-D extra Manning roughness [s m^-1/3] at
C                       tracer points (e.g. river ice cover and ice jams),
C                       added to wadManningN on each face (the larger of
C                       the two neighbours), times a factor 1 until
C                       wadManningT1, falling linearly to 0 at
C                       wadManningT2 [s, model time] (constant if
C                       wadManningT2 <= wadManningT1)
C     wadMaxFroude   :: cap the speed on each open interior face at
C                       wadMaxFroude*sqrt(g D) (D the face flow depth, as
C                       for the Manning drag), applied to u(n+1) before
C                       the outflow limiter, so volume and tracers stay
C                       conserved ; for fast thin sheets pouring off
C                       banks (breakup floods) ; off if <= 0
C     wadCarryVel    :: a face that opens takes the depth-mean u* of the
C                       upstream face of its wet cell (the advancing
C                       front then converges with resolution) ; with
C                       Nr > 1 it needs momImplVertAdv=.TRUE.
      _RL wadMinDepth
      _RL wadCritDepth
      _RL wadDragDepth
      _RL wadAdvDepth
      _RL wadGMDepth
      _RL wadKPPDepth, wadKPPCapDepth
      _RL wadSmoothWidth
      _RL wadMonFreq
      _RL wadMaxSpeed
      _RL wadManningN, wadDragMax, wadMaxFroude
      _RL wadManningT1, wadManningT2
      CHARACTER*(MAX_LEN_FNAM) wadManningFile
      LOGICAL wadConserveVol, wadCarryVel, wadUpwindFace
      LOGICAL wadDryForcing
      COMMON /WAD_PARAMS_R/
     &     wadMinDepth, wadCritDepth, wadDragDepth, wadAdvDepth,
     &     wadGMDepth,
     &     wadKPPDepth, wadKPPCapDepth,
     &     wadSmoothWidth, wadMonFreq, wadMaxSpeed,
     &     wadManningN, wadDragMax, wadMaxFroude,
     &     wadManningT1, wadManningT2
      COMMON /WAD_PARAMS_C/ wadManningFile
      COMMON /WAD_PARAMS_L/
     &     wadConserveVol, wadCarryVel, wadUpwindFace,
     &     wadDryForcing

C--   WAD fields
C     wadMaskW/S  :: 1 = face open, 0 = face closed by WAD
C     wadMaskW0/S0:: static (initial) maskW/S, before WAD masking
C     wadDryC     :: 1 where the column is at the film floor (dry)
C     wadVolAdj   :: cumulative volume added by the film floor, net of
C                    what was taken from wet neighbours [m^3]
C     wadVolAdjStep :: same, for the latest time step [m^3]
C     wadNflips   :: number of face open/close changes, latest step
C     wadRelaxW/S :: thin-water velocity relaxation factor at U/V points
C                    (1 = no relaxation, 0 = velocity set to zero)
C     wadAdvFacW/S:: momentum-advection factor at U/V points (0..1)
C     wadOutFac   :: outflow limiter factor of the latest step: outgoing
C                    velocities of a column are scaled by this so that
C                    it keeps at least wadMinDepth (1 = not limited)
C     wadNlimit   :: number of columns limited, latest step
C     wadVol0     :: initial volume of the interior (maskInC) domain [m^3]
C     wadCumIn    :: cumulative volume that entered the interior domain
C                    through its lateral boundary (e.g. OBCS) [m^3]
      _RS wadMaskW  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS wadMaskS  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS wadDryC   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS wadMaskW0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS wadMaskS0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
C     wadManEx    :: extra Manning roughness at tracer points
C                    (wadManningFile; 0 without it)
      _RS wadManEx  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /WAD_FIELDS_RS/
     &     wadMaskW, wadMaskS, wadDryC, wadMaskW0, wadMaskS0,
     &     wadManEx
      _RL wadRelaxW (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL wadRelaxS (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL wadAdvFacW(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL wadAdvFacS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL wadOutFac (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL wadVolAdj, wadVolAdjStep, wadVol0, wadCumIn
C     wadEvLim    :: net evaporation (EmPmR > 0) removed by the dry-column
C                    ramp and the outflow limiter, latest step [kg/m^2/s]
      _RL wadEvLim  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /WAD_FIELDS_RL/
     &     wadRelaxW, wadRelaxS, wadAdvFacW, wadAdvFacS,
     &     wadOutFac,
     &     wadVolAdj, wadVolAdjStep, wadVol0, wadCumIn, wadEvLim
C     wadUpW/S    :: sign (+1, -1, 0) of the surface velocity that chose
C                    the upwind side of each face (wadUpwindFace); set in
C                    the in-step WAD_UPDATE_MASKS and saved in the pickup
C                    (pickup_wadface) so that a restart uses the same side
      _RL wadUpW(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL wadUpS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /WAD_UPWIND_RL/ wadUpW, wadUpS
C     wadQswIn/Out :: Qsw before / after the dry-column scaling of the
C                     last call of WAD_DRY_FORCING; Qsw is not reloaded
C                     every step, so an unchanged Qsw (= wadQswOut) is
C                     scaled again from wadQswIn, not from itself
      _RS wadQswIn  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS wadQswOut (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /WAD_QSW_RS/ wadQswIn, wadQswOut
      INTEGER wadNflips, wadNlimit
      COMMON /WAD_FIELDS_I/
     &     wadNflips, wadNlimit
