#include "WAD_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: WAD_UPDATE_MASKS

C !INTERFACE:
      SUBROUTINE WAD_UPDATE_MASKS(
     I                       inStep, myTime, myIter, myThid )

C     !DESCRIPTION:
C     Open or close each W and S face according to the local water level
C     (Casulli 2009, Stelling & Duinmeijer 2003):
C       z_bf = MAX( R_low(i-1), R_low(i) )             bed crest at face
C       h_f  = MAX( Ro_surf+etaH (i-1), Ro_surf+etaH (i) ) - z_bf
C       face is open if h_f > wadCritDepth
C     then set maskW/S = static mask * face mask. The face mask is also
C     applied to hFacW/S in S/R UPDATE_R_STAR (r*) or UPDATE_SURF_DR
C     (surface-level NLFS), which must be called after this routine.
C     When called within the time step (inStep=T), also:
C      - if wadCarryVel, give a face that just opened the velocity u* of
C        the upstream face of its wet cell, if that flows towards it
C        (otherwise the advancing front restarts from rest at every new
C        cell, a drag that does not decrease with resolution),
C      - relax u* (in gU,gV) to zero in thin water (Oey 2005),
C      - set the momentum-advection factor wadAdvFacW/S used in
C        S/R MOM_VECINV at the next time step,
C      - reset the Adams-Bashforth momentum history on closed faces.

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SURFACE.h"
#include "DYNVARS.h"
#include "WAD.h"

C     !INPUT PARAMETERS:
C     inStep :: called from within the time step (not initialisation)
C     myTime :: current time in simulation
C     myIter :: current iteration number
C     myThid :: my Thread Id number
      LOGICAL inStep
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_WAD
C     !LOCAL VARIABLES:
C     zbf     :: bed crest at the face [m]
C     hf      :: flow depth over the crest [m]
C     newMask :: updated face mask
C     nFlips  :: number of face changes (tile, then global)
      INTEGER i, j, k, bi, bj
      _RL zbf, hf, newMask, rlxRange, advRange
      _RL sumT, sumH
      INTEGER iUp, jUp, ks
      _RL nFlipsTile(nSx,nSy)
      _RL nFlips
C     openW/S :: 1 where the face opened at this call (else 0)
C     srcU/V  :: u*,v* taken from the upstream face, for opened faces
      _RS openW(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS openS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL srcU (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL srcV (1-OLx:sNx+OLx,1-OLy:sNy+OLy)

      rlxRange = wadDragDepth - wadMinDepth
      advRange = wadAdvDepth  - wadCritDepth
      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        nFlipsTile(bi,bj) = 0. _d 0
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          openW(i,j,bi,bj) = zeroRS
          openS(i,j,bi,bj) = zeroRS
         ENDDO
        ENDDO

C--   West faces
        DO j=1,sNy
         DO i=1,sNx+1
          newMask = 0. _d 0
          wadRelaxW(i,j,bi,bj) = 0. _d 0
          wadAdvFacW(i,j,bi,bj) = 0. _d 0
          IF ( kSurfW(i,j,bi,bj).LE.Nr ) THEN
           zbf = MAX( R_low(i-1,j,bi,bj), R_low(i,j,bi,bj) )
           hf  = MAX( Ro_surf(i-1,j,bi,bj) + etaH(i-1,j,bi,bj),
     &                Ro_surf(i,j,bi,bj)   + etaH(i,j,bi,bj) ) - zbf
           IF ( hf.GT.wadCritDepth ) newMask = 1. _d 0
           wadRelaxW(i,j,bi,bj) = newMask
           IF ( rlxRange.GT.zeroRL ) wadRelaxW(i,j,bi,bj) = newMask
     &        *MIN( oneRL, MAX( zeroRL, (hf-wadMinDepth)/rlxRange ) )
           wadAdvFacW(i,j,bi,bj) = newMask
           IF ( advRange.GT.zeroRL ) wadAdvFacW(i,j,bi,bj) = newMask
     &        *MIN( oneRL, MAX( zeroRL, (hf-wadCritDepth)/advRange ) )
          ENDIF
C-    count each face once: face sNx+1 is face 1 of the next tile
          IF ( newMask.NE.wadMaskW(i,j,bi,bj) .AND. i.LE.sNx )
     &      nFlipsTile(bi,bj) = nFlipsTile(bi,bj) + 1. _d 0
          IF ( newMask.GT.wadMaskW(i,j,bi,bj) ) openW(i,j,bi,bj) = oneRS
          wadMaskW(i,j,bi,bj) = newMask
         ENDDO
        ENDDO

C--   South faces
        DO j=1,sNy+1
         DO i=1,sNx
          newMask = 0. _d 0
          wadRelaxS(i,j,bi,bj) = 0. _d 0
          wadAdvFacS(i,j,bi,bj) = 0. _d 0
          IF ( kSurfS(i,j,bi,bj).LE.Nr ) THEN
           zbf = MAX( R_low(i,j-1,bi,bj), R_low(i,j,bi,bj) )
           hf  = MAX( Ro_surf(i,j-1,bi,bj) + etaH(i,j-1,bi,bj),
     &                Ro_surf(i,j,bi,bj)   + etaH(i,j,bi,bj) ) - zbf
           IF ( hf.GT.wadCritDepth ) newMask = 1. _d 0
           wadRelaxS(i,j,bi,bj) = newMask
           IF ( rlxRange.GT.zeroRL ) wadRelaxS(i,j,bi,bj) = newMask
     &        *MIN( oneRL, MAX( zeroRL, (hf-wadMinDepth)/rlxRange ) )
           wadAdvFacS(i,j,bi,bj) = newMask
           IF ( advRange.GT.zeroRL ) wadAdvFacS(i,j,bi,bj) = newMask
     &        *MIN( oneRL, MAX( zeroRL, (hf-wadCritDepth)/advRange ) )
          ENDIF
          IF ( newMask.NE.wadMaskS(i,j,bi,bj) .AND. j.LE.sNy )
     &      nFlipsTile(bi,bj) = nFlipsTile(bi,bj) + 1. _d 0
          IF ( newMask.GT.wadMaskS(i,j,bi,bj) ) openS(i,j,bi,bj) = oneRS
          wadMaskS(i,j,bi,bj) = newMask
         ENDDO
        ENDDO

       ENDDO
      ENDDO

C-    Faces on tile edges (i=1, sNx+1 ...) are computed by both tiles
C     and give the same answer; exchange to fill the full overlap
      CALL EXCH_UV_XY_RS( wadMaskW, wadMaskS, .FALSE., myThid )
      CALL EXCH_UV_XY_RL( wadRelaxW, wadRelaxS, .FALSE., myThid )
      CALL EXCH_UV_XY_RL( wadAdvFacW, wadAdvFacS, .FALSE., myThid )

      CALL GLOBAL_SUM_TILE_RL( nFlipsTile, nFlips, myThid )
      _BEGIN_MASTER(myThid)
      wadNflips = NINT( nFlips )
      _END_MASTER(myThid)

C-    Opened faces on tile edges take u* from the overlap: it must be
C     valid there (DYNAMICS does not fill the full overlap)
      IF ( inStep .AND. wadCarryVel )
     &     CALL EXCH_UV_XYZ_RL( gU, gV, .TRUE., myThid )

C--   Apply the face mask to maskW/S (hFacW/S done in UPDATE_R_STAR)
      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)

#ifdef NONLIN_FRSURF
C-    The face thickness comes from the mean (r*: rStarFacW/S) or the
C     minimum (surface-level NLFS: hFac_surfW/S) of the two surface
C     heights and can drop to (or below) zero on an open face next to a
C     bed step: keep at least a film of wadMinDepth on open faces
C-    upwind side of each face: from the surface velocity within the time
C     step; at initialisation from wadUpW/S (WAD_INIT_VARIA: initial
C     velocity, or the pickup on a restart)
        IF ( wadUpwindFace .AND. inStep ) THEN
         DO j=2-OLy,sNy+OLy
          DO i=2-OLx,sNx+OLx
           ks = 1
           IF ( select_rStar.EQ.0 ) ks = MIN( kSurfW(i,j,bi,bj), Nr )
           wadUpW(i,j,bi,bj) = 0. _d 0
           IF ( uVel(i,j,ks,bi,bj).GT.zeroRL ) wadUpW(i,j,bi,bj) = 1.
           IF ( uVel(i,j,ks,bi,bj).LT.zeroRL ) wadUpW(i,j,bi,bj) = -1.
           ks = 1
           IF ( select_rStar.EQ.0 ) ks = MIN( kSurfS(i,j,bi,bj), Nr )
           wadUpS(i,j,bi,bj) = 0. _d 0
           IF ( vVel(i,j,ks,bi,bj).GT.zeroRL ) wadUpS(i,j,bi,bj) = 1.
           IF ( vVel(i,j,ks,bi,bj).LT.zeroRL ) wadUpS(i,j,bi,bj) = -1.
          ENDDO
         ENDDO
        ENDIF
        IF ( select_rStar.GT.0 ) THEN
         IF ( wadUpwindFace ) THEN
C-    face thickness from the upwind surface height (default: mean)
          DO j=2-OLy,sNy+OLy
           DO i=2-OLx,sNx+OLx
            IF ( wadMaskW(i,j,bi,bj).NE.zeroRS ) THEN
             IF ( wadUpW(i,j,bi,bj).GT.zeroRL ) THEN
              rStarFacW(i,j,bi,bj) = 1. _d 0 + etaH(i-1,j,bi,bj)
     &          /( rSurfW(i,j,bi,bj) - rLowW(i,j,bi,bj) )
             ELSEIF ( wadUpW(i,j,bi,bj).LT.zeroRL ) THEN
              rStarFacW(i,j,bi,bj) = 1. _d 0 + etaH(i,j,bi,bj)
     &          /( rSurfW(i,j,bi,bj) - rLowW(i,j,bi,bj) )
             ENDIF
            ENDIF
            IF ( wadMaskS(i,j,bi,bj).NE.zeroRS ) THEN
             IF ( wadUpS(i,j,bi,bj).GT.zeroRL ) THEN
              rStarFacS(i,j,bi,bj) = 1. _d 0 + etaH(i,j-1,bi,bj)
     &          /( rSurfS(i,j,bi,bj) - rLowS(i,j,bi,bj) )
             ELSEIF ( wadUpS(i,j,bi,bj).LT.zeroRL ) THEN
              rStarFacS(i,j,bi,bj) = 1. _d 0 + etaH(i,j,bi,bj)
     &          /( rSurfS(i,j,bi,bj) - rLowS(i,j,bi,bj) )
             ENDIF
            ENDIF
           ENDDO
          ENDDO
         ENDIF
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           IF ( wadMaskW(i,j,bi,bj).NE.zeroRS ) rStarFacW(i,j,bi,bj) =
     &       MAX( rStarFacW(i,j,bi,bj), wadMinDepth
     &            /( rSurfW(i,j,bi,bj) - rLowW(i,j,bi,bj) ) )
           IF ( wadMaskS(i,j,bi,bj).NE.zeroRS ) rStarFacS(i,j,bi,bj) =
     &       MAX( rStarFacS(i,j,bi,bj), wadMinDepth
     &            /( rSurfS(i,j,bi,bj) - rLowS(i,j,bi,bj) ) )
          ENDDO
         ENDDO
        ELSE
         IF ( wadUpwindFace ) THEN
C-    face thickness from the upwind surface height (CALC_SURF_DR uses
C     the lower one, which throttles the flux at the shoreline)
          DO j=2-OLy,sNy+OLy
           DO i=2-OLx,sNx+OLx
            ks = kSurfW(i,j,bi,bj)
            IF ( wadMaskW(i,j,bi,bj).NE.zeroRS .AND. ks.LE.Nr ) THEN
             IF ( wadUpW(i,j,bi,bj).GT.zeroRL ) THEN
              hFac_surfW(i,j,bi,bj) = h0FacW(i,j,ks,bi,bj)
     &          + ( Ro_surf(i-1,j,bi,bj) + etaH(i-1,j,bi,bj)
     &            - rSurfW(i,j,bi,bj) )*recip_drF(ks)
             ELSEIF ( wadUpW(i,j,bi,bj).LT.zeroRL ) THEN
              hFac_surfW(i,j,bi,bj) = h0FacW(i,j,ks,bi,bj)
     &          + ( Ro_surf(i,j,bi,bj) + etaH(i,j,bi,bj)
     &            - rSurfW(i,j,bi,bj) )*recip_drF(ks)
             ENDIF
            ENDIF
            ks = kSurfS(i,j,bi,bj)
            IF ( wadMaskS(i,j,bi,bj).NE.zeroRS .AND. ks.LE.Nr ) THEN
             IF ( wadUpS(i,j,bi,bj).GT.zeroRL ) THEN
              hFac_surfS(i,j,bi,bj) = h0FacS(i,j,ks,bi,bj)
     &          + ( Ro_surf(i,j-1,bi,bj) + etaH(i,j-1,bi,bj)
     &            - rSurfS(i,j,bi,bj) )*recip_drF(ks)
             ELSEIF ( wadUpS(i,j,bi,bj).LT.zeroRL ) THEN
              hFac_surfS(i,j,bi,bj) = h0FacS(i,j,ks,bi,bj)
     &          + ( Ro_surf(i,j,bi,bj) + etaH(i,j,bi,bj)
     &            - rSurfS(i,j,bi,bj) )*recip_drF(ks)
             ENDIF
            ENDIF
           ENDDO
          ENDDO
         ENDIF
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           ks = kSurfW(i,j,bi,bj)
           IF ( wadMaskW(i,j,bi,bj).NE.zeroRS .AND. ks.LE.Nr )
     &       hFac_surfW(i,j,bi,bj) = MAX( hFac_surfW(i,j,bi,bj),
     &                                    wadMinDepth*recip_drF(ks) )
           ks = kSurfS(i,j,bi,bj)
           IF ( wadMaskS(i,j,bi,bj).NE.zeroRS .AND. ks.LE.Nr )
     &       hFac_surfS(i,j,bi,bj) = MAX( hFac_surfS(i,j,bi,bj),
     &                                    wadMinDepth*recip_drF(ks) )
          ENDDO
         ENDDO
        ENDIF
#endif /* NONLIN_FRSURF */

        DO k=1,Nr
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           maskW(i,j,k,bi,bj) = wadMaskW0(i,j,k,bi,bj)
     &                        * wadMaskW(i,j,bi,bj)
           maskS(i,j,k,bi,bj) = wadMaskS0(i,j,k,bi,bj)
     &                        * wadMaskS(i,j,bi,bj)
          ENDDO
         ENDDO
        ENDDO
        IF ( inStep .AND. wadCarryVel ) THEN
C-    Opened faces: carry the depth-mean u* of the upstream face of the
C     wet cell onto all levels of the new face (a level-by-level copy
C     is unstable with Nr > 1). Source values are computed first, so
C     that successive opened faces do not chain.
         DO j=1,sNy
          DO i=1,sNx+1
           srcU(i,j) = 0. _d 0
           IF ( openW(i,j,bi,bj).NE.zeroRS ) THEN
            iUp = i+1
            IF ( Ro_surf(i-1,j,bi,bj)+etaH(i-1,j,bi,bj) .GE.
     &           Ro_surf(i,j,bi,bj)  +etaH(i,j,bi,bj) ) iUp = i-1
            sumT = 0. _d 0
            sumH = 0. _d 0
            DO k=1,Nr
             sumT = sumT + gU(iUp,j,k,bi,bj)*hFacW(iUp,j,k,bi,bj)*drF(k)
             sumH = sumH + hFacW(iUp,j,k,bi,bj)*drF(k)
            ENDDO
            IF ( sumH.GT.zeroRL ) THEN
             IF ( iUp.LT.i ) srcU(i,j) = MAX( sumT/sumH, zeroRL )
             IF ( iUp.GT.i ) srcU(i,j) = MIN( sumT/sumH, zeroRL )
            ENDIF
           ENDIF
          ENDDO
         ENDDO
         DO j=1,sNy+1
          DO i=1,sNx
           srcV(i,j) = 0. _d 0
           IF ( openS(i,j,bi,bj).NE.zeroRS ) THEN
            jUp = j+1
            IF ( Ro_surf(i,j-1,bi,bj)+etaH(i,j-1,bi,bj) .GE.
     &           Ro_surf(i,j,bi,bj)  +etaH(i,j,bi,bj) ) jUp = j-1
            sumT = 0. _d 0
            sumH = 0. _d 0
            DO k=1,Nr
             sumT = sumT + gV(i,jUp,k,bi,bj)*hFacS(i,jUp,k,bi,bj)*drF(k)
             sumH = sumH + hFacS(i,jUp,k,bi,bj)*drF(k)
            ENDDO
            IF ( sumH.GT.zeroRL ) THEN
             IF ( jUp.LT.j ) srcV(i,j) = MAX( sumT/sumH, zeroRL )
             IF ( jUp.GT.j ) srcV(i,j) = MIN( sumT/sumH, zeroRL )
            ENDIF
           ENDIF
          ENDDO
         ENDDO
         DO k=1,Nr
          DO j=1,sNy
           DO i=1,sNx+1
            IF ( openW(i,j,bi,bj).NE.zeroRS )
     &        gU(i,j,k,bi,bj) = srcU(i,j)*maskW(i,j,k,bi,bj)
           ENDDO
          ENDDO
          DO j=1,sNy+1
           DO i=1,sNx
            IF ( openS(i,j,bi,bj).NE.zeroRS )
     &        gV(i,j,k,bi,bj) = srcV(i,j)*maskS(i,j,k,bi,bj)
           ENDDO
          ENDDO
         ENDDO
        ENDIF

        IF ( inStep ) THEN
         DO k=1,Nr
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
C-    Thin-water relaxation of u*,v* (zero on closed faces)
            gU(i,j,k,bi,bj) = gU(i,j,k,bi,bj)*wadRelaxW(i,j,bi,bj)
            gV(i,j,k,bi,bj) = gV(i,j,k,bi,bj)*wadRelaxS(i,j,bi,bj)
C-    Forget the momentum tendency history on closed faces
#ifdef ALLOW_ADAMSBASHFORTH_3
            guNm(i,j,k,bi,bj,1) = guNm(i,j,k,bi,bj,1)*maskW(i,j,k,bi,bj)
            guNm(i,j,k,bi,bj,2) = guNm(i,j,k,bi,bj,2)*maskW(i,j,k,bi,bj)
            gvNm(i,j,k,bi,bj,1) = gvNm(i,j,k,bi,bj,1)*maskS(i,j,k,bi,bj)
            gvNm(i,j,k,bi,bj,2) = gvNm(i,j,k,bi,bj,2)*maskS(i,j,k,bi,bj)
#else
            guNm1(i,j,k,bi,bj) = guNm1(i,j,k,bi,bj)*maskW(i,j,k,bi,bj)
            gvNm1(i,j,k,bi,bj) = gvNm1(i,j,k,bi,bj)*maskS(i,j,k,bi,bj)
#endif
           ENDDO
          ENDDO
         ENDDO
        ENDIF

       ENDDO
      ENDDO

C-    Refresh the overlap after the changes above (faces on tile edges
C     are set by both tiles and must agree for conservation)
      IF ( inStep .AND. wadCarryVel )
     &     CALL EXCH_UV_XYZ_RL( gU, gV, .TRUE., myThid )

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics .AND. inStep ) THEN
        CALL DIAGNOSTICS_FILL_RS( wadMaskW, 'WADmaskW', 0, 1, 0,1,1,
     &                            myThid )
        CALL DIAGNOSTICS_FILL_RS( wadMaskS, 'WADmaskS', 0, 1, 0,1,1,
     &                            myThid )
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

C--   depth-dependent (Manning) bottom drag from the new face depths
      IF ( wadManningN.GT.zeroRL .OR. wadManningFile.NE.' ' )
     &     CALL WAD_BOTDRAG( myTime, myThid )

#endif /* ALLOW_WAD */

      RETURN
      END
