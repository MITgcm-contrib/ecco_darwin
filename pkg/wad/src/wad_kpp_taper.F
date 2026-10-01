#include "WAD_OPTIONS.h"
#ifdef ALLOW_KPP
# include "KPP_OPTIONS.h"
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: WAD_KPP_TAPER

C !INTERFACE:
      SUBROUTINE WAD_KPP_TAPER( bi, bj, myThid )

C !DESCRIPTION:
C     Limit the KPP non-local transport in thin water (r* columns):
C     1) Cap the non-local flux fraction. The non-local flux through the
C        top of level k is  Kz(k)*KPPghat(k-1)*Q  (Q the surface flux) and
C        must not exceed Q. With weak surface forcing KPPghat ~ 1/(ws*hbl)
C        is large, and in fast thin-water flow the shear-driven interior
C        part of Kz is large too, so the product can exceed 1: the term
C        then moves more heat or salt than the surface supplies (winter
C        mudflat, -1.5 C water under ice: T from -2.6 to +4.6 C).
C        KPPghat(k-1) is limited to 1/max(KPPdiffKzT(k), KPPdiffKzS(k)),
C        only in columns shallower than wadKPPCapDepth (deeper columns
C        keep stock KPP).
C     2) Optionally (wadKPPDepth > wadCritDepth) ramp KPPghat to zero in
C        columns shallower than wadKPPDepth, with the factor
C        fC = MIN(1, MAX(0, (D - wadCritDepth)/(wadKPPDepth - wadCritDepth)))
C     Called at the end of S/R KPP_CALC when useWAD.

C !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "WAD.h"
#ifdef ALLOW_KPP
# include "KPP.h"
#endif

C !INPUT PARAMETERS:
      INTEGER bi, bj
      INTEGER myThid
CEOP

#if (defined ALLOW_WAD) && (defined ALLOW_KPP) && (defined KPP_GHAT)
C !LOCAL VARIABLES:
      INTEGER i, j, k
      _RL fC, kpRange, kzMax

      IF ( .NOT.useWAD ) RETURN

C--   cap the non-local flux fraction at 1 in shallow columns
      IF ( wadKPPCapDepth.GT.zeroRL ) THEN
       DO k=2,Nr
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          kzMax = MAX( KPPdiffKzT(i,j,k,bi,bj),
     &                 KPPdiffKzS(i,j,k,bi,bj) )
          IF ( kzMax*KPPghat(i,j,k-1,bi,bj).GT.oneRL .AND.
     &         etaH(i,j,bi,bj) + Ro_surf(i,j,bi,bj) - R_low(i,j,bi,bj)
     &         .LT.wadKPPCapDepth )
     &      KPPghat(i,j,k-1,bi,bj) = oneRL/kzMax
         ENDDO
        ENDDO
       ENDDO
      ENDIF

      kpRange = wadKPPDepth - wadCritDepth
      IF ( kpRange.LE.zeroRL ) RETURN

      DO j=1-OLy,sNy+OLy
       DO i=1-OLx,sNx+OLx
        fC = MIN( oneRL, MAX( zeroRL,
     &     ( etaH(i,j,bi,bj) + Ro_surf(i,j,bi,bj) - R_low(i,j,bi,bj)
     &       - wadCritDepth )/kpRange ) )
        IF ( fC.LT.oneRL ) THEN
         DO k=1,Nr
          KPPghat(i,j,k,bi,bj) = KPPghat(i,j,k,bi,bj)*fC
         ENDDO
        ENDIF
       ENDDO
      ENDDO
#endif /* ALLOW_WAD and ALLOW_KPP and KPP_GHAT */

      RETURN
      END
