#include "SEDIMENT_OPTIONS.h"
#ifdef ALLOW_PTRACERS
# include "PTRACERS_OPTIONS.h"
#endif

CBOP
C !ROUTINE: SEDIMENT_ADD_TENDENCY
C !INTERFACE:
      SUBROUTINE SEDIMENT_ADD_TENDENCY(
     U                    gPtracer,
     I                    iMin,iMax,jMin,jMax, k, bi, bj,
     I                    iTracer, myTime, myIter, myThid )
C !DESCRIPTION:
C     Sediment carried by the added fluid (addMass, e.g. a river): the
C     tracer advection lets added fluid enter at the local concentration,
C     so, as for T and S (temp_addMass, salt_addMass in APPLY_FORCING),
C     add  addMass*(sedRiverC - C)/(rhoConst*rA*drF*hFacC)  to the
C     tendency of each sediment pTracer. Called from
C     PTRACERS_APPLY_FORCING.
C !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "FFIELDS.h"
#include "SEDIMENT.h"
#ifdef ALLOW_PTRACERS
# include "PTRACERS_SIZE.h"
# include "PTRACERS_FIELDS.h"
#endif
C !INPUT/OUTPUT PARAMETERS:
      _RL gPtracer(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin,iMax,jMin,jMax, k, bi, bj, iTracer
      _RL myTime
      INTEGER myIter, myThid
CEOP

#if (defined ALLOW_SEDIMENT) && (defined ALLOW_PTRACERS) && (defined ALLOW_ADDFLUID)
C !LOCAL VARIABLES:
      INTEGER i, j, n

      IF ( .NOT.useSEDIMENT .OR. selectAddFluid.LT.1 ) RETURN
      DO n=1,sedNcl
       IF ( sedPtrIdx(n).EQ.iTracer ) THEN
        DO j=jMin,jMax
         DO i=iMin,iMax
C-    only where fluid is added: withdrawn fluid (addMass < 0) leaves at
C     the local concentration, which the advection term already does
          gPtracer(i,j) = gPtracer(i,j)
     &      + MAX( addMass(i,j,k,bi,bj), zeroRL )*mass2rUnit
     &        *( sedRivC(i,j,bi,bj,n) - pTracer(i,j,k,bi,bj,iTracer) )
     &        *recip_rA(i,j,bi,bj)
     &        *recip_drF(k)*_recip_hFacC(i,j,k,bi,bj)
         ENDDO
        ENDDO
       ENDIF
      ENDDO
#endif

      RETURN
      END
