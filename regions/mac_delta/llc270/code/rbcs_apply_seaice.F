C $Header: /u/gcmpack/MITgcm/pkg/rbcs/rbcs_add_tendency.F,v 1.12 2015/10/24 21:20:09 jmc Exp $
C $Name:  $

#include "RBCS_OPTIONS.h"

CBOP
C !ROUTINE: RBCS_APPLY_SEAICE

C !INTERFACE: ==========================================================
      SUBROUTINE RBCS_APPLY_SEAICE(
     I                    myTime, myIter, myThid )

C !DESCRIPTION:
C     Add to relaxation to seaice area and/or heff

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
c#include "GRID.h"
#include "DYNVARS.h"
#ifdef ALLOW_PTRACERS
#include "PTRACERS_SIZE.h"
#include "PTRACERS_FIELDS.h"
#endif
#include "RBCS_SIZE.h"
#include "RBCS_PARAMS.h"
#include "RBCS_FIELDS.h"
chzh[
#ifdef ALLOW_SEAICE
# include "SEAICE_SIZE.h"
# include "SEAICE_PARAMS.h"
# include "SEAICE.h"
#endif
chzh]

C !INPUT/OUTPUT PARAMETERS: ============================================
C  myTime         :: current time
C  myIter         :: current timestep
C  myThid         :: my Thread Id number
      _RL myTime
      INTEGER myIter
      INTEGER myThid

C !LOCAL VARIABLES: ====================================================
C  i,j            :: loop indices
      INTEGER i,j,bi,bj
      _RL rbcsVanishingFac
      _RL rec_tauRlx
#ifdef ALLOW_PTRACERS
      INTEGER iTracer
#endif
CEOP

#ifdef ALLOW_RBCS

      IF ( rbcsVanishingTime.GT.0. _d 0 ) THEN
        rbcsVanishingFac =
     &      MAX( 0. _d 0 , 1. _d 0 - myTime / rbcsVanishingTime )
      ELSE
        rbcsVanishingFac = 1. _d 0
      ENDIF

      IF ( useRBCarea ) THEN
        rec_tauRlx = rbcsVanishingFac/tauRelaxA*SEAICE_deltaTtherm
	DO bj=myByLo(myThid),myByHi(myThid)
	DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=0,sNy+1
         DO i=0,sNx+1
          AREA(i,j,bi,bj) = AREA(i,j,bi,bj)
     &       - RBC_maskA(i,j,bi,bj)*rec_tauRlx
     &        *( AREA(i,j,bi,bj)- RBCarea(i,j,bi,bj) )
         ENDDO
        ENDDO
        ENDDO
        ENDDO
      ENDIF
      IF ( useRBCheff ) THEN
        rec_tauRlx = rbcsVanishingFac/tauRelaxH*SEAICE_deltaTtherm
	DO bj=myByLo(myThid),myByHi(myThid)
	DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=0,sNy+1
         DO i=0,sNx+1
          HEFF(i,j,bi,bj) = HEFF(i,j,bi,bj)
     &       - RBC_maskH(i,j,bi,bj)*rec_tauRlx
     &        *( HEFF(i,j,bi,bj)- RBCheff(i,j,bi,bj) )
         ENDDO
        ENDDO
        ENDDO
        ENDDO
      ENDIF

#endif /* ALLOW_RBCS */

      RETURN
      END
