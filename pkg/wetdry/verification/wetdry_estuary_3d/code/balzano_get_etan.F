#include "CPP_OPTIONS.h"

CBOP
C     !ROUTINE: BALZANO_GET_ETAN
C     !INTERFACE:
      SUBROUTINE BALZANO_GET_ETAN( i, j, bi, bj, etaOut, myThid )
C     !DESCRIPTION:
C     Return etaN(i,j,bi,bj) to OBCS_CALC, which cannot include
C     DYNVARS.h (its dummy arguments share names with DYNVARS fields).
C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DYNVARS.h"
C     !INPUT/OUTPUT PARAMETERS:
      INTEGER i, j, bi, bj
      _RL etaOut
      INTEGER myThid
CEOP
      etaOut = etaN(i,j,bi,bj)
      RETURN
      END
