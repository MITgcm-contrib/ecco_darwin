#include "DARWIN_OPTIONS.h"

C malibu_v01 [A12]: kelp column state (vertical transfers and frond length)
C is set to zero at the start and kept across restarts in pickup_kelp.<suff>
C (64-bit: growth_transfer, mort_transfer, Qmag_transfer, mag_length,
C mag_blade_frac [Nr levels each], then mag_cum_length [1 level]).
C A missing file on restart leaves the state at zero.

CBOP
C !ROUTINE: DARWIN_KELP_INIT
C !INTERFACE: ==========================================================
      SUBROUTINE DARWIN_KELP_INIT( myThid )
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "DARWIN_SIZE.h"
#include "DARWIN_FIELDS.h"
      INTEGER myThid
CEOP
#ifdef DARWIN_ALLOW_MACROALGAE
      INTEGER i, j, k, bi, bj
      CHARACTER*(MAX_LEN_FNAM) fn, filNam
      LOGICAL fileExist

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          mag_tot_length(i,j,bi,bj) = 0. _d 0
          mag_cum_length(i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
        DO k=1,Nr
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           growth_transfer(i,j,k,bi,bj) = 0. _d 0
           mort_transfer(i,j,k,bi,bj) = 0. _d 0
           Qmag_transfer(i,j,k,bi,bj) = 0. _d 0
           mag_length(i,j,k,bi,bj) = 0. _d 0
           mag_blade_frac(i,j,k,bi,bj) = 1. _d 0
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      IF ( nIter0.GT.PTRACERS_Iter0 ) THEN
       IF ( pickupSuff.EQ.' ' ) THEN
        WRITE(fn,'(A,I10.10)') 'pickup_kelp.', nIter0
       ELSE
        WRITE(fn,'(A,A10)') 'pickup_kelp.', pickupSuff
       ENDIF
       CALL MDS_CHECK4FILE( fn, '.data', 'DARWIN_KELP_INIT',
     &                      filNam, fileExist, .TRUE., myThid )
       IF ( fileExist ) THEN
        CALL READ_REC_3D_RL( fn, precFloat64, Nr, growth_transfer,
     &                       1, nIter0, myThid )
        CALL READ_REC_3D_RL( fn, precFloat64, Nr, mort_transfer,
     &                       2, nIter0, myThid )
        CALL READ_REC_3D_RL( fn, precFloat64, Nr, Qmag_transfer,
     &                       3, nIter0, myThid )
        CALL READ_REC_3D_RL( fn, precFloat64, Nr, mag_length,
     &                       4, nIter0, myThid )
        CALL READ_REC_3D_RL( fn, precFloat64, Nr, mag_blade_frac,
     &                       5, nIter0, myThid )
        CALL READ_REC_3D_RL( fn, precFloat64, 1, mag_cum_length,
     &                       5*Nr+1, nIter0, myThid )
       ENDIF
      ENDIF
#endif /* DARWIN_ALLOW_MACROALGAE */
      RETURN
      END

CBOP
C !ROUTINE: DARWIN_KELP_WRITE_PICKUP
C !INTERFACE: ==========================================================
      SUBROUTINE DARWIN_KELP_WRITE_PICKUP( suff, myIter, myThid )
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DARWIN_SIZE.h"
#include "DARWIN_FIELDS.h"
      CHARACTER*(*) suff
      INTEGER myIter, myThid
CEOP
#ifdef DARWIN_ALLOW_MACROALGAE
      CHARACTER*(MAX_LEN_FNAM) fn
      WRITE(fn,'(A,A)') 'pickup_kelp.', suff
      CALL WRITE_REC_3D_RL( fn, precFloat64, Nr, growth_transfer,
     &                      1, myIter, myThid )
      CALL WRITE_REC_3D_RL( fn, precFloat64, Nr, mort_transfer,
     &                      2, myIter, myThid )
      CALL WRITE_REC_3D_RL( fn, precFloat64, Nr, Qmag_transfer,
     &                      3, myIter, myThid )
      CALL WRITE_REC_3D_RL( fn, precFloat64, Nr, mag_length,
     &                      4, myIter, myThid )
      CALL WRITE_REC_3D_RL( fn, precFloat64, Nr, mag_blade_frac,
     &                      5, myIter, myThid )
      CALL WRITE_REC_3D_RL( fn, precFloat64, 1, mag_cum_length,
     &                      5*Nr+1, myIter, myThid )
#endif /* DARWIN_ALLOW_MACROALGAE */
      RETURN
      END
