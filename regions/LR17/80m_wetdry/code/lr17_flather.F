#include "OBCS_OPTIONS.h"

CBOP
C     !ROUTINE: LR17_FLATHER
C     !INTERFACE:
      SUBROUTINE LR17_FLATHER( myThid )
C     !DESCRIPTION:
C     Flather (1976) radiation condition on the LR17 open boundaries
C     (west column, south row), called from OBCS_CALC right after
C     OBCS_PRESCRIBE_READ has set the IBI values (re-interpolated in time
C     at every call, so the correction does not accumulate):
C       u_b = u_IBI + sqrt(g/D_b) ( eta_IBI - eta_interior )
C     added uniformly over depth to the prescribed normal velocity
C     (west: OBWu, positive into the domain; south: OBSv, positive into
C     the domain). With the IBI velocity alone the volume entering the
C     domain is IBI's depth-mean current times the model's own boundary
C     depth, which differs from IBI's, and nothing pulls the interior
C     sea level back to the boundary level: the mean level drifted by
C     +-0.7 m with the spring-neap cycle. D_b is the boundary cell's
C     total depth at the IBI level, at least 0.5 m.
C     No correction within 2 cells of a corner where the two open
C     boundaries meet: there the interior neighbour of one boundary lies
C     on the other, and the two corrections drove a 6 m/s jet and a 10 m
C     surface depression in the south-west corner (2026-10-04).
C     Off with lr17Flather = .FALSE. in data.lr17ob (optional file).
C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "OBCS_PARAMS.h"
#include "OBCS_GRID.h"
#include "OBCS_FIELDS.h"
C     !INPUT/OUTPUT PARAMETERS:
      INTEGER myThid
CEOP
#if (defined ALLOW_OBCS) && (defined NONLIN_FRSURF)
C     !LOCAL VARIABLES:
      INTEGER bi, bj, i, j, k, iOb, jOb, iUnit
      _RL depOb, du
      LOGICAL exst, flInit, lr17Flather, corner
      COMMON /LR17_FLATHER_L/ lr17Flather
      SAVE flInit
      DATA flInit /.FALSE./
      NAMELIST /LR17OB_PARM01/ lr17Flather
      CHARACTER*(MAX_LEN_MBUF) msgBuf

      IF ( .NOT.flInit ) THEN
       _BARRIER
       _BEGIN_MASTER(myThid)
       lr17Flather = .TRUE.
       INQUIRE( FILE='data.lr17ob', EXIST=exst )
       IF ( exst ) THEN
        CALL OPEN_COPY_DATA_FILE( 'data.lr17ob', 'LR17_FLATHER',
     O                            iUnit, myThid )
        READ(UNIT=iUnit, NML=LR17OB_PARM01)
#ifdef SINGLE_DISK_IO
        CLOSE(iUnit)
#else
        CLOSE(iUnit,STATUS='DELETE')
#endif
       ENDIF
       WRITE(msgBuf,'(A,L2)') 'LR17_FLATHER: lr17Flather=', lr17Flather
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
       flInit = .TRUE.
       _END_MASTER(myThid)
       _BARRIER
      ENDIF
      IF ( .NOT.lr17Flather ) RETURN

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
#ifdef ALLOW_OBCS_WEST
        DO j=1-OLy,sNy+OLy
         iOb = OB_Iw(j,bi,bj)
         IF ( iOb.NE.OB_indexNone ) THEN
          corner = .FALSE.
#ifdef ALLOW_OBCS_SOUTH
          IF ( iOb+1.LE.sNx+OLx ) THEN
           IF ( OB_Js(iOb+1,bi,bj).NE.OB_indexNone )
     &       corner = j.LE.OB_Js(iOb+1,bi,bj)+2
          ENDIF
#endif
          IF ( iOb+1.LE.sNx+OLx .AND. .NOT.corner ) THEN
           depOb = MAX( OBWeta(j,bi,bj) + Ro_surf(iOb,j,bi,bj)
     &                - R_low(iOb,j,bi,bj), 0.5 _d 0 )
           du = SQRT( gravity/depOb )
     &        *( OBWeta(j,bi,bj) - etaN(iOb+1,j,bi,bj) )
           DO k=1,Nr
            OBWu(j,k,bi,bj) = OBWu(j,k,bi,bj)
     &                      + du*maskW(iOb+1,j,k,bi,bj)
           ENDDO
          ENDIF
         ENDIF
        ENDDO
#endif
#ifdef ALLOW_OBCS_SOUTH
        DO i=1-OLx,sNx+OLx
         jOb = OB_Js(i,bi,bj)
         IF ( jOb.NE.OB_indexNone ) THEN
          corner = .FALSE.
#ifdef ALLOW_OBCS_WEST
          IF ( jOb+1.LE.sNy+OLy ) THEN
           IF ( OB_Iw(jOb+1,bi,bj).NE.OB_indexNone )
     &       corner = i.LE.OB_Iw(jOb+1,bi,bj)+2
          ENDIF
#endif
          IF ( jOb+1.LE.sNy+OLy .AND. .NOT.corner ) THEN
           depOb = MAX( OBSeta(i,bi,bj) + Ro_surf(i,jOb,bi,bj)
     &                - R_low(i,jOb,bi,bj), 0.5 _d 0 )
           du = SQRT( gravity/depOb )
     &        *( OBSeta(i,bi,bj) - etaN(i,jOb+1,bi,bj) )
           DO k=1,Nr
            OBSv(i,k,bi,bj) = OBSv(i,k,bi,bj)
     &                      + du*maskS(i,jOb+1,k,bi,bj)
           ENDDO
          ENDIF
         ENDIF
        ENDDO
#endif
       ENDDO
      ENDDO
#endif /* ALLOW_OBCS and NONLIN_FRSURF */
      RETURN
      END
