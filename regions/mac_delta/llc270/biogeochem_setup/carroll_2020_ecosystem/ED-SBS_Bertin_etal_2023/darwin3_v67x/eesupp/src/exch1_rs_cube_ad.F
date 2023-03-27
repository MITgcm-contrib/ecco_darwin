#include "CPP_EEOPTIONS.h"

C--  File exch1_RS_cube_ad.F: adjoint routines to exch1 cube routines
C--   Contents
C--   o EXCH1_RS_CUBE_AD
C--   o EXCH1_Z_RS_CUBE_AD
C--   o EXCH1_UV_RS_CUBE_AD
C--   o EXCH1_BG_RS_CUBE_AD

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: EXCH1_RS_CUBE_AD

C     !INTERFACE:
#ifdef AUTODIFF_TAMC_COMPATIBILITY
      SUBROUTINE EXCH1_RS_CUBE_AD(
     I                 withSigns,
     I                 myOLw, myOLe, myOLs, myOLn, myNz,
     I                 exchWidthX, exchWidthY,
     I                 cornerMode, myThid,
     U                 array )
#else
      SUBROUTINE EXCH1_RS_CUBE_AD(
     U                 array,
     I                 withSigns,
     I                 myOLw, myOLe, myOLs, myOLn, myNz,
     I                 exchWidthX, exchWidthY,
     I                 cornerMode, myThid )
#endif

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE EXCH1_RS_CUBE_AD
C     | o Reverse-mode edge exchanges for RS array on CS config.
C     *==========================================================*
C     |   - not yet coded -
C     *==========================================================*

C     !USES:
      IMPLICIT NONE

C     == Global data ==
#include "SIZE.h"
#include "EEPARAMS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     array       :: Array with edges to exchange.
C     withSigns   :: Flag controlling whether field sign depends on orientation
C                 :: (signOption not yet implemented but needed for SM exch)
C     myOLw,myOLe :: West  and East  overlap region sizes.
C     myOLs,myOLn :: South and North overlap region sizes.
C     exchWidthX  :: Width of data region exchanged in X.
C     exchWidthY  :: Width of data region exchanged in Y.
C                    Note --
C                    1. In theory one could have a send width and
C                    a receive width for each face of each tile. The only
C                    restriction would be that the send width of one
C                    face should equal the receive width of the sent to
C                    tile face. Dont know if this would be useful. I
C                    have left it out for now as it requires additional
C                    bookeeping.
C     cornerMode  :: Flag indicating whether corner updates are needed.
C     myThid      :: my Thread Id number

      INTEGER myOLw, myOLe, myOLs, myOLn, myNz
      _RS     array( 1-myOLw:sNx+myOLe,
     &               1-myOLs:sNy+myOLn,
     &               myNz, nSx, nSy )
      LOGICAL withSigns
      INTEGER exchWidthX
      INTEGER exchWidthY
      INTEGER cornerMode
      INTEGER myThid

C     !LOCAL VARIABLES:
C     == Local variables ==
C     theSimulationMode :: Holds working copy of simulation mode
C     theCornerMode     :: Holds working copy of corner mode
c     INTEGER theSimulationMode
c     INTEGER theCornerMode
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

c     theSimulationMode = REVERSE_SIMULATION
c     theCornerMode     = cornerMode

c     IF ( simulationMode.EQ.REVERSE_SIMULATION ) THEN
      WRITE(msgBuf,'(A)') 'ADEXCH1_RS_CUBE: AD mode not implemented'
      CALL PRINT_ERROR( msgBuf, myThid )
      STOP 'ABNORMAL END: EXCH1_RS_CUBE_AD: no AD code'
c     ENDIF

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP

C     !ROUTINE: EXCH1_Z_RS_CUBE_AD

C     !INTERFACE:
#ifdef AUTODIFF_TAMC_COMPATIBILITY
      SUBROUTINE EXCH1_Z_RS_CUBE_AD(
     I                 withSigns,
     I                 myOLw, myOLe, myOLs, myOLn, myNz,
     I                 exchWidthX, exchWidthY,
     I                 cornerMode, myThid,
     U                 array )
#else
      SUBROUTINE EXCH1_Z_RS_CUBE_AD(
     U                 array,
     I                 withSigns,
     I                 myOLw, myOLe, myOLs, myOLn, myNz,
     I                 exchWidthX, exchWidthY,
     I                 cornerMode, myThid )
#endif

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE EXCH1_Z_RS_CUBE_AD
C     | o Reverse-mode edge exchanges for RS array on CS config:
C     |   Fill overlap region through tile exchanges,
C     |   according to CS topology,
C     |   for a Zeta-located, scalar field RS arrays.
C     *==========================================================*
C     |   - not yet coded -
C     *==========================================================*

C     !USES:
      IMPLICIT NONE

C     == Global data ==
#include "SIZE.h"
#include "EEPARAMS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     array       :: Array with edges to exchange.
C     withSigns   :: Flag controlling whether field sign depends on orientation
C                 :: (signOption not yet implemented but needed for SM exch)
C     myOLw,myOLe :: West  and East  overlap region sizes.
C     myOLs,myOLn :: South and North overlap region sizes.
C     exchWidthX  :: Width of data region exchanged in X.
C     exchWidthY  :: Width of data region exchanged in Y.
C                    Note --
C                    1. In theory one could have a send width and
C                    a receive width for each face of each tile. The only
C                    restriction would be that the send width of one
C                    face should equal the receive width of the sent to
C                    tile face. Dont know if this would be useful. I
C                    have left it out for now as it requires additional
C                    bookeeping.
C     cornerMode  :: Flag indicating whether corner updates are needed.
C     myThid      :: my Thread Id number

      INTEGER myOLw, myOLe, myOLs, myOLn, myNz
      _RS     array( 1-myOLw:sNx+myOLe,
     &               1-myOLs:sNy+myOLn,
     &               myNz, nSx, nSy )
      LOGICAL withSigns
      INTEGER exchWidthX
      INTEGER exchWidthY
      INTEGER cornerMode
      INTEGER myThid

C     !LOCAL VARIABLES:
C     == Local variables ==
C     theSimulationMode :: Holds working copy of simulation mode
C     theCornerMode     :: Holds working copy of corner mode
c     INTEGER theSimulationMode
c     INTEGER theCornerMode
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

c     theSimulationMode = REVERSE_SIMULATION
c     theCornerMode     = cornerMode

c     IF ( simulationMode.EQ.REVERSE_SIMULATION ) THEN
      WRITE(msgBuf,'(A)') 'ADEXCH1_Z_RS_CUBE: AD mode not implemented'
      CALL PRINT_ERROR( msgBuf, myThid )
      STOP 'ABNORMAL END: EXCH1_Z_RS_CUBE_AD: no AD code'
c     ENDIF

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: EXCH1_UV_RS_CUBE_AD

C     !INTERFACE:
#ifdef AUTODIFF_TAMC_COMPATIBILITY
      SUBROUTINE EXCH1_UV_RS_CUBE_AD(
     I                 withSigns,
     I                 myOLw, myOLe, myOLs, myOLn, myNz,
     I                 exchWidthX, exchWidthY,
     I                 cornerMode, myThid,
     U                 Uarray, Varray )
#else
      SUBROUTINE EXCH1_UV_RS_CUBE_AD(
     U                 Uarray, Varray,
     I                 withSigns,
     I                 myOLw, myOLe, myOLs, myOLn, myNz,
     I                 exchWidthX, exchWidthY,
     I                 cornerMode, myThid )
#endif

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE EXCH1_UV_RS_CUBE_AD
C     | o Reverse-mode edge exchanges for RS vector on CS config.
C     *==========================================================*
C     |   - not yet coded -
C     *==========================================================*

C     !USES:
      IMPLICIT NONE

C     == Global data ==
#include "SIZE.h"
#include "EEPARAMS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     Uarray      :: (u-type) Array with edges to exchange.
C     Varray      :: (v-type) Array with edges to exchange.
C     withSigns   :: sign of Uarray,Varray depends on orientation
C     myOLw,myOLe :: West  and East  overlap region sizes.
C     myOLs,myOLn :: South and North overlap region sizes.
C     exchWidthX  :: Width of data region exchanged in X.
C     exchWidthY  :: Width of data region exchanged in Y.
C                    Note --
C                    1. In theory one could have a send width and
C                    a receive width for each face of each tile. The only
C                    restriction would be that the send width of one
C                    face should equal the receive width of the sent to
C                    tile face. Dont know if this would be useful. I
C                    have left it out for now as it requires additional
C                    bookeeping.
C     cornerMode  :: Flag indicating whether corner updates are needed.
C     myThid      :: my Thread Id number

      INTEGER myOLw, myOLe, myOLs, myOLn, myNz
      _RS     Uarray( 1-myOLw:sNx+myOLe,
     &                1-myOLs:sNy+myOLn,
     &                myNZ, nSx, nSy )
      _RS     Varray( 1-myOLw:sNx+myOLe,
     &                1-myOLs:sNy+myOLn,
     &                myNZ, nSx, nSy )
      LOGICAL withSigns
      INTEGER exchWidthX
      INTEGER exchWidthY
      INTEGER cornerMode
      INTEGER myThid

C     !LOCAL VARIABLES:
C     == Local variables ==
C     theSimulationMode :: Holds working copy of simulation mode
C     theCornerMode     :: Holds working copy of corner mode
c     INTEGER theSimulationMode
c     INTEGER theCornerMode
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

c     theSimulationMode = REVERSE_SIMULATION
c     theCornerMode     = cornerMode

c     IF ( simulationMode.EQ.REVERSE_SIMULATION ) THEN
      WRITE(msgBuf,'(A)') 'ADEXCH1_UV_RS_CUBE: AD mode not implemented'
      CALL PRINT_ERROR( msgBuf, myThid )
      STOP 'ABNORMAL END: EXCH1_UV_RS_CUBE_AD: no AD code'
c     ENDIF

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP

C     !ROUTINE: EXCH1_BG_RS_CUBE_AD

C     !INTERFACE:
#ifdef AUTODIFF_TAMC_COMPATIBILITY
      SUBROUTINE EXCH1_BG_RS_CUBE_AD(
     I                 withSigns,
     I                 myOLw, myOLe, myOLs, myOLn, myNz,
     I                 exchWidthX, exchWidthY,
     I                 cornerMode, myThid,
     U                 uField, vField )
#else
      SUBROUTINE EXCH1_BG_RS_CUBE_AD(
     U                 uField, vField,
     I                 withSigns,
     I                 myOLw, myOLe, myOLs, myOLn, myNz,
     I                 exchWidthX, exchWidthY,
     I                 cornerMode, myThid )
#endif

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE EXCH1_BG_RS_CUBE_AD
C     | o Forward-mode edge exchanges for RS vector on CS config:
C     |   Fill overlap region through tile exchanges,
C     |   according to CS topology,
C     |   for a 2-Components B-Grid vector field RS arrays.
C     *==========================================================*
C     |   - not yet coded -
C     *==========================================================*

C     !USES:
      IMPLICIT NONE

C     == Global data ==
#include "SIZE.h"
#include "EEPARAMS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     uField     :: 1rst component array with overlap to exchange.
C     vField     :: 2nd  component array with overlap to exchange.
C     withSigns  :: sign of uField,vField depends on orientation.
C     myOLw,myOLe :: West  and East  overlap region sizes.
C     myOLs,myOLn :: South and North overlap region sizes.
C     exchWidthX  :: Width of data region exchanged in X.
C     exchWidthY  :: Width of data region exchanged in Y.
C                    Note --
C                    1. In theory one could have a send width and
C                    a receive width for each face of each tile. The only
C                    restriction would be that the send width of one
C                    face should equal the receive width of the sent to
C                    tile face. Dont know if this would be useful. I
C                    have left it out for now as it requires additional
C                    bookeeping.
C     cornerMode  :: Flag indicating whether corner updates are needed.
C     myThid      :: my Thread Id number

      INTEGER myOLw, myOLe, myOLs, myOLn, myNz
      _RS uField( 1-myOLw:sNx+myOLe, 1-myOLs:sNy+myOLn,
     &            myNz, nSx, nSy )
      _RS vField( 1-myOLw:sNx+myOLe, 1-myOLs:sNy+myOLn,
     &            myNz, nSx, nSy )
      LOGICAL withSigns
      INTEGER exchWidthX
      INTEGER exchWidthY
      INTEGER cornerMode
      INTEGER myThid

C     !LOCAL VARIABLES:
C     == Local variables ==
C     theSimulationMode :: Holds working copy of simulation mode
C     theCornerMode     :: Holds working copy of corner mode
c     INTEGER theSimulationMode
c     INTEGER theCornerMode
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

c     theSimulationMode = REVERSE_SIMULATION
c     theCornerMode     = cornerMode

c     IF ( simulationMode.EQ.REVERSE_SIMULATION ) THEN
      WRITE(msgBuf,'(A)') 'ADEXCH1_BG_RS_CUBE: AD mode not implemented'
      CALL PRINT_ERROR( msgBuf, myThid )
      STOP 'ABNORMAL END: EXCH1_BG_RS_CUBE_AD: no AD code'
c     ENDIF

      RETURN
      END
