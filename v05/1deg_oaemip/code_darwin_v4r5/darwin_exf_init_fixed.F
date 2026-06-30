#include "DARWIN_OPTIONS.h"
#include "EXF_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: DARWIN_EXF_INIT_FIXED

C     !INTERFACE:
      SUBROUTINE DARWIN_EXF_INIT_FIXED( myThid )

C     !DESCRIPTION:
C     Initialize DARWIN constants for exf field loading

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "EXF_PARAM.h"
#include "EXF_INTERP_SIZE.h"
#include "DARWIN_SIZE.h"
#include "DARWIN_EXF_PARAMS.h"
#include "DARWIN_EXF_FIELDS.h"

C     !INPUT PARAMETERS:
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
C     === Local variables ===
C     msgBuf    :: Informational/error message buffer
C     errCount  :: error counter
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER errCount

#ifdef ALLOW_DARWIN
#ifdef ALLOW_EXF

      _BEGIN_MASTER( myThid )
      errCount = 0


      IF ( PARfile .NE. ' ' ) THEN
# ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('GETFFIELD_START PAR',
     &  myThid)
# endif
       CALL EXF_GETFFIELD_START( useExfYearlyFields,
     I                         'darwin', 'PAR', PARperiod,
     I                         PARstartdate1, PARstartdate2,
     U                         PARStartTime, errCount,
     I                         myThid )
      ENDIF

      IF ( ironfile .NE. ' ' ) THEN
# ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('GETFFIELD_START iron',
     &  myThid)
# endif
       CALL EXF_GETFFIELD_START( useExfYearlyFields,
     I                         'darwin', 'iron', ironperiod,
     I                         ironstartdate1, ironstartdate2,
     U                         ironStartTime, errCount,
     I                         myThid )
      ENDIF

      IF ( icefile .NE. ' ' ) THEN
# ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('GETFFIELD_START ice',
     &  myThid)
# endif
       CALL EXF_GETFFIELD_START( useExfYearlyFields,
     I                         'darwin', 'ice', iceperiod,
     I                         icestartdate1, icestartdate2,
     U                         iceStartTime, errCount,
     I                         myThid )
      ENDIF

      IF ( windfile .NE. ' ' ) THEN
# ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('GETFFIELD_START wind',
     &  myThid)
# endif
       CALL EXF_GETFFIELD_START( useExfYearlyFields,
     I                         'darwin', 'wind', windperiod,
     I                         windstartdate1, windstartdate2,
     U                         windStartTime, errCount,
     I                         myThid )
      ENDIF

      IF ( pCO2file .NE. ' ' ) THEN
# ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('GETFFIELD_START pCO2',
     &  myThid)
# endif
       CALL EXF_GETFFIELD_START( useExfYearlyFields,
     I                         'darwin', 'pCO2', pCO2period,
     I                         pCO2startdate1, pCO2startdate2,
     U                         pCO2StartTime, errCount,
     I                         myThid )
      ENDIF

      IF ( ventHe3file .NE. ' ' ) THEN
# ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('GETFFIELD_START ventHe3',
     &  myThid)
# endif
       CALL EXF_GETFFIELD_START( useExfYearlyFields,
     I                         'darwin', 'ventHe3', ventHe3period,
     I                         ventHe3startdate1, ventHe3startdate2,
     U                         ventHe3StartTime, errCount,
     I                         myThid )
      ENDIF

      IF ( DOCrunofffile .NE. ' ' ) THEN
# ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('GETFFIELD_START DOCrunoff',
     &  myThid)
# endif
       CALL EXF_GETFFIELD_START( useExfYearlyFields,
     I                         'darwin', 'DOCrunoff', DOCrunoffperiod,
     I                         DOCrunoffstartdate1, DOCrunoffstartdate2,
     U                         DOCrunoffStartTime, errCount,
     I                         myThid )
      ENDIF

      IF ( DONrunofffile .NE. ' ' ) THEN
# ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('GETFFIELD_START DONrunoff',
     &  myThid)
# endif
       CALL EXF_GETFFIELD_START( useExfYearlyFields,
     I                         'darwin', 'DONrunoff', DONrunoffperiod,
     I                         DONrunoffstartdate1, DONrunoffstartdate2,
     U                         DONrunoffStartTime, errCount,
     I                         myThid )
      ENDIF

      IF ( DOPrunofffile .NE. ' ' ) THEN
# ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('GETFFIELD_START DOPrunoff',
     &  myThid)
# endif
       CALL EXF_GETFFIELD_START( useExfYearlyFields,
     I                         'darwin', 'DOPrunoff', DOPrunoffperiod,
     I                         DOPrunoffstartdate1, DOPrunoffstartdate2,
     U                         DOPrunoffStartTime, errCount,
     I                         myThid )
      ENDIF

      IF ( DINrunofffile .NE. ' ' ) THEN
# ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('GETFFIELD_START DINrunoff',
     &  myThid)
# endif
       CALL EXF_GETFFIELD_START( useExfYearlyFields,
     I                         'darwin', 'DINrunoff', DINrunoffperiod,
     I                         DINrunoffstartdate1, DINrunoffstartdate2,
     U                         DINrunoffStartTime, errCount,
     I                         myThid )
      ENDIF

      IF ( IPrunofffile .NE. ' ' ) THEN
# ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('GETFFIELD_START IPrunoff',
     &  myThid)
# endif
       CALL EXF_GETFFIELD_START( useExfYearlyFields,
     I                         'darwin', 'IPrunoff', IPrunoffperiod,
     I                         IPrunoffstartdate1, IPrunoffstartdate2,
     U                         IPrunoffStartTime, errCount,
     I                         myThid )
      ENDIF

      IF ( DSirunofffile .NE. ' ' ) THEN
# ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('GETFFIELD_START DSirunoff',
     &  myThid)
# endif
       CALL EXF_GETFFIELD_START( useExfYearlyFields,
     I                         'darwin', 'DSirunoff', DSirunoffperiod,
     I                         DSirunoffstartdate1, DSirunoffstartdate2,
     U                         DSirunoffStartTime, errCount,
     I                         myThid )
      ENDIF

      IF ( POCrunofffile .NE. ' ' ) THEN
# ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('GETFFIELD_START POCrunoff',
     &  myThid)
# endif
       CALL EXF_GETFFIELD_START( useExfYearlyFields,
     I                         'darwin', 'POCrunoff', POCrunoffperiod,
     I                         POCrunoffstartdate1, POCrunoffstartdate2,
     U                         POCrunoffStartTime, errCount,
     I                         myThid )
      ENDIF

      IF ( POPrunofffile .NE. ' ' ) THEN
# ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('GETFFIELD_START POPrunoff',
     &  myThid)
# endif
       CALL EXF_GETFFIELD_START( useExfYearlyFields,
     I                         'darwin', 'POPrunoff', POPrunoffperiod,
     I                         POPrunoffstartdate1, POPrunoffstartdate2,
     U                         POPrunoffStartTime, errCount,
     I                         myThid )
      ENDIF

      IF ( PONrunofffile .NE. ' ' ) THEN
# ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('GETFFIELD_START PONrunoff',
     &  myThid)
# endif
       CALL EXF_GETFFIELD_START( useExfYearlyFields,
     I                         'darwin', 'PONrunoff', PONrunoffperiod,
     I                         PONrunoffstartdate1, PONrunoffstartdate2,
     U                         PONrunoffStartTime, errCount,
     I                         myThid )
      ENDIF

      IF ( DICrunofffile .NE. ' ' ) THEN
# ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('GETFFIELD_START DICrunoff',
     &  myThid)
# endif
       CALL EXF_GETFFIELD_START( useExfYearlyFields,
     I                         'darwin', 'DICrunoff', DICrunoffperiod,
     I                         DICrunoffstartdate1, DICrunoffstartdate2,
     U                         DICrunoffStartTime, errCount,
     I                         myThid )
      ENDIF

      IF ( errCount.GE.1 ) THEN
       WRITE(msgBuf,'(A,I3,A)')
     &     'DARWIN_EXF_INIT_FIXED: detected', errCount,' fatal error(s)'
       CALL PRINT_ERROR( msgBuf, myThid )
       CALL ALL_PROC_DIE( 0 )
       STOP 'ABNORMAL END: S/R DARWIN_EXF_INIT_FIXED'
      ENDIF

#ifdef DARWIN_OAE
C========================================================
C     OAE SITE to GRID CELL MAPPING (INITIALIZATION ONLY)
C========================================================
        WRITE(msgBuf,'(A)')
     &   'OAE: initializing site index mapping'
        CALL PRINT_MESSAGE(msgBuf,
     &       standardMessageUnit,
     &       SQUEEZE_RIGHT,myThid)

        DO isite = 1,nOAEsites
C========================================================
C     SEARCH OVER ALL TILES
C========================================================
         OAE_nCells(isite) = 0  
         DO bj = 1,nSy
           DO bi = 1,nSx
            DO j = 1,sNy
             DO i = 1,sNx
C--     Only ocean points
              IF ( hFacC(i,j,1,bi,bj) .GT. 0. _d 0 ) THEN
C--     Longitude wrap
                dlon = XC(i,j,bi,bj) - OAE_lon(isite)
                dlat = YC(i,j,bi,bj) - OAE_lat(isite)
                IF ( dlon .GT.  180. ) dlon = dlon - 360.
                IF ( dlon .LT. -180. ) dlon = dlon + 360.
C--     simple spherical metric (robust enough for 1° ECCO)
                dist2 =
     &            dlat*dlat + dlon*dlon
C--     keep minimum
                IF ( dist2 .LT. 1.0d20 ) THEN
                  minDist = dist2
                  OAE_i(isite)  = i
                  OAE_j(isite)  = j
                  OAE_bi(isite) = bi
                  OAE_bj(isite) = bj
                ENDIF
              ENDIF
            ENDDO
           ENDDO
          ENDDO
         ENDDO

C========================================================
C     SAFETY CHECK
C========================================================
        IF ( OAE_i(isite) .LT. 0 ) THEN
          WRITE(msgBuf,'(A,I4)')
     &     'OAE ERROR: no ocean cell found for site ', isite
          CALL PRINT_ERROR(msgBuf,myThid)
          STOP 'ABNORMAL END'
        ENDIF

C========================================================
C     PRINT RESULT
C========================================================
        WRITE(msgBuf,'(A,I4,A,2I6)')
     &    'OAE site mapped: ', isite,
     &    '  -> (i,j)=', OAE_i(isite), OAE_j(isite)

        CALL PRINT_MESSAGE(msgBuf,
     &      standardMessageUnit,
     &      SQUEEZE_RIGHT,myThid)


C========================================================
C     GAUSSIAN FOOTPRINT CONSTRUCTION
C========================================================
        WRITE(msgBuf,'(A)')
     &   'OAE: initializing gaussian mapping'
        CALL PRINT_MESSAGE(msgBuf,
     &       standardMessageUnit,
     &       SQUEEZE_RIGHT,myThid)

        icount = 0
        weight_sum = 0. _d 0

        DO bj = 1,nSy
           DO bi = 1,nSx

            DO j = 1,sNy
             DO i = 1,sNx

              dlon = XC(i,j,bi,bj) - OAE_lon(isite)
              dlat = YC(i,j,bi,bj) - OAE_lat(isite)

              IF ( dlon .GT. 180. ) dlon = dlon - 360.
              IF ( dlon .LT.-180. ) dlon = dlon + 360.

              dx = dlon * 111000. _d 0 * COS(YC(i,j,bi,bj)*deg2rad)
              dy = dlat * 111000. _d 0

              dist2 = dx*dx + dy*dy

              IF ( dist2 .LT. (3.0*OAE_sigma)**2 ) THEN
                icount = icount + 1

                IF ( icount .GT. maxOAECells ) THEN
                  CALL PRINT_ERROR(
     &              'OAE: increase maxOAECells',
     &              myThid )
                  STOP 'ABNORMAL END'
                ENDIF

                OAE_iCell(icount,isite) = i
                OAE_jCell(icount,isite) = j
                OAE_biCell(icount,isite) = bi
                OAE_bjCell(icount,isite) = bj
      
                OAE_weight(icount,isite) = EXP( -r2 /(2.0*OAE_sigma(isite)
     &                                  *OAE_sigma(isite)) )
              
              ENDIF
              IF ( hFacC(i,j,bi,bj) .EQ. 0.0 ) THEN
                 OAE_weight(icount,isite) = 0
              ENDIF

              weight_sum = weight_sum + OAE_weight(icount,isite)

                IF ( icount .GT. maxOAECells ) THEN
                  CALL PRINT_ERROR(
     &              'OAE: increase maxOAECells',
     &              myThid )
                  STOP 'ABNORMAL END'
                ENDIF

        ENDIF

        ENDDO
        ENDDO
        ENDDO
        ENDDO

        DO ic = 1, icount
            OAE_weight(ic,isite)
     &        = OAE_weight(ic,isite) / weight_sum
          ENDDO

        OAE_nCells(isite) = icount

        IF ( OAE_nCells(isite) .LE. 0 ) THEN
            CALL PRINT_ERROR(
     &        'OAE: no Gaussian cells found',
     &        myThid )
            STOP 'ABNORMAL END'
          ENDIF

          WRITE(msgBuf,'(A,I4,A,I6)')
     &      'OAE site ', isite,
     &      ' cells = ', OAE_nCells(isite)

          CALL PRINT_MESSAGE(msgBuf,
     &        standardMessageUnit,
     &        SQUEEZE_RIGHT,myThid)

C - end loop over sites
          ENDDO

#endif


      _END_MASTER( myThid )
      _BARRIER

#endif /* ALLOW_EXF */
#endif /* ALLOW_DARWIN */

      RETURN
      END
