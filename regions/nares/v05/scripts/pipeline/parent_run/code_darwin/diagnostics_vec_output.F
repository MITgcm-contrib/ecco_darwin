C----&------------------------------------------------------------------xxxxxxx|
C                           diagnostics_vec_output
C----&------------------------------------------------------------------xxxxxxx|
C Contributors: Michael Wood, Ian Fenty, April Shin

#include "DIAGNOSTICS_VEC_OPTIONS.h"

C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|
CBOP
C !ROUTINE: DIAGNOSTICS_VEC_OUTPUT

C !INTERFACE:
      SUBROUTINE DIAGNOSTICS_VEC_OUTPUT( myTime, myIter, myThid )

C     !DESCRIPTION:
C     Extracts virtual transects and slices of global model fields.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "DIAGNOSTICS_VEC.h"

C     !INPUT/OUTPUT PARAMETERS:
      _RL myTime
      INTEGER myThid
      INTEGER myIter

C     !LOCAL VARIABLES:
      
      INTEGER i,j, k, bi, bj, fld2D, fld3D, fld
      INTEGER fldSF
      INTEGER vec_id, sf_id
      INTEGER tag
      INTEGER total_valid
      REAL*8 buffer2D((VEC_points)*(nSx*nSy))
      REAL*8 buffer3D((VEC_points)*(nSx*nSy), Nr)
      REAL*8 bufferSF((sNx*sNy)*(nSx*nSy))
      LOGICAL DIFFERENT_MULTIPLE

      INTEGER, PARAMETER :: debug = 0
      CHARACTER*(MAX_LEN_MBUF) msgBuf

#ifndef ALLOW_USE_MPI
      INTEGER, PARAMETER :: mpiMyId = 0
#endif

CEOP
C----&------------------------------------------------------------------xxxxxxx|


#ifdef ALLOW_USE_MPI
       INTEGER pid, ierror, status(MPI_STATUS_SIZE)
       common / mpi_var / status, ierror
C  All MPI objects (e.g., MPI_Datatype, MPI_Comm) are of type INTEGER in Fortran.
       INTEGER full_dim(2), sub_dim(2), start_ind(2)
#endif

       COMMON / local_var / total_valid


       if (debug .eq. 1) then
        WRITE(msgBuf,'(A)') "========================================="
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') "  BEGIN DIAGNOSTICS_VEC_OUTPUT DEBUG "
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') "========================================="
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A,F20.3)') "     myTime: ", myTime
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A,F20.3)') "     startTime_vec: ", startTime_vec
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A,F20.3)') "     endTime_vec: ", endTime_vec
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A,I10)') "     myIter: ", myIter
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A,I5)') "     mpiMyId: ",mpiMyId
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
       endif

       IF ((myTime.le.endTime_vec).AND.
     &     (myTime.ge.startTime_vec)) then

C----&------------------------------------------------------------------xxxxxxx|
C      All Processes Tasks:
C         At each time step:
C          1. Extract field values on open boundary points
C          2. Accumulate extracted field values
C         At end of averaging period:
C          1. Generate time-averaged field values
C          2. Processes > 0 send time-averaged field values to master process
C             (Convert _RL array to R8 before sending. Send R8 array)
C----&------------------------------------------------------------------xxxxxxx|

C       Initializing local arrays
        DO i=1,(VEC_points)*(nSx*nSy)
          buffer2D(i) = 0.0
          DO k=1,Nr
            buffer3D(i,k) = 0.0
          ENDDO
        ENDDO

        DO i=1,(sNx*sNy)*(nSx*nSy)
          bufferSF(i) = 0.0
        ENDDO

        CALL TIMER_START('SET_SUBFIELDS  [DV_OUTPUT]',myThid)
C       Extracting field values on open boundary points at each time step
        CALL SET_SUBFIELDS(myThid, myIter)
        CALL TIMER_STOP('SET_SUBFIELDS  [DV_OUTPUT]',myThid)

C       Now looping through each boundary mask and communicating fields between master process and child processes
C       Each child sends their piece of the mask to the parent here
C       The parent node does the receiving in VEC_MASTER_PROC_TASKS
      DO vec_id=1, nVEC_mask

      IF (vec_nTimeSteps(vec_id) .eq. 0) then
C        Initialize averages to zero
         CALL SET_TO_VEC_INI_STATE(myThid,vec_id)
      ELSE
C        Cumulate diagnostics_vec fields and number of time steps (for time average)
         CALL CUMULATE_VEC_FLD_TIME(myThid,vec_id)
      ENDIF


C----&------------------------------------------------------------------xxxxxxx|
C      Output debug notes on check if at end of averaging period:

      if (debug .eq. 1) then
         WRITE(msgBuf,'(A)') "----------------------------------------"
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,A,I5)') " CHECKING IF OUTPUT REQUESTED FOR ",
     &     "VEC MASK ",vec_id
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
         WRITE(msgBuf,'(A)') "----------------------------------------"
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,F9.3)') "     vec_nTimeSteps: ",
     &   vec_nTimeSteps(vec_id)
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,F12.3)') "     vec_avg_periods: ",
     &   vec_avg_periods(vec_id)
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,F20.3)') "     vec_time_passed: ",
     &   vec_time_passed(vec_id)
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,F9.3)') "     deltaTClock: ",deltaTClock
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,L1)') "     DIFFERENT_MULTIPLE: ",
     &   DIFFERENT_MULTIPLE(vec_avg_periods(vec_id),myTime,deltaTClock)
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,L1)') "     Output Requested: ",
     &  ((vec_nTimeSteps(vec_id) .NE. 0) .AND.
     &  DIFFERENT_MULTIPLE(vec_avg_periods(vec_id),
     &      vec_time_passed(vec_id),deltaTClock))
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

       IF (((vec_nTimeSteps(vec_id) .NE. 0)
     &  .AND.DIFFERENT_MULTIPLE(vec_avg_periods(vec_id),
     &        vec_time_passed(vec_id),deltaTClock))
     &  .AND. (myTime .le. endTime_vec)
     &  .AND. (myTime .ge. startTime_vec)) THEN

C----&------------------------------------------------------------------xxxxxxx|
C      If at end of averaging period:

C       Compute time average of all sub fields
        CALL TIMEAVE_VEC_FLD(myThid, myTime, vec_id)

#ifdef ALLOW_USE_MPI
C        When not using MPI, proc 0 already has all of the information
C           In other words, there are no other process that would be sending info

C         Processes 1 and onward send sub fields to master process
C         if they have open boundary points
         if (mpiMyId .gt. 0) then
          if (vec_numPnts_allproc(vec_id,mpiMyId+1) .ge. 1) then

           DO fld2D=1, vec_nFlds2D(vec_id)
C           Convert _RL array to R8 array
            call VEC_PASS_RL_to_R8(buffer2D,buffer3D,vec_id,2,fld2D)

            tag = vec_tags(vec_id, 1, fld2D)

C           Send R8 array to master process
            call MPI_SEND(buffer2D, (VEC_points)*(nSx*nSy), 
     &                  MPI_REAL8, 0, tag, MPI_COMM_MODEL, ierror)

           ENDDO


           DO fld3D=1, vec_nFlds3D(vec_id)
            
C           Convert _RL array to R8 array
            call VEC_PASS_RL_to_R8(buffer2D,buffer3D,vec_id,3,fld3D)

            tag = vec_tags(vec_id, 2, fld3D)

C           Send R8 array to master process
           call MPI_SEND(buffer3D, 
     &      vec_levels3D(fld3D,vec_id)*(VEC_points)*(nSx*nSy),
     &      MPI_REAL8, 0, tag, MPI_COMM_MODEL, ierror)

           ENDDO

          end if
         end if

#endif

C         Processes 0 receives all of the points
         IF (mpiMyId .eq. 0) then
          DO fld2D=1, vec_nFlds2D(vec_id)
          call VEC_MASTER_PROC_TASKS(fld2D, vec_id, 2, myThid, myTime,
     &                                myIter)
          ENDDO
          DO fld3D=1, vec_nFlds3D(vec_id)
          call VEC_MASTER_PROC_TASKS(fld3D, vec_id, 3, myThid, myTime,
     &                                myIter)
          ENDDO
         ENDIF

C      Set averages back to zero
         CALL SET_TO_VEC_INI_STATE(myThid,vec_id)

      if (vec_iters_per_file(vec_id).ne.0) then
         if (vec_time_levels(vec_id,1).eq.vec_iters_per_file(vec_id))
     &    then
            vec_time_levels(vec_id,1) = 0
            vec_time_levels(vec_id,2) = myIter+1
         endif
      endif
      vec_time_levels(vec_id,1) = vec_time_levels(vec_id,1) + 1


C     Ending if end of averaging period
      ENDIF

      vec_nTimeSteps(vec_id) = vec_nTimeSteps(vec_id) + 1
      vec_time_passed(vec_id) = vec_time_passed(vec_id) + 
     &                           deltaTClock

C       !ENDDO for vec_id
        ENDDO




C ---- Surface Masks --------------------------------------------------------------------

C       Now looping through each surface mask and communicating fields between master process and child processes
C       Each child sends their piece of the mask to the parent here
C       The parent node does the receiving in SF_MASTER_PROC_TASKS
       DO sf_id=1, nSURF_mask

      IF (sf_nTimeSteps(sf_id) .eq. 0) then
C        Initialize averages to zero
         CALL SET_TO_SF_INI_STATE(myThid,sf_id)
      ELSE
C        Cumulate diagnostics_sf fields and number of time steps (for time average)
          CALL CUMULATE_SF_FLD_TIME(myThid,sf_id)
      ENDIF

C----&------------------------------------------------------------------xxxxxxx|
C      Output debug notes on check if at end of averaging period:

      if (debug .eq. 1) then
         WRITE(msgBuf,'(A)') "----------------------------------------"
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,A,I5)') " CHECKING IF OUTPUT REQUESTED FOR SF",
     &     " MASK ", sf_id
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
         WRITE(msgBuf,'(A)') "----------------------------------------"
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,F9.3)') "     sf_nTimeSteps: ",
     &   sf_nTimeSteps(sf_id)
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,F12.3)') "     sf_avg_periods: ",
     &   sf_avg_periods(sf_id)
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,F20.3)') "     sf_time_passed: ",
     &   sf_time_passed(sf_id)
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,F9.3)') "     deltaTClock: ",deltaTClock
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,L1)') "     DIFFERENT_MULTIPLE: ",
     &   DIFFERENT_MULTIPLE(sf_avg_periods(sf_id),myTime,deltaTClock)
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,L1)') "     Output Requested: ",
     &  ((sf_nTimeSteps(sf_id) .NE. 0) .AND.
     &  DIFFERENT_MULTIPLE(sf_avg_periods(sf_id),
     &      sf_time_passed(sf_id),deltaTClock))
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

C----&------------------------------------------------------------------xxxxxxx|

C   Check if at end of averaging period
      IF (((sf_nTimeSteps(sf_id) .NE. 0).AND.
     &  DIFFERENT_MULTIPLE(sf_avg_periods(sf_id),
     &      sf_time_passed(sf_id),deltaTClock))
     &  .AND. (myTime .le. endTime_vec)
     &  .AND. (myTime .ge. startTime_vec)) THEN

C       Compute time average of all sub fields
      CALL TIMEAVE_SF_FLD(myThid, myTime, sf_id)

#ifdef ALLOW_USE_MPI

C         Processes 1 and onward send sub fields to master process
C         if they have surface points
         if (mpiMyId .gt. 0) then
          if (sf_numPnts_allproc(sf_id,mpiMyId+1) .ge. 1) then

           DO fldSF=1, sf_nFlds(sf_id)
C           Convert _RL array to R8 array
            call SF_PASS_RL_to_R8(bufferSF,sf_id,fldSF)

            tag = sf_tags(sf_id, 1, fldSF)

C           Send R8 array to master process
            call MPI_SEND(bufferSF, (sNx*sNy)*(nSx*nSy), MPI_REAL8, 
     &                   0, tag, MPI_COMM_MODEL, ierror)
           ENDDO

          end if
         end if

#endif

         IF (mpiMyId .eq. 0) then
          DO fldSF=1, sf_nFlds(sf_id)
          call SF_MASTER_PROC_TASKS(fldSF, sf_id, myThid, myTime,
     &                                myIter)
          ENDDO
C         !ENDIF for mpiMyId .eq. 0
         ENDIF

C      Set averages back to zero
         CALL SET_TO_SF_INI_STATE(myThid,sf_id)

         if (sf_iters_per_file(sf_id).ne.0) then
         if (sf_time_levels(sf_id,1).eq.sf_iters_per_file(sf_id))
     &    then
            sf_time_levels(sf_id,1) = 0
            sf_time_levels(sf_id,2) = myIter+1
         endif
         endif
         sf_time_levels(sf_id,1) = sf_time_levels(sf_id,1) + 1

C     !ENDIF for end of averaging period
        ENDIF 

C      Move on to next time step
        sf_nTimeSteps(sf_id) = sf_nTimeSteps(sf_id) + 1
        sf_time_passed(sf_id) = sf_time_passed(sf_id) + 
     &                           deltaTClock

C       !ENDDO for sf_id
        ENDDO
      
      ENDIF

       if (debug .eq. 1) then
        WRITE(msgBuf,'(A)') "========================================="
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') "  END DIAGNOSTICS_VEC_OUTPUT DEBUG "
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') "========================================="
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        endif

      RETURN

      

      END








C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|

      SUBROUTINE VEC_MASTER_PROC_TASKS(fld, vec_id, nDim, myThid, 
     &                             myTime, myIter)

C     Appending subFields to final global_vec array

         IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "DIAGNOSTICS_VEC.h"

C     Input Parameters:
C     fld :: current field number
C     vec_id :: current mask number
C     nDim :: number of dimensions of field to output
C     myThid :: thread number
C     myTime :: current time in model
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER fld, vec_id, nDim, myThid, myIter
      _RL myTime
      REAL*8 buffer2D((VEC_points)*(nSx*nSy))
      REAL*8 buffer3D((VEC_points)*(nSx*nSy), Nr)
      INTEGER total_valid, i, k, count
      INTEGER, PARAMETER :: debug = 0
#ifndef ALLOW_USE_MPI
      INTEGER, PARAMETER :: mpiMyId = 0   
#else
      INTEGER pid, status(MPI_STATUS_SIZE), ierror, tag
      common / mpi_var / status, ierror
#endif
      common / local_var / total_valid
      



      if (debug .eq. 1) then
        WRITE(msgBuf,'(A)') "========================================="
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') "  BEGIN VEC_MASTER_PROC_TASKS DEBUG "
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') "========================================="
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        endif

        DO i=1,(VEC_points)*(nSx*nSy)
         buffer2D(i) = 0.0
         DO k=1,Nr
           buffer3D(i,k) = 0.0
         ENDDO
       ENDDO
 
       total_valid = 0
 
C-----Append Local Array--------------------------------------------------------
C     - First append process 0 local arrays to final global array if number of
C          valid elements in those local arrays > 0
C     - Indices referring to process IDs begin with 1 in vec_numPnts_allproc
C-------------------------------------------------------------------------------

      
      

      IF (vec_numPnts_allproc(vec_id,1).ge.1) then
       if (nDim .eq. 2) then

      if (debug .eq. 1) then
      WRITE(msgBuf,'(A,I5,A,A)') "mask number:",vec_id, ", field: ",
     &            vec_flds2D(fld,vec_id)
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

C       Append process 0's 2D time-averaged field value on VEC points to final array
        call VEC_PASS_RL_to_R8(buffer2D, buffer3D, vec_id, 2,fld)
        call APPEND_VEC_VALUE2D(vec_id, 0, buffer2D, myThid)

C       Add onto total number of VEC points
        total_valid = total_valid + vec_numPnts_allproc(vec_id, 1)

      if (debug .eq. 1) then
      WRITE(msgBuf,'(A,A,I5,A)') "   appending local field to global ",
     &          "array ( ",vec_numPnts_allproc(vec_id, 1),
     &          " valid elements)"
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

       else if (nDim .eq. 3) then

      if (debug .eq. 1) then
      WRITE(msgBuf,'(A,I5,A,A)') "mask number:",vec_id, ", field: ",
     &            vec_flds3D(fld,vec_id)
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

C       Append process 0's 3D time-averaged field value on VEC points to final array
        call VEC_PASS_RL_to_R8(buffer2D, buffer3D, vec_id, 3,fld)
        call APPEND_VEC_VALUE3D(vec_id, 0, fld, buffer3D, myThid)

C       Add onto total number of VEC points
        total_valid = total_valid + vec_numPnts_allproc(vec_id, 1)

      if (debug .eq. 1) then
      WRITE(msgBuf,'(A,A,I5,A)') "   appending local field to global ",
     &          "array ( ",vec_numPnts_allproc(vec_id, 1),
     &          " valid elements)"
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

      endif

      ENDIF

C-----End of appending local array----------------------------------------------


#ifdef ALLOW_USE_MPI
C------Receiving arrays through MPI and appending-------------------------------
C      - Receive all sub arrays from processes with pid>0 and valid elements>0.
C      - Append receive sub arrays to final array.
C-------------------------------------------------------------------------------

       DO pid=1, mpiNprocs-1

C       If number of valid elements is greater than or equal to 1
        IF (vec_numPnts_allproc(vec_id,pid+1) .ge. 1) then

          if (debug .eq. 1) then
            WRITE(msgBuf,'(A,I5,A,I5,A,I5)') "receiving ",
     &       vec_numPnts_allproc(vec_id,pid+1), " elements from proc ",
     &       pid, " on mask number",vec_id
            CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
          endif

         IF (nDim .eq. 2) then

          tag = vec_tags(vec_id, 1, fld)

C         MPI_RECV params: received_data, receive_count, receive_type, sender_ID, tag, comm, status, ierror
C         Master process receiving time-averaged 2D field values on VEC points from process pid
          call MPI_RECV(buffer2D, (VEC_points)*(nSx*nSy), MPI_REAL8,
     &                   pid, tag ,MPI_COMM_MODEL, status, ierror)

C         Appending the received time-averaged 2D field values on VEC points from process pid
          call APPEND_VEC_VALUE2D(vec_id, pid, buffer2D, myThid)

         ELSE IF (nDim .eq. 3) then

          tag = vec_tags(vec_id, 2, fld)

C         MPI_RECV params: received_data, receive_count, receive_type, sender_ID, tag, comm, status, ierror
C         Master process receiving time-averaged 3D field values on VEC points from process pid
          call MPI_RECV(buffer3D(:,1:vec_levels3D(fld,vec_id)),
     &           vec_levels3D(fld,vec_id)*(VEC_points)*(nSx*nSy), 
     &          MPI_REAL8, pid, tag ,MPI_COMM_MODEL, status, ierror)

C         append_vec_values3D params: vec_id, pid, fld, vec_values
          call APPEND_VEC_VALUE3D(vec_id,pid, fld, buffer3D, myThid)

         ENDIF
C        Keeping count of total number of VEC points for each mask
         total_valid = total_valid + vec_numPnts_allproc(vec_id, pid+1)

C       ENDIF for number of valid elements
        ENDIF

C      ENDDO for pid
       ENDDO
C------End of receiving arrays through MPI and appending------------------------
#endif

      if (debug .eq. 1) then
         WRITE(msgBuf,'(A,I5,A)') "   Sending ",
     &       total_valid, " elements to file "
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

C      Writing final output to binary file
       call WRITE_GLOBAL_BIN_VEC(total_valid,fld,vec_id,myIter,nDim)


       if (debug .eq. 1) then
        WRITE(msgBuf,'(A)') "========================================="
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') "  END VEC_MASTER_PROC_TASKS DEBUG "
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') "========================================="
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        endif

      END
C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|











C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|

      SUBROUTINE SF_MASTER_PROC_TASKS(fld, sf_id, myThid, myTime,
     &                             myIter)

C     Appending subFields to final global_vec array

         IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "DIAGNOSTICS_VEC.h"

C     Input Parameters:
C     fld :: current field number
C     sf_id :: current mask number
C     nDim :: number of dimensions of field to output
C     myThid :: thread number
C     myTime :: current time in model
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER fld, sf_id, myThid, myIter
      _RL myTime
      REAL*8 bufferSF((sNx*sNy)*(nSx*nSy))
      INTEGER total_valid, i, count
      INTEGER, PARAMETER :: debug = 0

#ifndef ALLOW_USE_MPI
      INTEGER, PARAMETER :: mpiMyId = 0   
#else
      INTEGER pid, status(MPI_STATUS_SIZE), ierror, tag
      common / mpi_var / status, ierror
#endif
      common / local_var / total_valid

      if (debug .eq. 1) then
      WRITE(msgBuf,'(A)') "========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "  BEGIN SF_MASTER_PROC_TASKS DEBUG "
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

      DO i=1,(sNx*sNy)*(nSx*nSy)
        bufferSF(i) = 0.0
      ENDDO

      total_valid = 0

C-----Append Local Array--------------------------------------------------------
C     - First append process 0 local arrays to final global array if number of
C          valid elements in those local arrays > 0
C     - Indices referring to process IDs begin with 1 in sf_numPnts_allproc
C-------------------------------------------------------------------------------

      if (debug .eq. 1) then
      WRITE(msgBuf,'(A,I5,A,A)') "mask number:", sf_id, ", field: ",
     &            sf_flds(fld,sf_id)
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

      IF (sf_numPnts_allproc(sf_id,1).ge.1) then

        call SF_PASS_RL_to_R8(bufferSF, sf_id, fld)
C       Append process 0's SF time-averaged field value on SF points to final array
        call APPEND_SF_VALUE(sf_id, 0, bufferSF, myThid)
C       Add onto total number of SF points
        total_valid = total_valid + sf_numPnts_allproc(sf_id, 1)

        if (debug .eq. 1) then
      WRITE(msgBuf,'(A,A,I5,A)') "   appending local field to global ",
     &          "array ( ",sf_numPnts_allproc(sf_id, 1),
     &          " valid elements)"
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
       endif

      ENDIF

C-----End of appending local array----------------------------------------------


#ifdef ALLOW_USE_MPI
C------Receiving arrays through MPI and appending-------------------------------
C      - Receive all sub arrays from processes with pid>0 and valid elements>0.
C      - Append receive sub arrays to final array.
C-------------------------------------------------------------------------------

       DO pid=1, mpiNprocs-1
C       If number of valid elements is greater than or equal to 1
        IF (sf_numPnts_allproc(sf_id,pid+1) .ge. 1) then

         if (debug .eq. 1) then
            WRITE(msgBuf,'(A,I5,A,I5,A,I5)') "   receiving ",
     &       sf_numPnts_allproc(sf_id,pid+1), " elements from proc ",
     &       pid+1, " on mask number",sf_id
            CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
          endif

          tag = sf_tags(sf_id, 1, fld)
C         !MPI_RECV params: received_data, receive_count, receive_type, sender_ID, tag, comm, status, ierror
C         Master process receiving time-averaged surface field values on SF points from process pid
          call MPI_RECV(bufferSF, (sNx*sNy)*(nSx*nSy), MPI_REAL8,
     &                   pid, tag ,MPI_COMM_MODEL, status, ierror)

C         Appending the received time-averaged SF field values on SF points from process pid
          call APPEND_SF_VALUE(sf_id, pid, bufferSF, myThid)

C        Keeping count of total number of SF points for each mask
         total_valid = total_valid + sf_numPnts_allproc(sf_id, pid+1)

C       !ENDIF for number of valid elements
        ENDIF

C      !ENDDO for pid
       ENDDO
C------End of receiving arrays through MPI and appending------------------------
#endif

      if (debug .eq. 1) then
         WRITE(msgBuf,'(A,I5,A)') "   Sending ",
     &       total_valid, " elements to file "
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

C      Writing final output to binary file
       call WRITE_GLOBAL_BIN_SF(total_valid,fld,sf_id,myIter)

       if (debug .eq. 1) then
      WRITE(msgBuf,'(A)') "========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "  END SF_MASTER_PROC_TASKS DEBUG "
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
         endif
      END
C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|









C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|
      SUBROUTINE APPEND_VEC_VALUE2D(vec_id, pid, vec_values, myThid)
C       Appends field value corresponding to open boundary points to final global array global_vec.


      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "DIAGNOSTICS_VEC.h"

         INTEGER myThid
         REAL*8 vec_values((VEC_points)*(nSx*nSy))
         INTEGER num_valid, i, vec_index, vec_id, pid
         INTEGER, PARAMETER :: debug = 0
         
         
C------------------------------------------------------------------------------|
C------------------------------------------------------------------------------|

         num_valid = vec_numPnts_allproc(vec_id,pid+1)

         do i=1, num_valid
C          Extracting vec_index from list for VEC point i in process pid
           vec_index = vec_mask_ind_list(vec_id,pid+1,i)

C          Inserting time-averaged field value on VEC point i in index vec_index of the final array, global_vec2D
           global_vec2D(vec_index) = vec_values(i)

         end do
         

       END
C------------------------------------------------------------------------------|
C------------------------------------------------------------------------------|








C------------------------------------------------------------------------------|
C------------------------------------------------------------------------------|
       SUBROUTINE APPEND_VEC_VALUE3D(vec_id, pid, fld, vec_values,
     &                               myThid)
C       Appends field value corresponding to open boundary points to final global array global_vec.

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "DIAGNOSTICS_VEC.h"

        INTEGER myThid
        REAL*8 vec_values((VEC_points)*(nSx*nSy), Nr)
        INTEGER num_valid, i, k, vec_index, vec_id, pid, fld
        INTEGER, PARAMETER :: debug = 0
        CHARACTER*(MAX_LEN_MBUF) msgBuf

C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|
        num_valid = vec_numPnts_allproc(vec_id,pid+1)

        if ((debug .eq. 1).and.(vec_id .eq. 1)) then
         WRITE(msgBuf,'(A,I5)') "From proc = ",pid
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
       endif

        do i=1, num_valid
C        Extracting vec_index from list for VEC point i in process pid
         vec_index = vec_mask_ind_list(vec_id,pid+1,i)

         if ((debug .eq. 1).and.(vec_id .eq. 1)) then
            WRITE(msgBuf,'(A,F10.3,A,I5,A,I5,A)') "  Passing theta = ",
     &      vec_values(i,1)," at location ",i,
     &      " from fldOnMsk_3Davg to location ",vec_index, " on mask 1"
            CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
          endif

C        Inserting time-averaged field value on VEC point i in index vec_index of
C        the final array, global_vec3D for each depth
         do k=1, vec_levels3D(fld,vec_id)
          global_vec3D(vec_index,k) = vec_values(i,k)
         enddo
       enddo

       END
C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|





C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|
      SUBROUTINE APPEND_SF_VALUE(sf_id, pid, sf_values,myThid)
C       Appends field value corresponding to surface boundary points to 
C       final global array global_vec.
      
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "DIAGNOSTICS_VEC.h"

         INTEGER myThid
         REAL*8 sf_values((sNx*sNy)*(nSx*nSy))
         INTEGER num_valid, i, sf_index, sf_id, pid
         INTEGER, PARAMETER :: debug = 0
         CHARACTER*(MAX_LEN_MBUF) msgBuf
C------------------------------------------------------------------------------|
C------------------------------------------------------------------------------|

         num_valid = sf_numPnts_allproc(sf_id,pid+1)

         if ((debug .eq. 1).and.(sf_id .eq. 1)) then
            WRITE(msgBuf,'(A,I5)') "From proc = ",pid
            CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
          endif

         do i=1, num_valid
C          Extracting vec_index from list for SF point i in process pid
           sf_index = sf_mask_ind_list(sf_id,pid+1,i)

           if ((debug .eq. 1).and.(sf_id .eq. 1).and.(i.ge.830)) then
            WRITE(msgBuf,'(A,F10.3,A,I5,A,I5,A)') "  Passing surf = ",
     &      sf_values(i)," at location ",i,
     &      " from fldOnMsk_SFavg to location ",sf_index, " on mask 1"
            CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
          endif

C          Inserting time-averaged field value on SF point i in index sf_index of the final array, global_SF
           global_SF(sf_index) = sf_values(i)
         end do


       END
C------------------------------------------------------------------------------|
C------------------------------------------------------------------------------|








C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|
       SUBROUTINE VEC_PASS_RL_to_R8(buffer2D,buffer3D,vec_id,nDim,fld)

         IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DIAGNOSTICS_VEC.h"

         REAL*8 buffer2D((VEC_points)*(nSx*nSy))
         REAL*8 buffer3D((VEC_points)*(nSx*nSy), Nr)
         INTEGER vec_id, nDim, fld
         INTEGER i,k
         INTEGER, PARAMETER :: debug = 0
C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|
         IF (nDim .eq. 2) then
           DO i=1, (VEC_points)*(nSx*nSy)
             buffer2D(i) = fldOnMsk_2Davg(vec_id,fld,i)
           ENDDO

         ELSE IF (nDim .eq. 3) then
           DO i=1, (VEC_points)*(nSx*nSy)
             DO k=1,vec_levels3D(fld,vec_id)
               buffer3D(i,k) = fldOnMsk_3Davg(vec_id,fld,i,k)
             ENDDO
           ENDDO

         ENDIF

       END
C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|






C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|
       SUBROUTINE SF_PASS_RL_to_R8(bufferSF, sf_id,fld)

         IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DIAGNOSTICS_VEC.h"

         REAL*8 bufferSF((sNx*sNy)*(nSx*nSy))
         INTEGER sf_id, fld
         INTEGER i
         INTEGER, PARAMETER :: debug = 0
C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|
           DO i=1, (sNx*sNy)*(nSx*nSy)
             bufferSF(i) = fldOnMsk_SFavg(sf_id,fld,i)
           ENDDO

       END
C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|





