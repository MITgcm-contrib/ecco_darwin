#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "DARWIN_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_DARWIN

#ifdef ALLOW_CARBON

CBOP
C !ROUTINE: DIC_ATMOS

C !INTERFACE: ==========================================================
      SUBROUTINE DIC_ATMOS( istate, myTime, myIter, myThid )

C !DESCRIPTION:
C  Calculate the atmospheric pCO2
C  dic_int1:
C  0=use default 278.d-6
C  1=use constant value - dic_pCO2, read in from data.dic
C  2=read in from file
C  3=interact with atmospheric box (use dic_pCO2 as initial atmos. value)

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "FFIELDS.h"
#include "DARWIN_SIZE.h"
#include "DARWIN_IO.h"
#include "DARWIN_FLUX.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_FIELDS.h"
#include "DIC_ATMOS.h"
#ifdef USE_EXFWIND
#include "EXF_FIELDS.h"
#endif

C !INPUT PARAMETERS: ===================================================
C  istate               :: 0=initial call, 1=subsequent calls
C  myTime               :: current time
C  myIter               :: current iteration number
C  myThid               :: my Thread Id number
      INTEGER istate
      _RL myTime
      INTEGER myIter, myThid

      LOGICAL  DIFFERENT_MULTIPLE
      EXTERNAL DIFFERENT_MULTIPLE

C !LOCAL VARIABLES: ====================================================
C   total_atmos_moles :: atmosphere total gas content (should be parameter)
      _RL total_atmos_moles
      INTEGER bi, bj, i,j,k
      INTEGER ntim, iUnit

      _RL total_flux
      _RL total_carbon

C for carbon budget ouput
      _RL total_ocean_carbon_old
      _RL total_atmos_carbon_old
      _RL atpco2
      _RL total_carbon_old, carbon_diff
      _RL year_diff_ocean, year_diff_atmos, year_total
      _RL start_diff_ocean, start_diff_atmos, start_total
C variables for reading CO2 input files
      _RL tmp
      _RL aWght, bWght

       CHARACTER*(MAX_LEN_FNAM) fn
       LOGICAL permCheckPoint
CEOP

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

c if coupled to atmsopheric model, use the
c Co2 value passed from the coupler
#ifndef USE_ATMOSCO2

      IF ( nThreads .GT. 1 .AND.
     &     ( dic_int1.EQ.2 .OR. dic_int1.EQ.3 ) ) THEN
C     Problem with I/O and global-sum for multi-threaded execution
C     Needs to be fixed before using this S/R in multi-threaded run
        STOP 'S/R DIC_ATMOS: multi-threaded not right'
      ENDIF
# ifdef USE_EXFCO2
          DO bj=myByLo(myThid),myByHi(myThid)
           DO bi=myBxLo(myThid),myBxHi(myThid)
            DO j=1-OLy,sNy+OLy
             DO i=1-OLx,sNx+OLx
               AtmospCO2(i,j,bi,bj)=apco2(i,j,bi,bj)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
# else

c default - set only once
      if (dic_int1.eq.0.and.istate.eq.0) then
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)

         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
             AtmospCO2(i,j,bi,bj)=278.0 _d -6
          ENDDO
         ENDDO

        ENDDO
       ENDDO
      endif

C user specified value (or default = 278 ppm)- set only once
      IF (dic_int1.EQ.1 .AND. istate.EQ.0 ) THEN
        DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
             AtmospCO2(i,j,bi,bj)=dic_pCO2
           ENDDO
          ENDDO
         ENDDO
        ENDDO
      ENDIF

C read from a file (note:
C                   dic_int2=number entries to read
C                   dic_int3=start timestep,
C                   dic_int4=timestep between file entries)
      IF (dic_int1.EQ.2) THEN
C linearly interpolate between file entries
        ntim=int((myIter-dic_int3)/dic_int4)+1
        aWght = FLOAT(myIter-dic_int3)
        bWght = FLOAT(dic_int4)
        aWght = 0.5 _d 0 + aWght/bWght - FLOAT(ntim-1)
        IF (aWght.GT.1. _d 0) THEN
          ntim=ntim+1
          aWght=aWght-1. _d 0
        ENDIF
        bWght = 1. _d 0 - aWght
        tmp=co2atmos(ntim)*bWght+co2atmos(ntim+1)*aWght
c       print*,'weights',ntim, aWght, bWght, tmp

        DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
           DO j=1-OLy,sNy+OLy
             DO i=1-OLx,sNx+OLx
               AtmospCO2(i,j,bi,bj)=tmp
             ENDDO
           ENDDO
         ENDDO
        ENDDO

      ENDIF

C interactive atmosphere
      IF (dic_int1.EQ.3) THEN

c      _BEGIN_MASTER(myThid)

C Mass dry atmosphere = (5.1352+/-0.0003)d18 kg (Trenberth & Smith,
C Journal of Climate 2005)
C and Mean molecular mass air = 28.97 g/mol (NASA earth fact sheet)
       total_atmos_moles= 1.77 _d 20
C for 278ppmv we need total_atmos_carbon=4.9206e+16

       if (istate.gt.0) then
        total_ocean_carbon_old=total_ocean_carbon
        total_atmos_carbon_old=total_atmos_carbon
       else
        total_ocean_carbon_old=0. _d 0
        total_atmos_carbon_old=0. _d 0
       endif

       total_flux= 0. _d 0
       total_ocean_carbon= 0. _d 0

       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
          DO j=1,sNy
           DO i=1,sNx
            IF (istate.gt.0) then
              total_flux=total_flux
     &                        + FluxCO2(i,j,bi,bj)*rA(i,j,bi,bj)
     &                         *maskC(i,j,1,bi,bj)*dTtracerLev(1)
            ENDIF
          ENDDO
         ENDDO
        ENDDO
       ENDDO

       _GLOBAL_SUM_RL(total_flux,myThid)
       _GLOBAL_SUM_RL(total_ocean_carbon,myThid)

       IF (istate.EQ.0) THEN
C use value read in dic_init_fixed
          total_atmos_carbon=total_atmos_carbon_ini
       ELSE
C calculate new atmos pCO2
          total_atmos_carbon=total_atmos_carbon - total_flux
C write out if time for a new pickup
          permCheckPoint = .FALSE.
          permCheckPoint =
     &      DIFFERENT_MULTIPLE(pChkptFreq,myTime,deltaTClock)
          IF (permCheckPoint) THEN
            DO i = 1,MAX_LEN_FNAM
                fn(i:i) = ' '
            ENDDO
            WRITE(fn,'(A,I10.10)') 'dic_atmos.',myIter
C     Going to really do some IO. Make everyone except master thread wait.
            _BARRIER
C write values to new pickup

            CALL MDSFINDUNIT( iUnit, myThid )
            open(iUnit,file=fn,status='new')
            write(iUnit,*) total_atmos_carbon, atpco2
            close(iUnit)

        ENDIF
       ENDIF

       atpco2=total_atmos_carbon/total_atmos_moles

c     print*,'QQpCO2', total_atmos_carbon, atpco2, total_ocean_carbon,
c    &                 total_flux

C--    Set AtmospCO2 for next iteration:
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
            AtmospCO2(i,j,bi,bj) = atpco2
          ENDDO
         ENDDO
        ENDDO
       ENDDO

      print*,'QQ atmos C, total, pCo2', total_atmos_carbon, atpco2
      total_carbon=total_atmos_carbon + total_ocean_carbon
      total_carbon_old=total_atmos_carbon_old + total_ocean_carbon_old
      carbon_diff=total_carbon-total_carbon_old
      print*,'QQ total C, current, old, diff', total_carbon,
     &                         total_carbon_old, carbon_diff
      carbon_diff=total_ocean_carbon-total_ocean_carbon_old
      tmp=carbon_diff-total_flux
      print*,'QQ ocean C, current, old, diff',total_ocean_carbon,
     &                   total_ocean_carbon_old, carbon_diff
      print*,'QQ air-sea flux, addition diff', total_flux, tmp

c if end of forcing cycle, find total change in ocean carbon
      if (istate.eq.0) then
       total_ocean_carbon_start=total_ocean_carbon
       total_ocean_carbon_year=total_ocean_carbon
       total_atmos_carbon_start=total_atmos_carbon
       total_atmos_carbon_year=total_atmos_carbon
      else
        permCheckPoint = .FALSE.
        permCheckPoint =
     &      DIFFERENT_MULTIPLE(externForcingCycle,myTime,deltaTClock)
        if (permCheckPoint) then
          year_diff_ocean=total_ocean_carbon-total_ocean_carbon_year
          year_diff_atmos=total_atmos_carbon-total_atmos_carbon_year
          year_total=(total_ocean_carbon+total_atmos_carbon) -
     &               (total_ocean_carbon_year+total_atmos_carbon_year)
          start_diff_ocean=total_ocean_carbon-total_ocean_carbon_start
          start_diff_atmos=total_atmos_carbon-total_atmos_carbon_start
          start_total=(total_ocean_carbon+total_atmos_carbon) -
     &            (total_ocean_carbon_start+total_atmos_carbon_start)
          print*,'QQ YEAR END'
          print*,'year diff: ocean, atmos, total', year_diff_ocean,
     &                year_diff_atmos, year_total
          print*,'start diff: ocean, atmos, total ', start_diff_ocean,
     &                start_diff_atmos, start_total
c
          total_ocean_carbon_year=total_ocean_carbon
          total_atmos_carbon_year=total_atmos_carbon
        endif
        endif

c     _END_MASTER(myThid)

      ENDIF

#endif  /*USE_EXFCO2*/

#endif

      RETURN
      END
#endif  /*ALLOW_CARBON*/

#endif  /*DARWIN*/
#endif  /*ALLOW_PTRACERS*/
c ==================================================================
