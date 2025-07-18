C $Header: /home/ubuntu/mnt/e9_copy/MITgcm_contrib/ecco_darwin/v4_llc270/code_darwin/darwin_generate_phyto.F,v 1.3 2020/01/08 16:05:31 dcarroll Exp $
C $Name:  $

#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "DARWIN_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_DARWIN

c ==========================================================
c SUBROUTINE DARWIN_GENERATE_PHYTO
c generate parameters for "functional group" of phyto (index np)
c using a "Monte Carlo" approach
c         Mick Follows, Scott Grant Fall/Winter 2005
c         Stephanie Dutkiewicz Spring/Summer 2005
c         Anna Hickman Summer 2008
c TWO_SPECIES_SETUP
c  1=large, 2=small
c NINE_SPECIES_SETUP
c  1=diatom, 2=other large, 3=syn, 4=hl pro, 5=ll pro, 6=trich, 
c  7=uni diaz, 8=small euk, 9=cocco
c ==========================================================
        SUBROUTINE DARWIN_GENERATE_PHYTO(myThid, np)

        implicit none
#include "EEPARAMS.h"
#include "DARWIN_SIZE.h"
#include "DARWIN_PARAMS.h"
#include "DARWIN.h"


c ANNA define WAVEBANDS variables
#ifdef WAVEBANDS
#include "SPECTRAL_SIZE.h"
#include "WAVEBANDS_PARAMS.h"
#endif



C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
        INTEGER myThid

C     === Functions ===
      _RL DARWIN_RANDOM
      EXTERNAL DARWIN_RANDOM
      _RL DARWIN_RANDOM_NORMAL
      EXTERNAL DARWIN_RANDOM_NORMAL

C     !LOCAL VARIABLES:
C     === Local variables ===
C     msgBuf    - Informational/error meesage buffer
      CHARACTER*(MAX_LEN_MBUF) msgBuf

        _RL RandNo
        _RL growthdays
        _RL mortdays
        _RL pday
        _RL year
        _RL month
        _RL fiveday
        _RL rtime
        _RL standin
        _RL dm
        _RL volp
        _RL PI
        INTEGER np
        INTEGER nz
        INTEGER signvar
        PARAMETER ( PI    = 3.14159265358979323844D0   )

CEOP
c
        standin=0. _d 0

c length of day (seconds)
        pday = 86400.0 _d 0

c each time generate another functional group add one to ngroups
        ngroups = ngroups + 1

c RANDOM NUMBERS
c phyto either "small" (physize(np)=0.0) or "big" (physize(np)=1.0)
c at this point independent of whether diatom or coccolithophor or not
        RandNo = darwin_random(myThid)
        if(RandNo .gt. 0.500 _d 0)then
          physize(np) = 1.0 _d 0
        else
          physize(np) = 0.0 _d 0
        end if
#ifdef TWO_SPECIES_SETUP
        if (np.eq.1) physize(np) = 1.0 _d 0
        if (np.eq.2) physize(np) = 0.0 _d 0
#endif
#ifdef NINE_SPECIES_SETUP
        if (np.lt.3.or.np.eq.6.or.np.eq.9) then
          physize(np) = 1.0 _d 0
        else
          physize(np) = 0.0 _d 0
        end if
#endif

c size of phytoplankton
        if(physize(np).eq. 1.0 _d 0)then
          dm = 10. _d 0  ! diameter (micrometer)
        else
          dm = 1. _d 0  ! diameter (micrometer)
        end if
c phytoplankton volume in micrometers cubed
        volp=4. _d 0/3. _d 0 *PI*(dm/2. _d 0)**3 _d 0
c
c common block variables (in m and m3)
        phyto_esd(np)=dm* 1. _d -6
        phyto_vol(np)=volp* 1. _d -18
c
c phyto either diatoms (diacoc=1.0) and use silica or cocolithophor
c (diacoc=2.0) and produce PIC or neither (diacoc=0.0)
c if they are large
        if (physize(np).eq.1.0 _d 0) then
          RandNo = darwin_random(myThid)
          if(RandNo .gt. 0.500 _d 0)then
            diacoc(np) = 1.0 _d 0
          else
            diacoc(np) = 0.0 _d 0
          end if
c         if(RandNo .gt. 0.670 _d 0)then
c           diacoc(np) = 1.0 _d 0
c         endif
c         if(RandNo .gt. 0.330 _d 0 .and. RandNo. le. 0.67 _d 0)then
c           diacoc(np) = 2.0 _d 0
c         endif
c         if (RandNo .le. 0.330 _d 0) then
c           diacoc(np) = 0.0 _d 0
c         endif
        else
           diacoc(np) = 0.0 _d 0
        endif
#ifdef TWO_SPECIES_SETUP
        diacoc(np) = 0.0 _d 0
#endif
#ifdef NINE_SPECIES_SETUP
        if (np.eq.1) then
          diacoc(np) = 1.0 _d 0
        else
          diacoc(np) = 0.0 _d 0
        endif
        if (np.eq.9) then
          diacoc(np) = 2.0 _d 0
        endif
#endif
c TEST ...........................................
c       diacoc(np) = 0.0 _d 0
c       write(msgBuf,'(A,I4,A)')
c    &     'generate Phyto: np = ',np,' FIXED - no DIAZO'
c       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
c    &                      SQUEEZE_RIGHT , mythid)
c TEST ...........................................



c phyto either diazotrophs (diazotroph=1.0) or not (diazotroph=0.0)
        RandNo = darwin_random(myThid)
        if(RandNo .gt. 0.6700 _d 0)then
          diazotroph(np) = 1.0 _d 0
        else
          diazotroph(np) = 0.0 _d 0
        end if
c TEST ...........................................
#ifndef ALLOW_DIAZ
        diazotroph(np) = 0.0 _d 0
        write(msgBuf,'(A,I4,A)')
     &     'generate Phyto: np = ',np,' FIXED - no DIAZO'
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT , mythid)
#endif
c TEST ...........................................
#ifdef TWO_SPECIES_SETUP
        diazotroph(np) = 0.0 _d 0
#endif
#ifdef NINE_SPECIES_SETUP
        if (np.gt.5.and.np.lt.8) then
           diazotroph(np) = 1.0 _d 0
        else
           diazotroph(np) = 0.0 _d 0
        end if
#endif


c growth rates
        RandNo = darwin_random(myThid)
c big/small phyto growth rates..
        if(physize(np) .eq. 1.0 _d 0)then
          growthdays = Biggrow +Randno*Biggrowrange
        else
          growthdays = Smallgrow +RandNo*Smallgrowrange
        end if
c but diazotrophs always slower due to energetics
        if(diazotroph(np) .eq. 1.0 _d 0) then
            growthdays = growthdays * diaz_growfac
        endif 
#ifdef TWO_SPECIES_SETUP
        if(physize(np) .eq. 1.0 _d 0)then
          growthdays = Biggrow 
        else
          growthdays = Smallgrow
        end if
#endif
#ifdef NINE_SPECIES_SETUP
        if(physize(np) .eq. 1.0 _d 0)then
          growthdays = Biggrow
        else
          growthdays = Smallgrow
        end if
c but diazotrophs always slower due to energetics
        if(diazotroph(np) .eq. 1.0 _d 0) then
            growthdays = growthdays * diaz_growfac
        endif
c cocco have slower growth than other large
        if (diacoc(np).eq.2. _d 0) then
           growthdays= growthdays * 1.3 _d 0
        endif
c diatom has faster thatn other large
        if (diacoc(np).eq.1. _d 0) then
           growthdays= growthdays * 1. _d 0 ! 0.95 _d 0
        endif
#endif
c now convert to a growth rate
        if (growthdays.gt.0. _d 0) then
         mu(np) = 1.0 _d 0/(growthdays*pday)
        else
         mu(np) = 0. _d 0
        endif

c mortality and export fraction rates
        RandNo = darwin_random(myThid)
c big/small phyto mortality rates..
        if(physize(np) .eq. 1.0 _d 0)then
          mortdays = Bigmort +Randno*Bigmortrange
          ExportFracP(np)=Bigexport
        else
          mortdays = Smallmort +RandNo*Smallmortrange
          ExportFracP(np)=Smallexport
        end if
#ifdef TWO_SPECIES_SETUP
        if(physize(np) .eq. 1.0 _d 0)then
          mortdays = Bigmort
        else
          mortdays = Smallmort
        end if
#endif
#ifdef NINE_SPECIES_SETUP
        if(physize(np) .eq. 1.0 _d 0)then
          mortdays = Bigmort
        else
          mortdays = Smallmort
        end if
#endif

c now convert to a mortality rate
        if (mortdays.gt.0. _d 0) then
          mortphy(np) = 1.0 _d 0/(mortdays*pday)
        else
          mortphy(np) = 0. _d 0
        endif



c nutrient source 
         if(diazotroph(np) .ne. 1.0 _d 0)then
           RandNo = darwin_random(myThid)
           if (physize(np).eq.1.0 _d 0) then   
             nsource(np) = 3
           else
             if(RandNo .gt. 0.670 _d 0)then
               nsource(np) = 1
             elseif(RandNo .lt. 0.33 _d 0)then
               nsource(np) = 2
             else
               nsource(np) = 3
             endif
c ANNA shift bias away from pros. Now equal chance of being HL, LL, Syn, Euk.
c ANNA i.e. now 50% chance of being Pro (nsource 1 or 2, with 50% change of each being HL)
c ANNA i.e. and 50% chance of being non-Pro (nsource 3, with 50% chance of non-pro being Syn)
c           if(RandNo .gt. 0.50 _d 0)then
c             nsource(np) = 3
c           elseif(RandNo .lt. 0.25 _d 0)then
c             nsource(np) = 2
c           else
c             nsource(np) = 1
c           endif 
          endif
        else
              nsource(np) = 0
        end if 
#ifdef TWO_SPECIES_SETUP
        nsource(np) = 3
#endif
#ifdef NINE_SPECIES_SETUP
        if (np.lt.4) then
          nsource(np) = 3
        end if
        nsource(4)=2
        nsource(5)=1
        if (np.gt.5.and.np.lt.8) then
          nsource(np) = 0
        end if
        if (np.gt.7) then
          nsource(np) = 3
        end if
#endif

c.....................................................
c ANNA make selections for WAVEBANDS 
c.....................................................
#ifdef WAVEBANDS
c for now, choice of four absorption spectra types
c pros get either 'HL' or 'LL'
c small others get 'syn' or 'euk'
c large get 'euk'
c each 'type', once assigned, gets given actual values in wavebands_init_vari.F

c ANNA_Q could use tricho abs and scattering spectra (Subramanian et al. 1999)
c ANNA_Q think diaz is turned off for now
c Diaz will be 0 if not defined, and will have nsource = 0. 
       if (tnabp.eq.4) then
         if (nsource(np).eq.0) then   !if diazotroph
          if (physize(np).eq.1.0d0) then !if BIG
          ap_type(np) = 1                !euk (assume diatom association)
          else                           !or
          ap_type(np) = 2                !syn (for now - tricho has billins)
          end if
         end if        

         RandNo = darwin_random(myThid)
         if (nsource(np).eq.3) then !if all three sources (NO3)
          if (physize(np).eq.1.0d0) then !if BIG 
          ap_type(np) = 1                !euk
          else                           !if SMALL
           if (RandNo.gt.0.500d0) then
           ap_type(np) = 1               !euk
           else                          !or
           ap_type(np) = 2               !Syn       
           end if
          end if
         endif
 
         RandNo = darwin_random(myThid)
         if (nsource(np).eq.2) then !if NH4 only
          if (RandNo.gt.0.500d0) then
          ap_type(np) = 3               !Pro HL   
          else                          !or
          ap_type(np) = 4               !Pro LL                 
          end if
         end if

         RandNo = darwin_random(myThid)
         if (nsource(np).eq.1) then !if NH4 & NO2
          if (RandNo.gt.0.500d0) then
          ap_type(np) = 3               !Pro HL   
          else                          !or
          ap_type(np) = 4               !Pro LL                 
          end if
         end if
        endif
c
        if (tnabp.eq.12) then
        if (nsource(np).eq.0) then   !if diazotroph
          if (physize(np).eq.1.0d0) then !if BIG
           if (diacoc(np).eq.1.0d0) then
             ap_type(np) = 5                !diatom association
           endif
           if (diacoc(np).eq.0.0d0) then
             ap_type(np) = 7                !tricho
           endif
           if (diacoc(np).eq.2.0d0) then
             ap_type(np) = 6                !coccolithopher(?)
           endif
          else                           !or
           ap_type(np) = 1                !unicellular (whould be 8 -
                                          !but currently zero)
          end if
         end if

         RandNo = darwin_random(myThid)
         if (nsource(np).eq.3) then !if all three sources (NO3)
          if (physize(np).eq.1.0d0) then !if BIG
           if (diacoc(np).eq.1.0d0) then
             ap_type(np) = 5                !diatom
           endif
           if (diacoc(np).eq.0.0d0) then
             ap_type(np) = 9                !Lg Euk
           endif
           if (diacoc(np).eq.2.0d0) then
             ap_type(np) = 6                !coccolithopher
           endif
          else                           !if SMALL
           if (RandNo.gt.0.500d0) then
           ap_type(np) = 1               !euk
           else                          !or
           ap_type(np) = 2               !Syn
           end if
          end if
         endif
        endif

#ifdef TWO_SPECIES_SETUP
         if (np.eq.1) ap_type(np) = 10
         if (np.eq.2) ap_type(np) = 10
#endif
#ifdef NINE_SPECIES_SETUP
         if (np.eq.1) ap_type(np) = 5
         if (np.eq.2) ap_type(np) = 9
         if (np.eq.3) ap_type(np) = 2
         if (np.eq.4) ap_type(np) = 3
         if (np.eq.5) ap_type(np) = 4
         if (np.eq.6) ap_type(np) = 7
         if (np.eq.7) ap_type(np) = 8
         if (np.eq.8) ap_type(np) = 1
         if (np.eq.9) ap_type(np) = 6
         ap_type(np) = 10
#endif


#else
c ANNA number of RandNo's carreid out MUST MATCH regardless of wavebands or not.
C ANNA the number of RandNo statements here MUST MATCH the number done above

c        RandNo = darwin_random(myThid)
c        RandNo = darwin_random(myThid)
c        RandNo = darwin_random(myThid)

#endif
c ANNA endif


c..........................................................
c generate phyto Temperature Function parameters  
c.......................................................
        phytoTempCoeff(np) = tempcoeff1
        phytoTempExp1(np) = tempcoeff3
        if(physize(np) .eq. 1.0 _d 0)then
          phytoTempExp2(np) = tempcoeff2_big
        else
          phytoTempExp2(np) = tempcoeff2_small
        endif

        RandNo = darwin_random(myThid)
#ifdef TEMP_RANGE
cswd    phytoTempOptimum(np) = 30.0 _d 0 - RandNo*28.0 _d 0 
        phytoTempOptimum(np) = tempmax - RandNo*temprange
        phytoDecayPower(np) = tempdecay
#else
        phytoTempOptimum(np) = 0. _d 0
        phytoDecayPower(np) = 0. _d 0
#endif
        
        write(msgBuf,'(A,I4,A,1P1G24.15E3)')
     &     'generate Phyto: np = ',np,' Topt       =',
     &     phytoTempOptimum(np)
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT , mythid)

c ...............................................
        write(msgBuf,'(A,I4,A,1P1G24.15E3)')
     &     'generate Phyto: np = ',np,' growthdays =', growthdays
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT , mythid)
c ...............................................

c stoichiometric ratios for each functional group of phyto 
c relative to phosphorus  - the base currency nutrient
c set Si:P
        if(diacoc(np) .eq. 1.0 _d 0)then
          R_SiP(np) =  val_R_SiP_diatom
        else
          R_SiP(np) = 0.0 _d 0
        end if
        if(diacoc(np) .eq. 2.0 _d 0)then
          R_PICPOC(np) =  val_R_PICPOC
        else
          R_PICPOC(np) = 0.0 _d 0
        end if
cswd %%%%%%%%%%% OCMIP STYLE for other phyto (not diatom or prochl)
        if (np.eq.2.or.np.eq.3.or.np.eq.8) then
C ECCO-Darwin JAMES paper value
          R_PICPOC(np) = 0.0419 _d 0
cBX GF optizm run ag4         R_PICPOC(np) = 0.133 _d 0
        endif
c set N:P and iron requirement according to diazotroph status
        if(diazotroph(np) .eq. 1.0 _d 0)then
          R_NP(np) = val_R_NP_diaz
          R_FeP(np) =  val_RFeP_diaz
        else
          R_NP(np) = val_R_NP
          R_FeP(np) = val_RFeP
        end if
c set C:P ratio
          R_PC(np) = val_R_PC
c set sinking rates according to allometry
        if(physize(np) .eq. 1.0 _d 0)then
           wsink(np) = BigSink
        else 
           wsink(np) = SmallSink
        end if 
c half-saturation coeffs 

        RandNo = darwin_random(myThid)
        if(physize(np) .eq. 1.0 _d 0)then
           ksatPO4(np) = BigPsat + RandNo*BigPsatrange
        else
c          ksatPO4(np) = SmallPsat + RandNo*SmallPsatrange
c          if (nsource(np).lt.3) then
c            ksatPO4(np) = ksatPO4(np)*prochlPsat
c           endif
           if (nsource(np).eq.3) then
             ksatPO4(np) = SmallPsat + RandNo*SmallPsatrange
           endif
           if (nsource(np).eq..0) then
c            ksatPO4(np) = SmallPsat + RandNo*SmallPsatrange
             ksatPO4(np) = UniDzPsat + RandNo*UniDzPsatrange 
           endif
           if (nsource(np).eq.2.or.nsource(np).eq.1) then
             ksatPO4(np) = ProcPsat + RandNo*ProcPsatrange
           endif
        endif
#ifdef TWO_SPECIES_SETUP
        if(physize(np) .eq. 1.0 _d 0)then
           ksatPO4(np) = BigPsat 
        else
           ksatPO4(np) = SmallPsat
        endif
#endif
#ifdef NINE_SPECIES_SETUP
        if(physize(np) .eq. 1.0 _d 0)then
           ksatPO4(np) = BigPsat
        else
           ksatPO4(np) = SmallPsat
        endif
        if (nsource(np).eq.2.or.nsource(np).eq.1) then
           ksatPO4(np) = ProcPsat 
        endif
        if (diacoc(np) .eq. 2.0 _d 0) then
           ksatPO4(np) = ksatPO4(np)/1.2 _d 0 
        endif
#endif

        ksatNO3(np) = ksatPO4(np)*R_NP(np)
        ksatNO2(np) = ksatNO3(np)*ksatNO2fac 
c Made ksatNH4 smaller since it is the preferred source
        ksatNH4(np) = ksatNO3(np)*ksatNH4fac
        ksatFeT(np) = ksatPO4(np)*R_FeP(np)
        ksatSi(np)  = val_ksatsi

#ifndef GEIDER
cNEW Light parameters:
c     ksatPAR {0.1 - 1.3}
c     0.35=Av High Light Adapted, 0.8=Av Low Light Adapted
c     kinhib  {0.0 - 3.0}
c     0.5 =Av High Light Adapted, 2.0=Av Low Light Adapted
c High Light Groups for Large size:
              if(physize(np) .eq. 1.0 _d 0)then
                 RandNo = darwin_random_normal(myThid)
                 ksatPAR(np) = abs(Bigksatpar+Bigksatparstd*RandNo)

                 RandNo = darwin_random_normal(myThid)
                 kinhib(np) = abs(Bigkinhib+Bigkinhibstd*RandNo)
              else
c QQ remove someday
                 RandNo = darwin_random(myThid)
c Low Light Groups for Small size:
                 RandNo = darwin_random_normal(myThid)
                 ksatPAR(np) = abs(smallksatpar+smallksatparstd*RandNo)

                 RandNo = darwin_random_normal(myThid)
                 kinhib(np) = abs(smallkinhib+smallkinhibstd*RandNo)
              endif
#ifdef TWO_SPECIES_SETUP
              if(physize(np) .eq. 1.0 _d 0)then
                 ksatPAR(np) = abs(Bigksatpar)
                 kinhib(np) = abs(Bigkinhib)
              else
                 ksatPAR(np) = abs(smallksatpar)
                 kinhib(np) = abs(smallkinhib)
              endif
#endif
#ifdef NINE_SPECIES_SETUP
              if(physize(np) .eq. 1.0 _d 0)then
                 ksatPAR(np) = abs(Bigksatpar)
                 kinhib(np) = abs(Bigkinhib)
              else
                 ksatPAR(np) = abs(smallksatpar)
                 kinhib(np) = abs(smallkinhib)
              endif
              if (np.eq.5) then
                kinhib(np) = 6.0 _d 0
              endif
              if (np.eq.9) then
                kinhib(np) = 0.5 _d 0
              endif
#endif
         write(msgBuf,'(A,I4,A,1P1G24.15E3)')
     &     'generate Phyto: np = ',np,' ksatPAR    =', ksatPAR(np)
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT , mythid)
         write(msgBuf,'(A,I4,A,1P1G24.15E3)')
     &     'generate Phyto: np = ',np,' kinhib     =', kinhib(np)
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT , mythid)
#endif

#ifdef GEIDER
        RandNo = darwin_random(myThid)
c big/small phyto growth rates..
        if(physize(np) .eq. 1.0 _d 0)then
          growthdays = Biggrow +Randno*Biggrowrange
        else
          growthdays = Smallgrow +RandNo*Smallgrowrange
        end if
c but diazotrophs always slower due to energetics
        if(diazotroph(np) .eq. 1.0 _d 0) then
            growthdays = growthdays * diaz_growfac
        endif
c cocco have slower growth than diatom
        if (diacoc(np).eq.2. _d 0) then
           growthdays= growthdays * 1.3 _d 0
        endif
#ifdef TWO_SPECIES_SETUP
        if(physize(np) .eq. 1.0 _d 0)then
          growthdays = Biggrow 
        else
          growthdays = Smallgrow 
        end if
#endif
#ifdef NINE_SPECIES_SETUP
        if(physize(np) .eq. 1.0 _d 0)then
          growthdays = Biggrow
        else
          growthdays = Smallgrow
        end if
c but diazotrophs always slower due to energetics
        if(diazotroph(np) .eq. 1.0 _d 0) then
            growthdays = growthdays * diaz_growfac
        endif
c cocco have slower growth than other large
        if (diacoc(np).eq.2. _d 0) then
           growthdays= growthdays * 1.3 _d 0
        endif
c diatom has faster thatn other large
        if (diacoc(np).eq.1. _d 0) then
           growthdays= growthdays * 1.0 _d 0 ! 0.95 _d 0
        endif
#endif
c now convert to a growth rate
        if (growthdays.gt.0. _d 0) then
         pcmax(np) = 1.0 _d 0/(growthdays*pday)
        else
         pcmax(np) = 0. _d 0
        endif
c
c photo-inhibition 
#ifdef WAVEBANDS
c only LL Pro are inhibited
         if (ap_type(np).eq.4) then
            inhibcoef_geid(np) = inhibcoef_geid_val
         else
            inhibcoef_geid(np) = 0. _d 0
         endif
#else
c no inhibition
         if(physize(np) .eq. 1.0 _d 0)then
           inhibcoef_geid(np) = 0. _d 0
         else
           inhibcoef_geid(np) = 0. _d 0  !inhibcoef_geid_val
         endif
#endif
c
        RandNo = darwin_random(myThid)

c big/small phyto PI slope (chl specific)
c        if(physize(np) .eq. 1.0 _d 0)then
c          alphachl(np) = Bigalphachl +Randno*Bigalphachlrange
c        else
c          alphachl(np) = Smallalphachl +RandNo*Smallalphachlrange
c        end if

c ANNA gieder via mQyield instead of alpha
c big/small phyto Maximum Quantum Yield
        if(physize(np) .eq. 1.0 _d 0)then
          mQyield(np) = BigmQyield +Randno*BigmQyieldrange
        else
          mQyield(np) = SmallmQyield +RandNo*SmallmQyieldrange
        end if
#ifdef TWO_SPECIES_SETUP
        if(physize(np) .eq. 1.0 _d 0)then
          mQyield(np) = BigmQyield
        else
          mQyield(np) = SmallmQyield
        end if
#endif
#ifdef NINE_SPECIES_SETUP
        if(physize(np) .eq. 1.0 _d 0)then
          mQyield(np) = BigmQyield
        else
          mQyield(np) = SmallmQyield
        end if
#endif
#ifdef WAVEBANDS
c ANNA for wavebands only, re-set mQyield to be constant for all np's
c ANNA i.e. let alpha vary only with aphy_chl_ps
c ANNA value is mean of vals for big and small.
          mQyield(np) = 4.0 _d -5
#endif

        RandNo = darwin_random(myThid)
c big/small phyto C:Chl max
        if(physize(np) .eq. 1.0 _d 0)then
          chl2cmax(np) = Bigchl2cmax +Randno*Bigchl2cmaxrange
        else
          chl2cmax(np) = Smallchl2cmax +RandNo*Smallchl2cmaxrange
        end if
#ifdef TWO_SPECIES_SETUP
        if(physize(np) .eq. 1.0 _d 0)then
          chl2cmax(np) = Bigchl2cmax 
        else
          chl2cmax(np) = Smallchl2cmax 
        end if
#endif
#ifdef NINE_SPECIES_SETUP
        if(physize(np) .eq. 1.0 _d 0)then
          chl2cmax(np) = Bigchl2cmax
        else
          chl2cmax(np) = Smallchl2cmax
        end if
#endif
c ANNA chl2cmin added
c       chl2cmin(np) = 0.003  _d 0 * 12. _d 0  ! mg Chl a/mmol C

         write(msgBuf,'(A,I4,A,1P1G24.15E3)')
     &     'generate Phyto: np = ',np,' pcmax    =', pcmax(np)
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT , mythid)
c         write(msgBuf,'(A,I4,A,1P1G24.15E3)')
c     &     'generate Phyto: np = ',np,' alphachl    =', alphachl(np)
c ANNA CHANGED TO MQYIELD from ALPHACHL
c ANNA STEPH msgBuf changed for mQyield?
         write(msgBuf,'(A,I4,A,1P1G24.15E3)')
     &     'generate Phyto: np = ',np,' mQyield    =', mQyield(np)
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT , mythid)
         write(msgBuf,'(A,I4,A,1P1G24.15E3)')
     &     'generate Phyto: np = ',np,' chl2cmax     =', chl2cmax(np)
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT , mythid)
#endif

#ifdef DAR_DIAG_CHL
        if(physize(np) .eq. 1.0 _d 0)then
          Geider_alphachl(np) = Geider_Bigalphachl
          Geider_chl2cmax(np) = Geider_Bigchl2cmax
          Geider_chl2cmin(np) = Geider_Bigchl2cmin
        else
          Geider_alphachl(np) = Geider_smallalphachl
          Geider_chl2cmax(np) = Geider_smallchl2cmax
          Geider_chl2cmin(np) = Geider_smallchl2cmin
        end if

         write(msgBuf,'(A,I4,A,1P1G24.15E3)')
     &     'generate Phyto: np = ',np,' Geider_alphachl =',
     &     Geider_alphachl(np)
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT , mythid)
         write(msgBuf,'(A,I4,A,1P1G24.15E3)')
     &     'generate Phyto: np = ',np,' Geider_chl2cmax =',
     &     Geider_chl2cmax(np)
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT , mythid)
#endif

        RETURN
        END
#endif  /*DARWIN*/
#endif  /*ALLOW_PTRACERS*/

c ===========================================================