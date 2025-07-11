C $Header: /home/ubuntu/mnt/e9_copy/MITgcm_contrib/ecco_darwin/v4_llc270/code_darwin/DARWIN_SIZE.h,v 1.2 2019/08/26 05:33:09 dcarroll Exp $
C $Name:  $

c DARWIN_SIZE.h 
c description: Size specification for generalized ecosystem model
c      Mick Follows, Scott Grant Fall/Winter 2005
c      Stephanie Dutkiewicz Spring/Summer 2006
c
c
c npmax = no of "functional groups" of phytoplankton
c nzmax = no of "functional groups" of zooplankton
c
c
         INTEGER npmax
         INTEGER nzmax
         PARAMETER(npmax=5,nzmax=2)

C iPO4          :: index of PO4  in Ptracer
C iNO3          :: index of NO3  in Ptracer
C iFeT          :: index of FeT  in Ptracer
C iSi           :: index of Si   in Ptracer
C iDOP          :: index of DOP  in Ptracer
C iDON          :: index of DON  in Ptracer
C iDOFe         :: index of DOFe in Ptracer
C iPOP          :: index of POP  in Ptracer
C iPON          :: index of PON  in Ptracer
C iPOFe         :: index of POFe in Ptracer
C iPOSi         :: index of POSi in Ptracer
C iNH4          :: index of NH4  in Ptracer
C iNO2          :: index of NO2  in Ptracer
C iCa           :: index of Ca   in Ptracer
C iZoo          :: index of first zooplankton
C iPhy          :: index of first phytoplankton
C iChl          :: index of first phytoplankton Chl (if using dynamic chl)
C nDarwin       :: total number of ptracers used by DARWIN
C nCompZooMax   :: maximum number of components each zooplankton can have (P,N,...)
C strideCompZoo :: increment between components of zooplankton
C               ::   iZooComp(nz) = iZoo+(nz-1)*strideTypeZoo+(iComp-1)*strideCompZoo
C strideTypeZoo :: increment between types of zooplankton
C               ::   iZooComp(nz) = iZoo+(nz-1)*strideTypeZoo+(iComp-1)*strideCompZoo
C the only reasonable choices for the above two are:
C   strideCompZoo = 1
C   strideTypeZoo = nCompZooMax
C (the traditional choice) or
C   strideCompZoo = nzmax
C   strideTypeZoo = 1
C remember to bring the fields in data.ptracers in the right order !

      INTEGER iPO4
      INTEGER iNO3
      INTEGER iFeT
      INTEGER iSi
      INTEGER iDOP
      INTEGER iDON
      INTEGER iDOFe
      INTEGER iPOP
      INTEGER iPON
      INTEGER iPOFe
      INTEGER iPOSi
      INTEGER iNH4
      INTEGER iNO2
      INTEGER iZoo
      INTEGER iPhy
      INTEGER nDarwin
      INTEGER nCompZooMax
      INTEGER strideCompZoo
      INTEGER strideTypeZoo
      INTEGER iTot
#ifdef ALLOW_CARBON
      INTEGER iDIC
      INTEGER iDOC
      INTEGER iPOC
      INTEGER iPIC
      INTEGER iALK
      INTEGER iO2
      INTEGER iZoC
      INTEGER iCA
#endif
#ifdef DYNAMIC_CHL
      INTEGER iChl
#endif
      PARAMETER (nCompZooMax=4)
      PARAMETER (strideCompZoo=1)
      PARAMETER (strideTypeZoo=nCompZooMax)
      PARAMETER (iPO4  =1)
      PARAMETER (iNO3  =2)
      PARAMETER (iFeT  =3)
      PARAMETER (iSi   =4)
      PARAMETER (iDOP  =5)
      PARAMETER (iDON  =6)
      PARAMETER (iDOFe =7)
      PARAMETER (iZoo  =8)
      PARAMETER (iPOP  =iZoo+nzmax*nCompZooMax)
      PARAMETER (iPON  =iPOP  +1)
      PARAMETER (iPOFe =iPON  +1)
      PARAMETER (iPOSi =iPOFe +1)
      PARAMETER (iNH4  =iPOSi +1)
      PARAMETER (iNO2  =iNH4  +1)
      PARAMETER (iPhy  =iNO2  +1)
#ifdef DYNAMIC_CHL
      PARAMETER (iChl =iPhy +npmax)
      PARAMETER (iTOT =iPhy +npmax+npmax)
#else
      PARAMETER (iTOT =iPhy +npmax)
#endif
#ifdef ALLOW_CARBON
      PARAMETER (iDIC  =iTOT)
      PARAMETER (iDOC  =iDIC+1)
      PARAMETER (iPOC  =iDOC+1)
      PARAMETER (iPIC  =iPOC+1)
      PARAMETER (iALK  =iPIC+1)
      PARAMETER (iO2   =iALK+1)
      PARAMETER (iZoC  =iO2+1)
      PARAMETER (iCA = iZoC+nzmax)
      PARAMETER (nDarwin=iCA-1)
#else
      PARAMETER (nDarwin=iTOT-1)
#endif
C iZooP(nzmax)  :: index of phosphorus content of each zooplankton type
C iZooN(nzmax)  :: index of nitrogen   content of each zooplankton type
C iZooFe(nzmax) :: index of iron       content of each zooplankton type
C iZooSi(nzmax) :: index of silica     content of each zooplankton type
C these are computed in darwin_init_fixed from iZoo, strideCompZoo and strideTypeZoo
      INTEGER iZooP (nzmax)
      INTEGER iZooN (nzmax)
      INTEGER iZooFe(nzmax)
      INTEGER iZooSi(nzmax)
#ifdef ALLOW_CARBON
      INTEGER iZooC(nzmax)
#endif
      COMMON/DARWIN_SIZE/ iZooP, iZooN, iZooFe, iZooSi
#ifdef ALLOW_CARBON
     &                  ,iZooC
#endif

