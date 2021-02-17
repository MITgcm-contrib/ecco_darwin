#ifdef ALLOW_CARBON
c air-sea flux info
C permil : is conversion factor for mol/m3 to mol/kg
C          assumes uniform (surface) density
C Pa2Atm : for conversion of atmospheric pressure
C          when coming from atmospheric model
       COMMON /CARBON_NEEDS/
     &              AtmospCO2, AtmosP, pH, pCO2, FluxCO2,
     &              wind, Kwexch_Pre,
     &                FluxO2,FluxO2TH,FluxO2AB,
     &              FluxN2O,FluxN2OTH,FluxAR
      _RL  AtmospCO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  AtmosP(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  pH(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  pCO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  FluxCO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  wind(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  Kwexch_Pre(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

CC ADDING NEW FLUXES

      _RL  FluxO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  FluxO2TH(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  FluxO2AB(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  FluxN2O(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  FluxN2OTH(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy) 
      _RL  FluxAR(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)



C permil : is conversion factor for mol/m3 to mol/kg
C          assumes uniform (surface) density
C Pa2Atm : for conversion of atmospheric pressure
C          when coming from atmospheric model
       COMMON /GLOBAL_SURF_MEAN/
     &                          gsm_alk,gsm_s,gsm_t,gsm_dic,
     &                          gsm_c14,permil,Pa2Atm
      _RL  gsm_alk
      _RL  gsm_s
      _RL  gsm_t
      _RL  gsm_DIC
      _RL  gsm_C14
      _RL  permil
      _RL  Pa2Atm

       COMMON /CARBON_CHEM/
     &                     ak0,ak1,ak2,akw,akb,aks,akf,
     &                     ak1p,ak2p,ak3p,aksi, fugf,
     &                     ff,ft,st,bt, Ksp_TP_Calc
      _RL  ak0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  ak1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  ak2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  akw(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  akb(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  aks(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  akf(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  ak1p(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  ak2p(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  ak3p(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  aksi(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C Fugacity Factor added following Val Bennington 
      _RL  fugf(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  ff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  ft(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  st(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  bt(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  Ksp_TP_Calc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

       COMMON /OXYGEN_CHEM/
     &              oA0,oA1,oA2,oA3,oA4,oA5,
     &              oB0,oB1,oB2,oB3,
     &              oC0
      _RL oA0,oA1,oA2,oA3,oA4,oA5
      _RL oB0,oB1,oB2,oB3
      _RL oC0

      COMMON /ARGON_CHEM/

     &              A0Ar,A1Ar,A2Ar,A3Ar,A4Ar,A5Ar,
     &              B0Ar,B1Ar,B2Ar,B3Ar,
     &              C0Ar
      _RL A0Ar,A1Ar,A2Ar,A3Ar,A4Ar,A5Ar
      _RL B0Ar,B1Ar,B2Ar,B3Ar
      _RL C0Ar





C schmidt number coefficients
      COMMON /DIC_SCHMIDT_NO/
     &                    sca1, sca2, sca3, sca4,
     &                    sox1, sox2, sox3, sox4,
C MM ADDING THE SAME FOR N2O and ARGON     
C - MM - OCTOBER 2015
     &                   sn2o1,sn2o2,sn2o3,sn2o4,
     &                   sar1,sar2,sar3,sar4      
      _RL  sca1
      _RL  sca2
      _RL  sca3
      _RL  sca4
      _RL  sox1
      _RL  sox2
      _RL  sox3
      _RL  sox4
CMM  NEW COEFFS

      _RL sn2o1
      _RL sn2o2
      _RL sn2o3
      _RL sn2o4

      _RL sar1
      _RL sar2
      _RL sar3
      _RL sar4




      COMMON /DIC_DIAGS_NEEDS/
     &     SURave, SURCave, SUROave, pCO2ave, pHave,     
     &     dic_timeave, SUROTHave,SUROABave,
     &     SURN2Oave,SURN2OTHave,SURARave,KwO2ave,
     &     dO2satave,NCPave
      _RL SURave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL SURCave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL SUROave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL pCO2ave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL pHave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DIC_timeave(nSx,nSy,nR)
C MM
C ADDING NEW TRACERS
       _RL SUROTHave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
       _RL SUROABave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
       _RL SURN2Oave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
       _RL SURN2OTHave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy) 
       _RL SURARave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
       _RL KwO2ave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
       _RL dO2satave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy) 
       _RL NCPave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,1:Nr,nSx,nSy)

   






#endif
