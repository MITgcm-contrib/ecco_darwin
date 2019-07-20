#ifdef ALLOW_CARBON
c air-sea flux info
C permil : is conversion factor for mol/m3 to mol/kg
C          assumes uniform (surface) density
C Pa2Atm : for conversion of atmospheric pressure
C          when coming from atmospheric model
       COMMON /CARBON_NEEDS/
     &              AtmospCO2, AtmosP, pH, pCO2, FluxCO2,
     &              wind, Kwexch_Pre
      _RL  AtmospCO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  AtmosP(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  pH(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  pCO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  FluxCO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  wind(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  Kwexch_Pre(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

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
     &                     ff,ft,st,bt, Ksp_TP_Calc,CO3,
     &                     omegaC
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
      _RL  CO3(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  omegaC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx,nSy)

       COMMON /CO2_FLUX_BUDGET/
     &  budgetPvPH,budgetPvPCO2,budgetPvCO3,
     &  budgetTempPH,budgetTempPCO2,budgetTempCO3,
     &  budgetSaltPH,budgetSaltPCO2,budgetSaltCO3,
     &  budgetAlkPH,budgetAlkPCO2,budgetAlkCO3,
     &  budgetDicPH,budgetDicPCO2,budgetDicCO3,
     &  fluxCO2_1,Kwexch1,budgetTemp1,budgetSalt1,budgetAlk1,
     &  budgetDic1,budgetPhos1,budgetSi1,dFluxCO2,dFluxCO2Pv,
     &  dFluxCO2Temp,dFluxCO2Salt,dFluxCO2Alk,dFluxCO2Dic,
     &  dFluxCO2Sf,budgetTStep1

      _RL budgetPvPH(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL budgetPvPCO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL budgetPvCO3(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL budgetTempPH(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL budgetTempPCO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL budgetTempCO3(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL budgetSaltPH(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL budgetSaltPCO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL budgetSaltCO3(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL budgetAlkPH(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL budgetAlkPCO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL budgetAlkCO3(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL budgetDicPH(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL budgetDicPCO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL budgetDicCO3(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      
      _RL fluxCO2_1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy) 
      _RL Kwexch1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL budgetTemp1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL budgetSalt1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL budgetAlk1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL budgetDic1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL budgetPhos1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL budgetSi1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
     
      _RL dFluxCO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL dFluxCO2Pv(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL dFluxCO2Temp(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL dFluxCO2Salt(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL dFluxCO2Alk(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL dFluxCO2Dic(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL dFluxCO2Sf(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL budgetTStep1

       COMMON /OXYGEN_CHEM/
     &              oA0,oA1,oA2,oA3,oA4,oA5,
     &              oB0,oB1,oB2,oB3,
     &              oC0
      _RL oA0,oA1,oA2,oA3,oA4,oA5
      _RL oB0,oB1,oB2,oB3
      _RL oC0

C schmidt number coefficients
      COMMON /DIC_SCHMIDT_NO/
     &                    sca1, sca2, sca3, sca4,
     &                    sox1, sox2, sox3, sox4
      _RL  sca1
      _RL  sca2
      _RL  sca3
      _RL  sca4
      _RL  sox1
      _RL  sox2
      _RL  sox3
      _RL  sox4

      COMMON /DIC_DIAGS_NEEDS/
     &     SURave, SURCave, SUROave, pCO2ave, pHave,     
     &     dic_timeave
      _RL SURave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL SURCave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL SUROave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL pCO2ave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL pHave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DIC_timeave(nSx,nSy)

#endif
