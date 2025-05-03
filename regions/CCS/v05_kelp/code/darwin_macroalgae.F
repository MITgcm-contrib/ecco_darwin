#include "DARWIN_OPTIONS.h"

CBOP
C     ROUTINE: DARWIN_MACROALGAE
C     INTERFACE: ==========================================================
      SUBROUTINE DARWIN_MACROALGAE(
     I     t,hor_vel,parl,ptr,
     I     Hw,
     I     Tw,deltaT,
     I     drF,hFacC,
     U     gtr,
     O     gT,gQ,gH,gE,muMAG,
     O     mortMAG,UptakeNout,
     O     mortout,dEout)
C   
C    ==========================================================

C     t:         temperature  [C]
C     hor_vel:   horizontal velocity [m s-1]
C     parl: photosynthetic radiation [mu Ein s-1 m-2]
C     ptr(ndarwin): tracer concentrations  
C     gtr(ndarwin): tracer tendencies (inout)
C     macroalgage relevant tracers in ptr
C     ptr(B) .... Kelp biomass [gDryWeight/m3]
C     ptr(C) .... Kelp biomass in carbon [gC/m3]
C     ptr(Q) .... Kelp nitrogen quota [mgN/gDryWeight]
C     gtr = tendency dptr/dt due to physical processes [ptr s-1]
C     DW: dry weight
C     Based on Broch et al and
C     Strong-Wright paper, https://github.com/jagoosw/SugarKelp.jl

      IMPLICIT NONE
      
#include "DARWIN_PARAMS.h"
#include "DARWIN_SIZE.h"
#include "DARWIN_INDICES.h"
#include "DARWIN_DIAGS.h"

C     Inputs/outputs
      _RL ptr(nDarwin)            ! Darwin tracers
      _RL gTr(nDarwin)            ! Darwin tendencies
      _RL t,hor_vel,parl          ! Darwin water temperature, horizontal velocity and PAR
      _RL deltaT                  ! Darwin time step duration
      _RL drF,hFacC
C     Additional variables needed for physiological calculation in MACMODS
      _RL vQ                      ! Quota-limited uptake
C     Variables related to ocean flow and wave motion affecting nutrient uptake
      _RL wp_tmp                  ! Local water temperature
      _RL Tw                      ! Local wave period
      _RL Tw_s                    ! Maximum wave period
      _RL Hw                      ! Local significant wave height
      _RL n_length                ! number of vertical levels of the kelp forest  
      _RL visc                    ! kinematic viscosity 
      _RL magu                    ! Local horizontal flow
      _RL magu_m                  ! Local maximum horizontal flow
      _RL DBL                     ! Diffusive boundary layer thickness
      _RL lambda_NO3              ! NO3 uptake related term
      _RL vval                    ! Iterations for oscillatory flow
      _RL n                       ! Iterative number for oscillatory flow over vertical levels of the kelp forest
      _RL vNuTw_NO3               ! NO3 uptake factor
      _RL Uptake_NO3              ! Uptake of NO3
      _RL UptakeNout              ! Uptake of N for Darwin tendency
      _RL dEout                   ! Exudation for Darwin tendency
      _RL mortout                 ! Mortality for Darwin tendency
      _RL growthout               ! Growth for Darwin tendency
      _RL ratioCN                 ! Ratio C:N
      _RL UptakeN                 ! Uptake of nitrogen
      _RL Dm                      ! Molecular diffusivity of nitrate at 18degC
      _RL oscillatory             ! Oscillatory flow
      _RL flow                    ! Uni-directional flow
      _RL beta                    ! Oscillatory and Uni-directional flow
      _RL Bmag                    ! local biomass of Kelp 
      _RL Qmag                    ! local nitrogen cell quota
      _RL Cmag                    ! local biomass of Kelp in C
      _RL lambda                  ! Nutrient uptake related term
C     Additional MACMODS variables
      _RL dE                      ! Exudation rate
C     Tendency terms
      _RL dBmagdt                 ! B tendency
      _RL dCmagdt                 ! C tendency
      _RL dQmagdt                 ! Q tendency
      _RL dQgrowthdt              ! Q growth related term     
C     Key variables to store forcing functions from ECCO-Darwin
C     mag = macroalgae
      _RL t_mag                   ! water temperature in degC
      _RL NO3_mag                 ! NO3 conc. in mmol/m3
      _RL chl_plk                 ! chlrophyll from plankton   
      _RL gQ,gT,gE,gH             ! Growth limitation terms
      _RL muMAG,mortMAG,WP        ! Growth and mortality rates, wave mortality
      _RL Mwave, mwh              ! Wave mortality parameters
      _RL A                       ! Carrying capacity related term
      _RL dNs                     ! Nitrogen uptake related term
      _RL N_new                   ! Nitrogen uptake compared to vertical flux of nitrogen
            
C     Additional parameters
      _RL mp_growth_lim_type ! Growth limitation Type
      _RL break_Duarte_Ferr  ! Duarte frond breakage Type
      _RL break_Rodrig       ! Rodrigues frond breakage Type
      _RL mp_break_type      ! Frond breakage Type
      _RL mp_N_uptake_type   ! Type of nutrient uptake
      _RL mp_N_flux_limit    ! Limit N03 uptake by vertica flux of NO3
      _RL mp_spp_kcap_m3     ! Carrying capcity converted to m3
      _RL mp_spp_Vmax_m3     ! Maximum uptake rate of nitrate converted to m3
      _RL c0                 ! 0
      _RL c1                 ! 1
      _RL c4                 ! 4
      _RL pi                 ! PI
      _RL mw_n, mw_c         ! Molecular weight N and C
      _RL spd                ! Number of seconds in a day
      _RL R_ON               ! O:N ratio for respiration and consumption
C     Physical constants used in MACMODS 
      _RL rho_sw             ! density of salt water (g/cm^3)
C     Additional MACMODS parameters with assigned values 
      parameter(break_Duarte_Ferr=0)		
      parameter(break_Rodrig=1)
      parameter(c0=0.0)
      parameter(c1=1.0)
      parameter(c4=4.0)
      parameter(pi=3.14159265358979323846)
      parameter(mw_n=14.00672)
      parameter(mw_c=12.01070)
      parameter(spd=86400.0)
      parameter(R_ON= 16.0/170.0)
C     Additional physical constants for MACMODS
      parameter(rho_sw=1.026)             ! density of salt water (g/cm^3)
C     Assign value for frond breakage type      
      mp_break_type = break_Duarte_Ferr
C     Assign forcings
      NO3_mag=max(c0, ptr(iNO3)) ! mmol/m3
      t_mag=t                    ! deg C   
C     growth Limitation model
C     if = 0 ==>  multiply all the limitaton factors
C     if > 0 ==>  independent light and nutrient limitation growth model
      mp_growth_lim_type=0
C     Setting the kind of nutrient uptake
      mp_N_uptake_type=1
C     Assign state variables to local variables     
      Bmag = ptr(iMagB)          ! [gDW/m3] 
      Qmag = ptr(iMagQ)          ! [mg-N/gDW]
      Cmag = ptr(iMagC)          ! [gC/m3] 
C     ---------------------------------------------------------- 
C     MACROALGAE MODEL         
C     ---------------------------------------------------------- 

C     Parameters of macroalgae in data.darwin already converted to Darwin units
C     following these conversion factors:
C     Conversion parameters from https://doi.org/10.3389/fmars.2022.752951 and
C     https://doi.org/10.1038/s43247-023-00833-2 to Darwin units      
C     mp_spp_Gmax_cap = mp_spp_Gmax_cap/spd          ! from d-1 to s-1
C     mp_spp_Vmax = mp_spp_Vmax * mp_spp_dry_sa      ! umol N/g-DW/d * g-DW/m2 = umol N/m2/d
C     mp_spp_Vmax = mp_spp_Vmax/spd                  ! from d-1 to s-1 [umol N/m2/s]
C     mp_spp_Vmax = mp_spp_Vmax/1E3                  ! from [umol N/m2/s] to [mmol N/m2/s] (Darwin unit)
C     mp_spp_E = mp_spp_E/spd                        ! from d-1 to s-1
C     mp_spp_death = mp_spp_death/spd                ! from d-1 to s-1
C     mp_spp_kcap_rate = mp_spp_kcap_rate/spd        ! from d-1 to s-1
C     mp_spp_Ks_NO3 = mp_spp_Ks_NO3/1E3              ! from [umol N/m3] to [mmol N/m3] (Darwin unit)
C     mp_spp_PARs = mp_spp_PARs*4.57                 ! from [W/m2] to [microEinstein/m2/s]
C     mp_spp_PARc = mp_spp_PARc*4.57                 ! from [W/m2] to [microEinstein/m2/s]
      
C     ----------------------------------------------------------       
C     Lambda NO3 function         
C     ----------------------------------------------------------   
C     Unit conversions, minima
      magu = hor_vel		                                ! [m/s]     from ECCO-Darwin CCS model
      magu_m = max(1.0,magu)                                    ! [m/s]
      Tw_s = max(0.01,Tw)                                       ! [s]       mean wave period
C     some parameters and variables
      n_length = 25 
      visc = 1E-6                                               ! [m2/s]    kinematic viscosity 
      Dm = (18*3.65*1E-11 + 9.72*1E-10)                         ! [m2/s]    molecular diffusivity of nitrate at 18degC 

      DBL = 10 * (visc / (sqrt(mp_spp_CD) * abs(magu_m)))       ! [m]       thickness of diffusive boundary layer 
      vval = 0.0
C     1. Oscillatory flow
      DO n = 1,n_length
        vval = vval + (1-exp((-Dm * n**2 * pi**2 *Tw_s)/(2.0*DBL**2))) 
     & /(n**2 * pi**2)                                          ! [no unit]
      ENDDO
      oscillatory = ((4.0*DBL)/Tw_s) * vval                     ! [m/s]
C     2. Uni-directional Flow
      flow = Dm / DBL                                           ! [m/s]
      beta = flow + oscillatory                                 ! [m/s]
      
C     ----------------------------------------------------------  
C     Nutrient Uptake
C     ----------------------------------------------------------              
      lambda_NO3 = c1 + (mp_spp_Vmax / (beta*mp_spp_Ks_NO3)) - 
     & (NO3_mag/mp_spp_Ks_NO3)
      lambda=lambda_NO3                                          ! [no unit]
            
C     Quota-limited uptake: maximum uptake when Q is minimum and
C     approaches zero as Q increases towards maximum; Possible that Q
C     is greater than Qmax. Set any negative values to zero.      
      IF (mp_N_uptake_type .EQ.1) THEN
       vQ = (mp_spp_Qmax-Qmag)/(mp_spp_Qmax-mp_spp_Qmin)              ! [0-1]
       vQ = max(c0,vQ)                                                ! [0-1]
       vQ = min(c1,vQ)                                                ! [0-1]
      ELSEIF (mp_N_uptake_type .EQ.0) THEN
       vQ = 1.0 - 1.0/(1.0+(max(0.0,mp_spp_Qmax-Qmag))*               ! [0-1]
     & 50.0/(mp_spp_Qmax-mp_spp_Qmin))
       vQ = max(c0,vQ)                                                ! [0-1]
       vQ = min(c1,vQ)                                                ! [0-1]
      ENDIF
      
C     Below is what we call "Uptake Factor." It varies betwen 0
C     and 1 and includes kinetically limited uptake and
C     mass-transfer-limited uptake (oscillatory + uni-directional flow)
      vNuTw_NO3 = NO3_mag / (mp_spp_Ks_NO3 * ((NO3_mag/mp_spp_Ks_NO3)+ 
     & 0.5 *(lambda+sqrt(lambda**2 + c4 * (NO3_mag/mp_spp_Ks_NO3))))) ! [0-1]   
      vNuTw_NO3 = max(c0,vNuTw_NO3)                                   ! [0-1]
      vNuTw_NO3 = min(c1,vNuTw_NO3)                                   ! [0-1]

C     Uptake Rate [mg N/gDW/s]
C     Nutrient Uptake Rate = Max Uptake * v[Ci,u,Tw] * vQ
      mp_spp_Vmax_m3 = mp_spp_Vmax / (drF*hFacC)                      ! [mmolN/m3/s]
      Uptake_NO3 = mp_spp_Vmax_m3 * vNuTw_NO3 * vQ                    ! [mmolN/m3/s]
      UptakeNout = Uptake_NO3                                         ! [mmolN/m3/s]
      Uptake_NO3 = Uptake_NO3 * mw_n                                  ! [mg/m3/s]          
      UptakeN = max(c0,Uptake_NO3 / Bmag)                             ! [mg/gDW/s]     
      
C     ---------------------------------------------------------- 
C     Growth
C     ----------------------------------------------------------
C     Growth, nitrogen movement from Q to B
C     Output:
C     Growth, [s-1]
C     gQ, quota-limited growth
C     from Wheeler and North 1980 Fig. 2
C     gT, temperature-limited growth
C     piecewise approach taken from Broch and Slagstad 2012 (for sugar
C     kelp) and optimized for Macrocystis pyrifera
C     gE, light-limited growth
C     from Dean and Jacobsen 1984
C     gH, carrying capacity-limited growth

C     Nutrient (quota-based) limitation
C     gQ = (Q(i,j) - params(mp_spp_Qmin)) / Q(i,j)                  ! Droop equation      
      gQ = ((Qmag - mp_spp_Qmin) / Qmag) * 
     & (mp_spp_Qmax/(mp_spp_Qmax-mp_spp_Qmin))                      ! Droop scaled from 0-1
C     gQ = (Q(i,j) - Qmin) / (Qmax-Qmin)                            ! Freider et al.
      gQ = max(c0,gQ)                                               ! [0-1]
      gQ = min(c1,gQ)                                               ! [0-1]
      
C     Temperature-dependent growth          
      IF (t_mag .GE. mp_spp_Topt1) THEN          
         IF (t_mag .LE. mp_spp_Topt2) THEN
            gT = 1.0                                           ! [0-1]
         ELSE  
            gT = exp(-mp_spp_K2*(t_mag-mp_spp_Topt2)**2)       ! [0-1]
         ENDIF          
      ELSE              
         gT = exp(-mp_spp_K1*(t_mag-mp_spp_Topt1)**2)          ! [0-1]      
      ENDIF 
      gT = max(0.0,gT)                                         ! [0-1]
      gT = min(1.0,gT)                                         ! [0-1]
    
C     Light limitation
      IF (parl .LT. mp_spp_PARc) THEN              
          gE = c0                                              ! [0-1]
      ELSEIF (parl .GT. mp_spp_PARs) THEN              
          gE = c1                                              ! [0-1]
      ELSE              
          gE = ((parl-mp_spp_PARc)/(mp_spp_PARs-mp_spp_PARc))*
     & exp(-((parl-mp_spp_PARc)/(mp_spp_PARs-mp_spp_PARc))+c1) ! [0-1]  
      ENDIF

C     Carrying capacity
C     gH -> density-limited growth (ranges from the max growth rate to the death rate, excluding wave mortality)
C     This expression follows Xiao et al (2019 and ignores wave mortality when
C     thinking about the death rate
C     Convert carrying capacity from gDW/m2 to gDW/m3 according to grid cell thickness   
      mp_spp_kcap_m3 = mp_spp_kcap/(drF*hFacC) 
      A = mp_spp_kcap_rate/mp_spp_Gmax_cap                          ! [no unit]
      gH = min(c1,A*(Bmag/mp_spp_kcap_m3)**(-0.75))                 ! [0-1]
      
C     Growth rate
C     below is confusing programatically - gH should not include Gmax - but this is what was used in gMACMODS < 2022-10
C     if gH > params[mp_spp_Gmax_cap]:
C         gH = params[mp_spp_Gmax_cap]
C     Growth = gH * gT * gE * gH
C     interacting limitation terms growth model
      IF (mp_growth_lim_type .EQ. 0) THEN
    	 muMAG =  mp_spp_Gmax_cap * gH * gT * gE * gQ                ! [per s]
      ELSE 
C     independent light and nutrient limitation growth model
    	 muMAG =  mp_spp_Gmax_cap * gH * gT * min(gE,gQ)             ! [per s]             
      ENDIF

C     ---------------------------------------------------------- 
C     Compute mortality rate due to wave activity      
C      ----------------------------------------------------------
C     Wave height is computed by the ERA wave model
C     d_wave = frond loss due to waves; dependent on Hs, significant
C     wave height [m]; Rodrigues et al. 2013 demonstrates linear relationship
C     between Hs and frond loss rate in Macrocystis [d-1] (continuous)
C     Duarte and Ferreira (1993) find a linear relationship between wave
C     power and mortality in red seaweed
C     Duarte and Ferreira compute daily wave mortality [per d] as a function of wave power 
      mwh  = Hw                                                            ! Assign wave height [m]
      wp_tmp = Tw                                                          ! Assign wave period [s]
      IF (mp_break_type.EQ.break_Duarte_Ferr ) THEN
         WP = (((rho_sw * (9.8**2)) /(64.0*pi)) * (mwh**2) * wp_tmp) / 1E3 ! [kW]   
         Mwave=((2.3*1E-4) * WP + (2.2*1E-3)) * mp_wave_mort_factor        ! [per d]
      ELSEIF (mp_break_type.EQ.break_Rodrig) THEN
C     Death rate scales with Hw (Hw, surface wave height)
         Mwave=(mp_spp_death*spd * Hw) * mp_wave_mort_factor               ! [per d]
      ENDIF
C     M_wave + general Mortality rate; [per s]       
      mortMAG = mp_spp_death + Mwave/spd                                   ! [per s]           
      
C     Compute exudation
      dE = mp_spp_E * (Qmag-mp_spp_Qmin) * deltaT                          ! [mg-N/gDW]
     
C     ---------------------------------------------------------- 
C     Output terms
C     ----------------------------------------------------------      
C     dilution total N d_Q
      dQgrowthdt = Qmag * (1.0/(1.0 + muMAG * deltaT) - 1.0)   ! [mg-N/gDW]
      dQmagdt = (UptakeN * deltaT + dQgrowthdt - dE)           ! [mg-N/gDW]
      dBmagdt = Bmag * muMAG * deltaT - 
     & Bmag * mortMAG * deltaT                                 ! [gDW/m3] 
      ratioCN = (Cmag*1E3)/(Bmag*Qmag)                         ! [(mgC/m3)/(mgN/m3)]
      IF (Bmag.GT.c0) THEN
        gtr(iMagB) = dBmagdt/deltaT                            ! [gDW/m3/s]
        gtr(iMagC) = (dBmagdt/deltaT) * mp_spp_carbon          ! [gC/m3/s]
        gtr(iMagQ) = dQmagdt/deltaT                            ! [mg-N/gDW/s]
C       Remove NO3 used by macroalgae to NO3 tracer        
        gtr(iNO3) = gtr(iNO3) - UptakeNout                       ! [mmol-N/m3/s] 
C       Add macroalgae exudation to DOC and DON      
        dEout = (dE/deltaT) * Bmag / mw_n                        ! [mmol-N/m3/s]
        gtr(iDOC) = gtr(iDOC) + (dEout * ratioCN)                ! [mmol-C/m3/s]
        gtr(iDON) = gtr(iDON) + dEout                            ! [mmol-N/m3/s]
C       Add macroalgae mortality to POC and PON 
        mortout = Cmag * mortMAG * 1E3 / mw_c                    ! [mmolC/m3/s] 
        gtr(iPOC) = gtr(iPOC) + mortout                          ! [mmolC/m3/s]   
        gtr(iPON) = gtr(iPON) + (mortout/ratioCN)                ! [mmolN/m3/s]   
C       Production of alkalinity by NO3 uptake by macroalgae
        gtr(iALK) = gtr(iALK) + UptakeNout                       ! [meq/m3/s] 
C       Consumption of DIC by macroalgae growth
        growthout = Cmag * muMAG * 1E3 / mw_c                    ! [mmol-C/m3/s] 
        gtr(iDIC) = gtr(iDIC) - growthout                        ! [mmol-C/m3/s] 
C       Production oxygen by macroalgae growth
        gtr(iO2) = gtr(iO2) + (R_ON * UptakeNout)                ! [mmol-O2/m3/s]
      ELSE
        gtr(iMagB) = c0
        gtr(iMagC) = c0
        gtr(iMagQ) = c0
        UptakeNout = c0
        dEout = c0
        mortout = c0
        growthout = c0
      ENDIF      

                      

C     Extreme carbon limitations 
C     If expected C<Cmin then increase C, 
C     so that C=Cmin & compute the loss of C and N from structural mass)
C     IMPLEMENT THIS

             
      RETURN
      END SUBROUTINE
      
