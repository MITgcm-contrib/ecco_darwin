#include "DARWIN_OPTIONS.h"

#ifdef DARWIN_ALLOW_MACROALGAE
CBOP
C     ROUTINE: DARWIN_MACROALGAE
C     INTERFACE: ==========================================================
      SUBROUTINE DARWIN_MACROALGAE(       
     I     k,k_bottom,deltaT,
     I     drF,hFacC,
     I     t,hor_vel,parl,
     I     Hw,Tw,
     I     ptr,
     I     cum_length,
     I     length,
     I     total_length,
     I     biomass_up,biomass_down,
     U     gtr,
     U     excess_growth,excess_mort, 
     I     qin,
     O     prev_qmag, bladefrac,
     O     gT,gQ,gE,muMAG,
     O     mortMAG,UptakeNout,
     O     mortout,dEout,growthout)     
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
      INTEGER k, k_bottom             ! Vertical index
      _RL ptr(nDarwin)            ! Darwin tracers
      _RL gtr(nDarwin)            ! Darwin tendencies
      _RL t,hor_vel,parl          ! Darwin water temperature, horizontal velocity and PAR
      _RL deltaT                  ! Darwin time step duration
      _RL drF,hFacC
      _RL excess_mort             ! mortality to apply to upper level if biomass present
C malibu_v01 [A9]: quota of incoming biomass, its amount, mixed quota
      _RL qin, eg_in, qmix, dQmix
C malibu_v01 [A3]: dry weight per blade area, Macrocystis (G-MACMODS mag_species.py)
      _RL mp_spp_dry_sa
      parameter(mp_spp_dry_sa=58.0)
C malibu_v01 [A21]: kelp N:P (C:N:P = 550:35:1, CCS kelp_setup scripts)
      _RL R_PN_mag
      parameter(R_PN_mag=1.0/35.0)
      _RL excess_growth           ! excess of growth  due to carrying capacity     
      _RL biomass_up              ! biomass in above grid cell 
      _RL biomass_down            ! biomass in lowest grid cell (base)
C     Additional variables needed for physiological calculation in MACMODS
      _RL vQ                      ! Quota-limited uptake
C     Variables related to ocean flow and wave motion affecting nutrient uptake
      _RL wp_tmp                  ! Local water temperature
      _RL Tw                      ! Local wave period
      _RL Tw_s                    ! Maximum wave period
      _RL Hw                      ! Local significant wave height
      _RL n_length                ! number of iterations of wave-driven currents
      _RL visc                    ! kinematic viscosity 
      _RL magu                    ! Local horizontal flow
      _RL magu_m                  ! Local maximum horizontal flow
      _RL DBL                     ! Diffusive boundary layer thickness
      _RL vval                    ! Iterations for oscillatory flow
      _RL n                       ! Iterative number for oscillatory flow over vertical levels of the kelp forest
      _RL vNuTw_NO3               ! NO3 uptake factor
      _RL Uptake_NO3              ! Uptake of NO3
      _RL vNuTw_NH4               ! NH4 uptake factor
      _RL Uptake_NH4              ! Uptake of NH4
      _RL vNuTw_Urea              ! Urea uptake factor
      _RL Uptake_Urea             ! Uptake of Urea
      _RL UptakeNm3               ! Uptake of N in [mg/m3/s]
      _RL UptakeN                 ! Uptake of N in [mgN/gDW/s]
      _RL UptakeNout              ! Uptake of N in [mmol/m3/s] 
      _RL dEout                   ! Exudation for Darwin tendency
      _RL mortout                 ! Mortality for Darwin tendency
      _RL growthout               ! Growth for Darwin tendency
      _RL prev_qmag               ! Save Q for transfer
C malibu_v01 [K3]: blade fraction of this cell, for radtrans
      _RL bladefrac
      _RL total_length            ! total vertical frond length
      _RL cum_length              ! cumulative length
      _RL length                  ! vertical grid frond length
      _RL ratioCN                 ! Ratio C:N
      _RL Dm                      ! Molecular diffusivity of nitrate at 18degC
      _RL oscillatory             ! Oscillatory flow
      _RL flow                    ! Uni-directional flow
      _RL beta                    ! Oscillatory and Uni-directional flow
      _RL Bmag                    ! local biomass of Kelp 
      _RL Qmag                    ! local nitrogen cell quota
      _RL Cmag                    ! local biomass of Kelp in C
      _RL lambda_NO3                  ! Nutrient uptake related term
      _RL lambda_NH4                  ! Nutrient uptake related term
      _RL lambda_Urea                  ! Nutrient uptake related term      
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
      _RL NH4_mag                 ! NH4 conc. in mmol/m3
      _RL Urea_mag                ! Urea conc. in mmol/m3
      _RL gQ,gT,gE,gH             ! Growth limitation terms
      _RL muMAG,mortMAG,WP        ! Growth and mortality rates, wave mortality
      _RL Mwave, mwh              ! Wave mortality parameters            
C     Additional parameters
      _RL mp_growth_lim_type ! Growth limitation Type
      _RL break_Duarte_Ferr  ! Duarte frond breakage Type
      _RL break_Rodrig       ! Rodrigues frond breakage Type
      _RL mp_break_type      ! Frond breakage Type
      _RL mp_N_uptake_type   ! Type of nutrient uptake
      _RL mp_spp_Vmax_NO3_m3 ! Maximum uptake rate of nitrate converted to m3
      _RL mp_spp_Vmax_NH4_m3 ! Maximum uptake rate of ammonium converted to m3
      _RL mp_spp_Vmax_Urea_m3! Maximum uptake rate of urea converted to m3      
      _RL c0                 ! 0
      _RL c1                 ! 1
      _RL c4                 ! 4
      _RL pi                 ! PI
      _RL mw_n, mw_c         ! Molecular weight N and C
      _RL spd                ! Number of seconds in a day
      _RL R_ON               ! O:N ratio for respiration and consumption
      _RL growth,mort        ! local growth and mortality
      _RL frac_length        ! Fractional frond height (local length/total length)
      _RL blade2stipe        ! blade weight to stipe weight ratio
      _RL Bstipe             ! Stipe biomass
      _RL Bblade             ! Blade biomass
      _RL blade2tot          ! blade weight to total weight ratio     
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
C malibu_v01 [A5]: kg/m3 for wave power in kW/m (was 1.026)
      parameter(rho_sw=1025.0)            ! density of salt water (kg/m^3)
C     Assign value for frond breakage type      
      mp_break_type = break_Duarte_Ferr
C     Assign forcings
      NO3_mag=max(c0, ptr(iNO3))          ! mmol/m3
      NH4_mag=max(c0, ptr(iNH4))          ! mmol/m3
      Urea_mag=max(c0, 0.2*ptr(iDON))     ! mmol/m3      
      t_mag=t                             ! deg C   
      magu = hor_vel		          ! [m/s]     from ECCO-Darwin CCS model
C     growth Limitation model
C     if = 0 ==>  multiply all the limitaton factors
C     if > 0 ==>  independent light and nutrient limitation growth model
      mp_growth_lim_type=0
C     Setting the kind of nutrient uptake
      mp_N_uptake_type=1
      
C     ---------------------------------------------------------- 
C     MACROALGAE MODEL         
C     ---------------------------------------------------------- 

C     Parameters of macroalgae in data.darwin already converted to Darwin units
C     following these conversion factors:
C     Conversion parameters from https://doi.org/10.3389/fmars.2022.752951 and
C     https://doi.org/10.1038/s43247-023-00833-2 to Darwin units      
C     mp_spp_Gmax_cap = mp_spp_Gmax_cap/spd          ! from d-1 to s-1
C     mp_spp_Vmax = mp_spp_Vmax/spd                  ! from d-1 to s-1 [umol N/m2/s]
C     mp_spp_Vmax = mp_spp_Vmax/1E3                  ! from [umol N/m2/s] to [mmol N/m2/s] (Darwin unit)
C     mp_spp_E = mp_spp_E/spd                        ! from d-1 to s-1
C     mp_spp_death = mp_spp_death/spd                ! from d-1 to s-1
C     mp_spp_Ks_NO3 = mp_spp_Ks_NO3/1E3              ! from [umol N/m3] to [mmol N/m3] (Darwin unit)     
C     mp_spp_PARs = mp_spp_PARs*4.57                 ! from [W/m2] to [microEinstein/m2/s]
C     mp_spp_PARc = mp_spp_PARc*4.57                 ! from [W/m2] to [microEinstein/m2/s]
C     Convert carrying capacity from m to dry biomass in gDW/m3 

C     ----------------------------------------------------------       
C     Lambda NO3 function         
C     ----------------------------------------------------------   
C     Unit conversions, minima
C malibu_v01 [A6]: floor of 1 m/day as in G-MACMODS (was 1 m/s)
      magu_m = max(1.0/spd,magu)                                ! [m/s]
      Tw_s = max(0.01,Tw)                                       ! [s]       mean wave period
C     some parameters and variables
      n_length = 25                                             ! number of iterations of wave-driven currents;
                                                                ! selected based on number of iterations attaining 95% of
                                                                ! maximum value for n_length = 1000
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
      lambda_NO3 = c1 + (mp_spp_Vmax_NO3 / (beta*mp_spp_Ks_NO3)) - 
     & (NO3_mag/mp_spp_Ks_NO3)
      lambda_NH4 = c1 + (mp_spp_Vmax_NH4 / (beta*mp_spp_Ks_NH4)) - 
     & (NH4_mag/mp_spp_Ks_NH4)
      lambda_Urea = c1 + (mp_spp_Vmax_Urea / (beta*mp_spp_Ks_Urea)) - 
     & (Urea_mag/mp_spp_Ks_Urea)
                                               ! [no unit]      
C     Quota-limited uptake: maximum uptake when Q is minimum and
C     approaches zero as Q increases towards maximum; Possible that Q
C     is greater than Qmax. Set any negative values to zero.      
C malibu_v01 [A15]: 0 = linear, 1 = "50" curve, as in mag_calc.f90
      IF (mp_N_uptake_type .EQ.0) THEN
       vQ = (mp_spp_Qmax-ptr(iMagQ) )/(mp_spp_Qmax-mp_spp_Qmin)       ! [0-1]
       vQ = max(c0,vQ)                                                ! [0-1]
       vQ = min(c1,vQ)                                                ! [0-1]
      ELSEIF (mp_N_uptake_type .EQ.1) THEN
       vQ = 1.0 - 1.0/(1.0+(max(0.0,mp_spp_Qmax-ptr(iMagQ) ))*        ! [0-1]
     & 50.0/(mp_spp_Qmax-mp_spp_Qmin))
       vQ = max(c0,vQ)                                                ! [0-1]
       vQ = min(c1,vQ)                                                ! [0-1]
      ENDIF
      
C     Below is what we call "Uptake Factor." It varies betwen 0
C     and 1 and includes kinetically limited uptake and
C     mass-transfer-limited uptake (oscillatory + uni-directional flow)
      vNuTw_NO3 = NO3_mag / (mp_spp_Ks_NO3 * 
     & ((NO3_mag/mp_spp_Ks_NO3)+ 
     & 0.5 *(lambda_NO3+sqrt(lambda_NO3**2 + c4 * 
     & (NO3_mag/mp_spp_Ks_NO3)))))                                    ! [0-1]   
      vNuTw_NO3 = max(c0,vNuTw_NO3)                                   ! [0-1]
      vNuTw_NO3 = min(c1,vNuTw_NO3)                                   ! [0-1]
      
      IF(mp_spp_allN)THEN
         vNuTw_NH4 = NH4_mag / (mp_spp_Ks_NH4 * 
     &    ((NH4_mag/mp_spp_Ks_NH4)+ 
     &    0.5 *(lambda_NH4+sqrt(lambda_NH4**2 + c4 * 
     &    (NH4_mag/mp_spp_Ks_NH4)))))                                 ! [0-1]   
         vNuTw_NH4 = max(c0,vNuTw_NH4)                                ! [0-1]
         vNuTw_NH4 = min(c1,vNuTw_NH4)                                ! [0-1]
      
         vNuTw_Urea = Urea_mag / (mp_spp_Ks_Urea * 
     &    ((Urea_mag/mp_spp_Ks_Urea)+ 
     &    0.5 *(lambda_Urea+sqrt(lambda_Urea**2 + c4 * 
     &    (Urea_mag/mp_spp_Ks_Urea)))))                               ! [0-1]   
         vNuTw_Urea = max(c0,vNuTw_Urea)                              ! [0-1]
         vNuTw_Urea = min(c1,vNuTw_Urea)                              ! [0-1]
      ELSE
         vNuTw_Urea = 0.0
         vNuTw_NH4 = 0.0       
      ENDIF                                        
                                           
C     Uptake Rate [mg N/gDW/s]
C     Nutrient Uptake Rate = Max Uptake * v[Ci,u,Tw] * vQ
C malibu_v01 [A3]: Vmax per m2 blade -> per gDW (was / cell thickness)
      mp_spp_Vmax_NO3_m3 = mp_spp_Vmax_NO3 / mp_spp_dry_sa      ! [mmolN/gDW/s]
      mp_spp_Vmax_NH4_m3 = mp_spp_Vmax_NH4 / mp_spp_dry_sa      ! [mmolN/gDW/s]
      mp_spp_Vmax_Urea_m3 = mp_spp_Vmax_Urea / mp_spp_dry_sa    ! [mmolN/gDW/s]
      
C     Compute weight distribution in blades and stipesaccording to Nyman et al 1993
C     equation 5b, slightly modified to avoid curvature at low fractional frond height (<= 0.3)     
      IF (mp_spp_vertical)THEN
C malibu_v01: guard total_length = 0 (no kelp in the column)
         frac_length = c0
         IF (total_length.GT.c0) frac_length =
C malibu_v01 [A14]: height of the cell centre above the holdfast / total
C malibu_v01 [A14]: (cum_length is summed from the surface down)
     &     min(c1,max(c0,(total_length-cum_length+0.5*length)
     &     /total_length))
C malibu_v01 [A14]: +23.39 f^2: minimum 0.34 at f = 0.29 (was -23.39)
         blade2stipe = 2.32 - 13.6*frac_length + 23.39*frac_length**2
         IF (frac_length.LE.0.3) THEN 
                 blade2stipe = 0.34 
         ENDIF
         Bstipe = ptr(iMagB) / ((blade2stipe)+1)                      ! Stipe biomass [gDW/m3]
         Bblade = ptr(iMagB) - Bstipe                                 ! Blade biomass [gDW/m3]
C malibu_v01: guard Bmag = 0
         blade2tot = c0
         IF (ptr(iMagB).GT.c0) blade2tot =
     &     min(c1,max(c0,Bblade/ptr(iMagB)))                ! blade weight to total weight ratio
      ELSE
         Bblade = 0.0
         Bstipe = 0.0
         blade2tot = 0.0                                  
      ENDIF 
      
      IF(ptr(iMagB).GT.0.0)THEN  
C malibu_v01 [A3]: uptake = Vmax/dry_sa * v * vQ * B (local B), capped by
C malibu_v01 [A3]: the nutrient; kelp gains exactly the N the water loses
        Uptake_NO3 = min(mp_spp_Vmax_NO3_m3 * vNuTw_NO3 * vQ
     &   * ptr(iMagB), NO3_mag/deltaT)                    ! [mmolN/m3/s]
        Uptake_NH4 = min(mp_spp_Vmax_NH4_m3 * vNuTw_NH4 * vQ
     &   * ptr(iMagB), NH4_mag/deltaT)                    ! [mmolN/m3/s]
        Uptake_Urea = min(mp_spp_Vmax_Urea_m3 * vNuTw_Urea * vQ
     &   * ptr(iMagB), Urea_mag/deltaT)                   ! [mmolN/m3/s]
        UptakeNm3 = (Uptake_NO3 + Uptake_NH4 + Uptake_Urea) * mw_n  ! [mgN/m3/s]
        UptakeN = max(c0, UptakeNm3 / ptr(iMagB))       ! [mgN/gDW/s]
      ELSE
        Uptake_NO3 = c0
        Uptake_NH4 = c0
        Uptake_Urea = c0
        UptakeN = c0
      ENDIF     
      
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
C     gQ = (Q(i,j) - params(mp_spp_Qmin)) / Q(i,j)              ! Droop equation      
C malibu_v01: guard Qmag = 0
      gQ = c0
      IF (ptr(iMagQ).GT.c0)
     & gQ = ((ptr(iMagQ)  - mp_spp_Qmin) / ptr(iMagQ) ) * 
     & (mp_spp_Qmax/(mp_spp_Qmax-mp_spp_Qmin))                  ! Droop scaled from 0-1
C     gQ = (Q(i,j) - Qmin) / (Qmax-Qmin)                        ! Freider et al.
      gQ = max(c0,gQ)                                           ! [0-1]
      gQ = min(c1,gQ)                                           ! [0-1]
      
C     Temperature-dependent growth          
      IF (t_mag .GE. mp_spp_Topt1) THEN          
         IF (t_mag .LE. mp_spp_Topt2) THEN
            gT = 1.0                                            ! [0-1]
         ELSE  
            gT = exp(-mp_spp_K2*(t_mag-mp_spp_Topt2)**2)        ! [0-1]
         ENDIF          
      ELSE              
         gT = exp(-mp_spp_K1*(t_mag-mp_spp_Topt1)**2)           ! [0-1]      
      ENDIF 
      gT = max(0.0,gT)                                          ! [0-1]
      gT = min(1.0,gT)                                          ! [0-1]
    
C     Light limitation
      IF (parl .LT. mp_spp_PARc) THEN              
          gE = c0                                               ! [0-1]
      ELSEIF (parl .GT. mp_spp_PARs) THEN              
          gE = c1                                               ! [0-1]
      ELSE              
          gE = ((parl-mp_spp_PARc)/(mp_spp_PARs-mp_spp_PARc))*
     & exp(-((parl-mp_spp_PARc)/(mp_spp_PARs-mp_spp_PARc))+c1)  ! [0-1]  
      ENDIF
                          
C     Growth rate
C     interacting limitation terms growth model
      IF (mp_growth_lim_type.EQ.0.AND.ptr(iMagB).GT.0.0) THEN
    	 muMAG =  mp_spp_Gmax_cap * gT * gE * gQ                ! [per s]
      ELSEIF(mp_growth_lim_type.EQ.1.AND.ptr(iMagB).GT.0.0)THEN
C     independent light and nutrient limitation growth model
    	 muMAG =  mp_spp_Gmax_cap * gT * min(gE,gQ)             ! [per s]   
      ELSE
         muMAG = 0.0          
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
         WP = (((rho_sw * (9.8**2)) /(64.0*pi)) * 
     &    (mwh**2) * wp_tmp) / 1E3 ! [kW]   
         Mwave=((2.3*1E-4) * WP + (2.2*1E-3)) * mp_wave_mort_factor        ! [per d]
      ELSEIF (mp_break_type.EQ.break_Rodrig) THEN
C     Death rate scales with Hw (Hw, surface wave height)
         Mwave=(mp_spp_death*spd * Hw) * mp_wave_mort_factor               ! [per d]
      ENDIF
C     M_wave + general Mortality rate; [per s]                                                  
C malibu_v01 [A23]: background mortality acts on the whole plant, as in G-MACMODS
C malibu_v01 [A23]: (Arzeno-Soltero et al. 2023; was mp_spp_death*blade2tot, so stipes never
C malibu_v01 [A23]: died and losses were ~half the SBC LTER winter rate, 0.022 d-1)
C malibu_v01 [SM1]: background mortality of this time step (seasonal when mp_spp_death_mon is set,
C     else mp_spp_death; set in DARWIN_FORCING)
      mortMAG = mp_spp_death_now + Mwave/spd                              ! [per s]
C     Compute exudation
      IF(ptr(iMagB).GT.0.0)THEN
        dE = max(c0,mp_spp_E * (ptr(iMagQ) -mp_spp_Qmin) * deltaT)         ! [mg-N/gDW]
      ELSE
        dE = 0.0
      ENDIF
      
C     ---------------------------------------------------------- 
C     Vertical transfer
C     ----------------------------------------------------------      
C     Carrying capacity for vertical growth of macroalgae
C malibu_v01 [A9]: incoming growth (already in this cell's units)
      eg_in = max(c0, excess_growth)
      IF (.NOT.mp_spp_vertical) eg_in = c0
      qmix  = ptr(iMagQ)
      dQmix = c0
      IF (ptr(iMagB)+eg_in .GT. c0 .AND. eg_in.GT.c0) THEN
        qmix  = (ptr(iMagB)*ptr(iMagQ) + eg_in*qin)
     &        / (ptr(iMagB)+eg_in)
        dQmix = qmix - ptr(iMagQ)
      ENDIF
      IF (mp_spp_vertical) THEN   
        IF (length.GE.(drF*hFacC).AND.
     &   total_length.LT.mp_spp_maxlength.AND.
     &   k.GT.1) THEN
C       If local length exceeds grid cell thickness, transfer growth and mortality in upper levels.
C       Growth is computed locally with associated uptake of N, DIC, production of ALK and release of N and C with exudation.
C       Growth is then applied to first grid cell with length lower than grid cell thickness or grid cell at the surface.
C       Mortality is computed locally but without release of N and C. Mortality is applied to first grid cell with length lower 
C       than grid cell thickness or grid cell at the surface, where N and C is released.                 
          growth = ptr(iMagB)  * muMAG * deltaT                           ! Growth is computed for feedback with C, N and ALK but not 
C malibu_v01 [A9]: outgoing quota: local growth (diluted Q) + pass-through
          prev_qmag = c0
          IF (growth+eg_in.GT.c0) prev_qmag =
     &      (growth*ptr(iMagQ)/(c1+muMAG*deltaT) + eg_in*qin)
     &      / (growth+eg_in)
                                                                          ! applied to biomass here [gDW/m3]
          excess_growth = (((ptr(iMagB)  * muMAG * deltaT) +              ! Add transfer of growth from lower levels [gDW/m3]
     &     excess_growth)) 
          mort = c0                                                       ! Local mortality is set to zero but added to mortality to transfer 
                                                                          ! bellow
          excess_mort = min(ptr(iMagB) ,
     &     (ptr(iMagB)  * mortMAG * deltaT) + 
     &     excess_mort)                                                   ! Local mortality is computed here and added to mortality from lower 
                                                                          ! levels and to transfer to upper levels                         
          dBmagdt = 0.0                                                   ! Biomass is not changing [gDW/m3]
          dQgrowthdt = ptr(iMagQ)  * (1.0/(1.0 + muMAG * deltaT) - 1.0)   ! [mg-N/gDW]   
          dQmagdt = (UptakeN * deltaT + dQgrowthdt - dE)                  ! [mg-N/gDW]   
        ELSEIF(total_length.GE.mp_spp_maxlength) THEN
C       If total frond length exceeds max frond length, then growth is stopped at every level. Only mortality is computed 
C       and applied to first grid cell with length lower than grid cell thickness or grid cell at the surface.
          growth = c0                                                     ! No growth [gDW/m3]                                                   
C malibu_v01 [A11]: incoming growth stays here (was discarded)
          excess_growth = c0                                              ! No vertical transfer of growth [gDW/m3]
          prev_qmag = c0                                                  ! No vertical transfer of quota [mg-N/gDW]  
          mort = c0                                                       ! Local mortality is set to zero but added to mortality to transfer 
                                                                          ! bellow
          excess_mort = min(ptr(iMagB) ,
     &     (ptr(iMagB)  * mortMAG * deltaT) + 
     &     excess_mort)                                                   ! Local mortality is computed here and added to mortality from lower
                                                                          ! levels and to transfer to upper levels                         
          IF(biomass_up.EQ.0.0)THEN
C         Check if local cell is the uppermost one with biomass so mortality can be applied                  
           mort = excess_mort                                             ! Local mortality account for local mortality plus the transfer of 
                                                                          ! mortality
           excess_mort = 0.0                                              ! Reset transfer of mortality
           dBmagdt = -mort + eg_in                                        ! Biomass only changes with mortality      
          ELSE
           dBmagdt = eg_in                                                ! No change of biomass if there is biomass in grid cell above that has
                                                                          ! to be empty first
          ENDIF 
          dQgrowthdt = 0.0                                                ! No change of quota associated with growth [mg-N/gDW]
          dQmagdt = (UptakeN * deltaT + dQgrowthdt - dE) + dQmix
        ELSE
C       Local length does not exceed grid cell thickness or max length not reached or we are at the surface.
C       Accounts for all the transfer of growth and mortality from lower grid cells.                         
          growth = ptr(iMagB)  * muMAG * deltaT                           ! Local growth [gDW/m3]         
          mort = min(ptr(iMagB) ,(ptr(iMagB)  * mortMAG * deltaT) +       ! Mortality accounts for local mortality and mortality from lower grid 
                                                                          ! cells
     &     excess_mort)
C malibu_v01 [A9]: incoming biomass mixes its quota in (was: Q overwritten)
          dBmagdt = growth - mort + eg_in                         ! [gDW/m3]
          dQgrowthdt = qmix * (1.0/(1.0 + muMAG * deltaT) - 1.0)  ! [mg-N/gDW]
          dQmagdt = (UptakeN * deltaT + dQgrowthdt - dE) + dQmix  ! [mg-N/gDW]
          excess_mort = c0                                                ! Reset transfer of mortality
          excess_growth = c0                                              ! Reset transfer of growth
          prev_qmag = c0                                                  ! Reset transfer of quota
        ENDIF              
      ELSE
C     If vertical transfer is not activated
         growth = ptr(iMagB)  * muMAG * deltaT                            ! [gDW/m3]  
         mort = ptr(iMagB)  * mortMAG * deltaT                            ! [gDW/m3]  
         excess_growth = c0                                               ! No transfer
         excess_mort = c0                                                 ! No transfer
         prev_qmag = c0                                                   ! No transfer
         dBmagdt = growth - mort                                          ! No transfer
         dQgrowthdt = ptr(iMagQ)  * (1.0/(1.0 + muMAG * deltaT) - 1.0)    ! [mg-N/gDW]
         dQmagdt = (UptakeN * deltaT + dQgrowthdt - dE)                   ! [mg-N/gDW]   
      ENDIF
             
C     ---------------------------------------------------------- 
C     Apply tendencies
C     ----------------------------------------------------------                                            
        gtr(iMagB) = dBmagdt/deltaT                                       ! [gDW/m3/s]                                   
        gtr(iMagC) = (dBmagdt/deltaT) * mp_spp_carbon                     ! [gC/m3/s]    
        gtr(iMagQ) = dQmagdt/deltaT                                       ! [mg-N/gDW/s]
        
C     ---------------------------------------------------------- 
C     Apply feedback to C, N and ALK (only if there is biomass)
C     ----------------------------------------------------------              
        IF(ptr(iMagB).GT.0.0.AND.
     &   ptr(iMagQ).GT.0.0)THEN      
C malibu_v01 [A7]: molar C:N (was mass C:N used on mmol)
          ratioCN = (ptr(iMagC) *1E3)/(ptr(iMagB) *ptr(iMagQ) )
     &      * mw_n/mw_c                                    ! [molC/molN]
        ELSE
          ratioCN = 16.6
        ENDIF
                
        UptakeNout = (Uptake_NO3 + Uptake_NH4 + Uptake_Urea)              ! [mmol-N/m3/s]
        dEout = (dE/deltaT) * ptr(iMagB)  / mw_n                          ! [mmol-N/m3/s]
        mortout = mort/deltaT * mp_spp_carbon * 1E3 / mw_c                ! [mmol-C/m3/s]
        growthout = growth/deltaT * mp_spp_carbon * 1E3 / mw_c            ! [mmol-C/m3/s]
 
        IF(ptr(iMagB).GT.0.0)THEN        
C       Remove N used by macroalgae        
        gtr(iNO3) = gtr(iNO3) - Uptake_NO3                       ! [mmol-N/m3/s] 
        gtr(iNH4) = gtr(iNH4) - Uptake_NH4                       ! [mmol-N/m3/s] 
        gtr(iDON) = gtr(iDON) - Uptake_Urea                       ! [mmol-N/m3/s] 
C       Add macroalgae exudation to DOC and DON      
C malibu_v01 [A8]: G-MACMODS exudes N only (was + dEout*ratioCN to DOC)
        gtr(iDON) = gtr(iDON) + dEout                            ! [mmol-N/m3/s]
C       Add macroalgae mortality to POC and PON 
        gtr(iPOC) = gtr(iPOC) + mortout                          ! [mmolC/m3/s]   
        gtr(iPON) = gtr(iPON) + (mortout/ratioCN)                ! [mmolN/m3/s]   
C       Production of alkalinity by NO3 uptake by macroalgae
        gtr(iALK) = gtr(iALK) + Uptake_NO3                       ! [meq/m3/s] 
C       Consumption of DIC by macroalgae growth
        gtr(iDIC) = gtr(iDIC) - growthout                        ! [mmol-C/m3/s] 
C       Production oxygen by macroalgae growth
C malibu_v01 [A4]: 1 O2 per C fixed + 2 O2 per NO3 assimilated
C malibu_v01 [A4]: (was R_ON = 16/170 O2 per N taken up)
        gtr(iO2) = gtr(iO2) + growthout + 2.0*Uptake_NO3     ! [mmol-O2/m3/s]
C malibu_v01 [A21]: kelp P at N:P = 35 (uptake, exudation, mortality)
        gtr(iPO4) = gtr(iPO4) - R_PN_mag*UptakeNout
        gtr(iDOP) = gtr(iDOP) + R_PN_mag*dEout
        gtr(iPOP) = gtr(iPOP) + R_PN_mag*(mortout/ratioCN)
        ENDIF
             
C malibu_v01 [K3]: blade fraction for the light code (1 without allometry)
      bladefrac = c1
      IF (mp_spp_vertical .AND. ptr(iMagB).GT.c0) bladefrac = blade2tot
      RETURN
      END SUBROUTINE
#endif
      
