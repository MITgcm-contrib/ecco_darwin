#include "DARWIN_OPTIONS.h"

CBOP
C !ROUTINE: DARWIN_MACROALGAE
C !INTERFACE: ==========================================================
      SUBROUTINE DARWIN_MACROALGAE(
C    I     t,hor_vel,parl,ptr,lat,dayyr,
     I     t,hor_vel,parl,ptr,
     I     Hw,Tw,
     U     gtr,
C    O     Aalg,Nres,Cres,
     O     gT,gQ,gH,gE,muMAG,
     O     mortMAG)
C     I     ,dBmagdt,dNmagdt,

C t:         temperature  [C]
C hor_vel:   horizontal velocity [m s-1]
C parl: photosynthetic radiation [mu Ein s-1 m-2]
C ptr(ndarwin): tracer concentrations 
C lat: latitude [deg]
C dayyr: day in a year [0-365/366]
 
C gtr(ndarwin): tracer tendencies (inout)

C macroalgage relevant tracers in ptr
C ptr(MaNum) .... number of macroalgage [km^-3]
C ptr(MaFrA) .... concentration of macroalgage structural (frond) area [dm^2/km^-3]
C ptr(MaC) .... concentration of C in macroalgage (storage) [uM]
C ptr(MaN) .... concentration of N in macroalgage (storage) [uM]
C
C gtr = tendency dptr/dt due to physical processes [ptr s-1]

C Based on Broch et al and
C Strong-Wright paper, https://github.com/jagoosw/SugarKelp.jl

      IMPLICIT NONE
#include "DARWIN_PARAMS.h"
#include "DARWIN_SIZE.h"
#include "DARWIN_INDICES.h"
#include "DARWIN_DIAGS.h"

C Inputs/outputs

      _RL ptr(nDarwin)
      _RL gTr(nDarwin)
      _RL t,hor_vel,parl
      _RL N_NO3    

C Additional variables needed for physiological calculation in MACMODS

      _RL vQ
   
C Variables related to ocean flow and wave motion affecting nutrient uptake

      _RL wp_tmp
      _RL Tw
      _RL Hw
      _RL TW_s
      _RL n_length ! number of vertical levels of the kelp forest  
      _RL visc
      _RL magu
      _RL magu_m
      _RL DBL
      _RL lambda_NO3
      _RL lambda_mcd
      _RL vval
      _RL n
      _RL vNuTw_NO3
      _RL NO3_u
      _RL Uptake_NO3
      _RL Q
      _RL UptakeN
      _RL Dm
      _RL oscillatory
      _RL flow
      _RL beta
      _RL Bmag ! biomass of Kelp 
      _RL Nmag
C      
      _RL temp_lim
      _RL lambda

C Additional MACMODS variables
      _RL dE ! Exudation rate
C      _RL vQ
      _RL Upc1,Upc2
C Tendency terms
      _RL dBmagdt
      _RL dNmagdt
      _RL dCmagdt            

C Key variables to store forcing functions from ECCO-Darwin
C mag = macroalgae

      _RL t_mag ! water temperature in degC
      _RL par_mag ! PAR in Watts/m2      
      _RL NO3_mag ! NO3 conc. in mmol/m3
      _RL chl_plk ! chlrophyll from plankton (test variable)

C Factors affecting growth  
        
C Growth, nitrogen movement from Ns to Nf = umax*gQ*gT*gE*gH; [per day]
C Output:
C Growth, [h-1]

C gQ, quota-limited growth
C from Wheeler and North 1980 Fig. 2

C gT, temperature-limited growth
C piecewise approach taken from Broch and Slagstad 2012 (for sugar
C kelp) and optimized for Macrocystis Pyrifera

C gE, light-limited growth
C from Dean and Jacobsen 1984

C gH, carrying capacity-limited growth

C muMAG, macroalgal growth      

      _RL gQ,gT,gE,gH,muMAG,mortMAG,WP
      _RL Mwave, mwh

      _RL par_watts,chlmin,mp_farm_depth
      _RL A,B,Bnew
      _RL dNs
      
C Additional parameters

      _RL mp_growth_lim_type 
      _RL break_Duarte_Ferr
      _RL break_Rodrig
      _RL mp_break_type ! frond breakage Type
      _RL mp_N_uptake_type ! type of nutrient uptake
      _RL mp_wave_mort_factor

      _RL c0    
      _RL c1 
      _RL c2 
      _RL c3  
      _RL c4 
      _RL c10 
      _RL c1000 
      _RL p001 
      _RL p5 
      _RL pi 

C Additional MACMODS parameters

      _RL mw_n		! molecular weight N
      _RL sphr     ! number of seconds in an hour
      _RL spd      ! number of seconds in a day
      _RL dpy      ! number of days in a year
      _RL spy      ! number of seconds in a year
      _RL hrps     ! number of hours in a second
      _RL dps      ! number of days in a second
      _RL ypd      ! number of years in a day
      _RL yps      ! number of years in a second
      _RL cmperm   ! cm per meter
      _RL mpercm   ! meters per cm      

C Physical constants used in MACMODS
 
      _RL vonkar       ! von Karman constant
      _RL K_Boltz      ! Boltzmann constant (eV/K)
      _RL rho_sw       ! density of salt water (g/cm^3)
      _RL epsC         ! small C concentration (mmol C/m^3)
      _RL epsTinv      ! small inverse time scale (1/year) (1/sec)
      _RL molw_Fe      ! molecular weight of iron (gFe / mol Fe)
      _RL R13C_std     ! actual 13C/12C PDB standard ratio (Craig, 1957) = 1123.72e-5_r8
      _RL R14C_std     ! actual 14C/12C NOSAMS standard ratio = 11.76e-13_r8
         
C time constant
      _RL day_in_s
      _RL h_in_s
	
C number of seconds in day	
      parameter(day_in_s=86400.)
      
C number of s in hour
      parameter(h_in_s=3600.)

C Additional MACMODS parameters with assigned values 

      parameter(break_Duarte_Ferr = 0)		! more-complex breakage
      parameter(break_Rodrig       = 1  )

      parameter(c0=0.0)
      parameter(c1=1.0)
      parameter(c2=2.0)
      parameter(c3=3.0)
      parameter(c4=4.0)
      parameter(c10=10.0)
      parameter(c1000  = 1000.0)
      parameter(p001   =    0.001)
      parameter(p5     =    0.5)
      parameter(pi     =    3.14159265358979323846)

      parameter(mw_n      = 14.00672)	! molecular weight N
      parameter(sphr      = 3600.0)     ! number of seconds in an hour
      parameter(spd       = 86400.0)    ! number of seconds in a day
      parameter(dpy       = 365.0)      ! number of days in a year
      parameter(spy       = dpy*spd)    ! number of seconds in a year
      parameter(hrps      = c1 / sphr)  ! number of hours in a second
      parameter(dps       = c1 / spd)   ! number of days in a second
      parameter(ypd       = c1 / dpy)   ! number of years in a day
      parameter(yps       = c1 / spy)   ! number of years in a second
      parameter(cmperm    = 100.0)      ! cm per meter
      parameter(mpercm    = 0.01)       ! meters per cm

C Additional physical constants for MACMODS

      parameter(vonkar    =   0.4)               ! von Karman constant
      parameter(K_Boltz   =   8.617330350*1E-5) ! Boltzmann constant (eV/K)
      parameter(rho_sw    =   1.026)             ! density of salt water (g/cm^3)
      parameter(epsC      =   1.0*1E-8)          ! small C concentration (mmol C/m^3)
      parameter(epsTinv   =   3.17*1E-8)         ! small inverse time scale (1/year) (1/sec)
      parameter(molw_Fe   =  55.845)             ! molecular weight of iron (gFe / mol Fe)
      parameter(R13C_std  =   1.0)               ! actual 13C/12C PDB standard ratio (Craig, 1957) = 1123.72e-5_r8
      parameter(R14C_std =    1.0)               ! actual 14C/12C NOSAMS standard ratio = 11.76e-13_

C Assign value for frond breakage type      

      mp_break_type = break_Duarte_Ferr

      NO3_mag=ptr(iNO3) ! mmol/m3

C Exclude any value less than 0
      IF (NO3_mag .LT. 0 ) THEN 
              
         NO3_mag=0.0  
             
      ENDIF

      t_mag=t ! deg C
      PAR_mag=parl	! Watts/m2

C Now fake value of chl from plakton for test 
      chl_plk=0.1	! mg/m3      

C growth Limitation model
C if = 0 ==>  multiply all the limitaton factors
C if > 0 ==>  independent light and nutrient limitation growth model

      mp_growth_lim_type=0

C Setting the kind of nuts uptake

      mp_N_uptake_type=1

C Lambda NO3 function
         
C        Tw = 10;	! wave period in seconds (typical number)
      n_length = 1 
      visc = 1E-6* 86400
      Dm = (18*3.65*1E-11 + 9.72*1E-10) * 86400

      magu = hor_vel		! m/s (from ECCO-Darwin CCS model)

C Unit conversions, minima

      magu_m = max(1.0,magu*86400) ! m/day
      Tw_s = max(0.01,Tw)/86400.0  ! whaa...?

      DBL = 10 * (visc / (sqrt(mp_spp_CD) * abs(magu_m)));

C For now, set to 1
      n_length=1
         
C Flow
C		1. Oscillatory flow
      DO n = 1,n_length    
          
        vval = (1-exp((-Dm * n**2 * pi**2 *Tw_s)/(2.0*DBL**2)))
     & /(n**2 * pi**2) 
         
      ENDDO
C        oscillatory = ((4.0*DBL)/Tw_s) * sum(vval)		! original line
      oscillatory = ((4.0*DBL)/Tw_s) * vval

C       2. Uni-directional Flow
      flow = Dm / DBL

      beta = flow + oscillatory

      lambda_NO3 = c1 + (mp_spp_Vmax / (beta*mp_spp_Ks_NO3)) - 
     & (NO3_u/mp_spp_Ks_NO3)
      lambda=lambda_NO3

C Temperature-dependent growth  
        
      IF (t_mag .GE. mp_spp_Topt1) THEN 
         
         IF (t_mag .LE. mp_spp_Topt2) THEN
            temp_lim = 1.0
         ELSE 
            temp_lim = exp(-mp_spp_K2*(t_mag-mp_spp_Topt2)**2)
         ENDIF 
         
      ELSE
              
        temp_lim = exp(-mp_spp_K1*(t_mag-mp_spp_Topt1)**2)
        
      ENDIF 
      temp_lim = max(0.0,temp_lim)
      temp_lim = min(1.0,temp_lim)
      gT=temp_lim

C     Bmag= mp_spp_seed
      Q = mp_spp_Qmin + NO3_mag*(mp_spp_Qmax-mp_spp_Qmin)/35
      
C Convert to micromol/m3
      Q = mp_spp_Qmin + (NO3_mag*1E3)*(mp_spp_Qmax-mp_spp_Qmin)/35

      Bmag = ptr(iMagB)
      Nmag = ptr(iMagN)
      Q = Nmag/Bmag 

      IF (mp_N_uptake_type .EQ.0) THEN

       vQ = (mp_spp_Qmax-Q)/(mp_spp_Qmax-mp_spp_Qmin)
       vQ = max(c0,vQ)
       vQ = min(c1,vQ)

      ELSEIF (mp_N_uptake_type .EQ.0) THEN

       vQ = 1.0 - 1.0/(1.0+(max(0.0,mp_spp_Qmax-Q))*
     & 50.0/(mp_spp_Qmax-mp_spp_Qmin))
       vQ = max(c0,vQ)
       vQ = min(c1,vQ)

      ENDIF

C Below is what we call "Uptake Factor." It varies betwen 0
C and 1 and includes kinetically limited uptake and
C mass-transfer-limited uptake (oscillatory + uni-directional flow)

      NO3_u = NO3_mag*1000.0        
      vNuTw_NO3 = NO3_u / (mp_spp_Ks_NO3 * ((NO3_u/mp_spp_Ks_NO3)+ 
     & 0.5 *(lambda+sqrt(lambda**2 + c4 * (NO3_u/mp_spp_Ks_NO3)))))     
      vNuTw_NO3 = max(c0,vNuTw_NO3)        
      vNuTw_NO3 = min(c1,vNuTw_NO3)   

C Uptake Rate [mg N/g(dry)/d]
C Nutrient Uptake Rate = Max Uptake * v[Ci,u,Tw] * vQ
C converted from umol N/m2/d -> mg N/g(dry)/d by 14.0067 / 1e3

      Uptake_NO3 = mp_spp_Vmax * vNuTw_NO3 * vQ ! [umol/m2/d]         
      Uptake_NO3 = Uptake_NO3 * mw_n / 1E3         
      UptakeN = Uptake_NO3
      UptakeN = UptakeN / mp_spp_dry_sa

C Growth, nitrogen movement from Ns to Nf = umax*gQ*gT*gE*gH; [per day]
C Output:
C   Growth, [h-1]
C   gQ, quota-limited growth
C       from Wheeler and North 1980 Fig. 2
C   gT, temperature-limited growth
C       piecewise approach taken from Broch and Slagstad 2012 (for sugar
C       kelp) and optimized for Macrocystis pyrifera
C   gE, light-limited growth
C       from Dean and Jacobsen 1984
C   gH, carrying capacity-limited growth

C 		gQ = (Q(i,j) - params(mp_spp_Qmin)) / Q(i,j) % Droop equation
      gQ = (Q - mp_spp_Qmin) / Q * mp_spp_Qmax/(mp_spp_Qmax-mp_spp_Qmin) ! Droop scaled from 0-1

C 		gQ = (Q(i,j) - Qmin) / (Qmax-Qmin) ! Freider et al.
      gQ = max(c0,gQ)
      gQ = min(c1,gQ)

      par_watts = par_mag

      IF (par_watts .LT. mp_spp_PARc) THEN
              
          gE = c0
          
      ELSEIF (par_watts .GT. mp_spp_PARs) THEN
              
          gE = c1
          
      ELSE
              
          gE = (par_watts-mp_spp_PARc)/(mp_spp_PARs-mp_spp_PARc)*
     & exp(-(par_watts-mp_spp_PARc)/(mp_spp_PARs-mp_spp_PARc)+c1)
     
      ENDIF

      Bnew = Bmag *mp_spp_line_sep ! converting from g/m2 to g/m
      A = mp_spp_kcap_rate/mp_spp_kcap**(-0.75)
      gH = A*Bnew**(-0.75)

C Compute growth

      IF (mp_growth_lim_type .EQ. 0) THEN
              
                 ! interacting limitation terms growth model
    	 muMAG =  min(gH,mp_spp_Gmax_cap) * gT * gE * gQ
         
      ELSE 
         ! independent light and nutrient limitation growth model
         ! this is the preferred formula
         ! used by MACMODS group
    	 muMAG =  min(gH,mp_spp_Gmax_cap) * gT * min(gE,gQ)
         
      ENDIF

C Compute mortality rate due to wave activity      
C Wave height is computed by the ERA wave model

C d_wave = frond loss due to waves; dependent on Hs, significant
C wave height [m]; Rodrigues et al. 2013 demonstrates linear relationship
C between Hs and frond loss rate in Macrocystis [d-1] (continuous)
C Duarte and Ferreira (1993) find a linear relationship between wave
C power and mortality in red seaweed        

C Assign here the value of wave height

      mwh  = Hw    ! Assign wave height 
      wp_tmp = Tw  ! Assign wave period

      IF (mp_break_type.EQ.break_Duarte_Ferr ) THEN

         WP = (rho_sw * (9.8**2)) /(64.0*pi) * (mwh**2) *wp_tmp / 1E3     
         Mwave=(2.3*1E-4) * WP + (2.2*1E-3) * mp_wave_mort_factor

      ELSEIF (mp_break_type.EQ.break_Rodrig) THEN

C Death rate scales with Hw (Hw, surface wave height)

         Mwave=mp_spp_death * Hw * mp_wave_mort_factor 

      ENDIF

C M_wave + general Mortality rate; [d-1]       
      mortMAG = mp_spp_death + Mwave
      dNs = UptakeN * tstep_mag * Bmag		! [mg-N/ m2]
      UptakeN = dNs/tstep_mag/Bmag		! [mg-N/g-dry/day] - redefining uptake

C Compute exudation

      dE = mp_spp_E * (Q-mp_spp_Qmin) * tstep_mag

C Growth via solver

      IF (vQ .LE. c0) THEN
              
        Upc1 = c0
        
      ELSE 
              
        Upc1 = UptakeN/vQ
        
      ENDIF 

      Upc2 = 50.0/(mp_spp_Qmax-mp_spp_Qmin)

C Compute nitrogen in the macroalgae = Biomass * Quota

      Nmag = Q*Bmag

C Intermediate step to compute the tendencies : dBmagdt and dNmagdt 
C (1/(3600*24)) to convert from per day into per second

      dBmagdt = Bmag * (muMAG - mortMAG) * (1/(3600.*24.)) 

      dNmagdt= (Bmag*Upc1* (1.0- 1.0/(1.0+ (mp_spp_Qmax-Nmag/Bmag)* 
     & Upc2))-(Nmag*mortMAG)-(Nmag - mp_spp_Qmin*Bmag)*mp_spp_E)*
     & (1/(3600.*24.))   

      gtr(iMagB)=dBmagdt
      gtr(iMagN)=dNmagdt
      gtr(iMagC)=(dBmagdt * (30./100.))

C Extreme carbon limitations 
C If expected C<Cmin then increase C, 
C so that C=Cmin & compute the loss of C and N from structural mass)
C IMPLEMENT THIS

C source terms for ocean nutrients [their units/s ]
C IMPLEMENT THIS
             
      RETURN
      END SUBROUTINE
      
