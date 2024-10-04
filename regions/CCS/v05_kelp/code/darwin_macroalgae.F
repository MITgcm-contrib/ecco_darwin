#include "DARWIN_OPTIONS.h"

CBOP

        SUBROUTINE DARWIN_MACROALGAE(
C     I  t,hor_vel,parl,ptr,lat,dayyr,
     I  t,hor_vel,parl,ptr,
     I  Hw,Tw,
     U gtr,
C     O Aalg,Nres,Cres,
     O gT,gQ,gH,gE,muMAG,
     O mortMAG)
CI ,dBmagdt,dNmagdt,



C   t:         temperature  [C]
C   hor_vel:   horizontal velocity [m s-1]
C   parl: photosynthetic radiation [mu Ein s-1 m-2]
C   ptr(ndarwin): tracer concentrations 
C   lat: latitude [deg]
C   dayyr: day in a year [0-365/366]
 
C   gtr(ndarwin): tracer tendencies (inout)

C
C macroalgage relevant tracers in ptr
C ptr(MaNum) .... number of macroalgage [km^-3]
C ptr(MaFrA) .... concentration of macroalgage structural (frond) area [dm^2/km^-3]
C ptr(MaC) .... concentration of C in macroalgage (storage) [uM]
C ptr(MaN) .... concentration of N in macroalgage (storage) [uM]
C
C gtr = tendency dptr/dt due to physical processes [ptr s-1]
C

C Based on Broch et al and
C Strong-Wright paper https://github.com/jagoosw/SugarKelp.jl


      IMPLICIT NONE

#include "DARWIN_SIZE.h"
#include "DARWIN_INDICES.h"
#include "DARWIN_DIAGS.h"


C in outs

      _RL ptr(nDarwin)
      _RL gTr(nDarwin)
      _RL t,hor_vel,parl


C      _RL Aalg 
C      _RL NRes
C      _RL CRes 
C      _RL dayyr
C      _RL lat

C   internal variables
C      _RL NumAlg
C      _RL WgtStr
C      _RL f_area
C      _RL f_photo
C      _RL f_temp
C      _RL mu
C      _RL nu
C      _RL s_n
C      _RL exp_f
C      _RL ps
C      _RL p
C      _RL e
C      _RL r
C      _RL Cupt
C      _RL Nupt
C      _RL Nenv

C     temp. variable  of NO3
      _RL N_NO3

C      _RL beta
C      _RL J
C      _RL lambda
C      _RL tk
C      _RL pmax

C tendencies
C      _RL t_FrA
C      _RL t_N
C      _RL t_C

C      _RL LENGTH_DAY
      
C constants (names start with c_)
C      _RL c_ka
C      _RL c_m1
C      _RL c_m2
C      _RL c_A0
C      _RL c_eps
C      _RL c_Jmax
C      _RL c_U065
C      _RL c_Nmax
C      _RL c_Nmin
C      _RL c_J
c      _RL c_Nupt
C      _RL c_Isat
C      _RL c_alph
C      _RL c_Cmin
c      _RL c_Tar
C      _RL c_Tr1
C      _RL c_Tr2
C      _RL c_r1
C      _RL c_gamm
C      _RL c_kx
C      _RL c_Cmol
C      _RL c_Nmol
C      _RL c_Cstr
C      _RL c_Nstr
C      _RL c_a1
C      _RL c_a2
C      _RL c_Tap      
C      _RL c_Tp1
C      _RL c_Tapl
C      _RL c_Tpl
CC      _RL c_Taph
C      _RL c_Tph
C      _RL c_P1
      
C-MM ADDING PARAMETERS OF MACMODS

C PHYSIOLOGICAL PARAMETERS FOR
C MACROCYCSTIS PYRIFERA - GIANT KELP

      _RL mp_spp_Vmax       
      _RL mp_spp_Ks_NO3       
      _RL mp_spp_kcap         
      _RL mp_spp_Gmax_cap     
      _RL mp_spp_PARs      
      _RL mp_spp_PARc     
      _RL mp_spp_Q0           
      _RL mp_spp_Qmin          
      _RL mp_spp_Qmax          
      _RL mp_spp_BtoSA          
      _RL mp_spp_line_sep      
      _RL mp_spp_kcap_rate    
      _RL mp_spp_Topt1         
      _RL mp_spp_K1            
      _RL mp_spp_Topt2         
      _RL mp_spp_K2             
      _RL mp_spp_CD             
      _RL mp_spp_dry_sa        
      _RL mp_spp_dry_wet     
      _RL mp_spp_E           
      _RL mp_spp_seed         
      _RL mp_spp_death        

C-MM End of physiological parameters
C declaration


C Adding some variables needed for
C physiological calculation      
C in MACMODS
        
      _RL Qmin
      _RL Qmax
      _RL Topt1
      _RL Topt2
      _RL K1
      _RL K2
      _RL CD           ! Drag coffcient
      _RL VmaxNO3
      _RL PARc
      _RL PARs 
      _RL vQ
   
C Now variable relative to water flow
C and wave motion affecting nutrient uptake

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
      _RL KsNO3
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

C More MACMODS variables
      _RL tstep_mag
      _RL dE ! Exudation rate
C      _RL vQ
      _RL Upc1,Upc2
C Tendency terms
      _RL dBmagdt
      _RL dNmagdt
      _RL dCmagdt
            


C Key variables to store forcing functions
C from ECCO-Darwin
C mag = macroalgae

      _RL t_mag ! water temperature in degC
      _RL par_mag ! PAR in Watts/m2      
      _RL NO3_mag ! NO3 conc. in mmol/m3
      _RL chl_plk ! chlrophyll from plankton (test variable)

C  Factors affecting growth          
C Growth, nitrogen movement from Ns to Nf = umax*gQ*gT*gE*gH; [per day]
C Output:
C!   Growth, [h-1]
C   gQ, quota-limited growth
C!       from Wheeler and North 1980 Fig. 2
C   gT, temperature-limited growth
C       piecewise approach taken from Broch and Slagstad 2012 (for sugar
C       kelp) and optimized for Macrocystis pyrifera
C!   gE, light-limited growth
C!       from Dean and Jacobsen 1984
C   gH, carrying capacity-limited growth
C   muMAG, macroalgal growth      


      _RL gQ,gT,gE,gH,muMAG,mortMAG,WP
      _RL Mwave, mwh


      _RL par_watts,chlmin,mp_farm_depth
      _RL A,B,Bnew
      _RL dNs
      

C more parameters to be set later down for  

      _RL mp_growth_lim_type 
      _RL break_Duarte_Ferr
      _RL break_Rodrig
      _RL mp_break_type ! Frond Breakage Type
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

C Some more parameters for MACMODS

       _RL mw_n       ! molecular weight N
       _RL sphr           ! number of seconds in an hour
       _RL spd           ! number of seconds in a day
       _RL dpy             ! number of days in a year
       _RL spy           ! number of seconds in a year
       _RL hrps        ! number of hours in a second
       _RL dps          ! number of days in a second
       _RL ypd          ! number of years in a day
       _RL yps          ! number of years in a second
       _RL cmperm          ! cm per meter
       _RL mpercm           ! meters per cm      
C
C-MM Physical constants used in MACMODS

       _RL vonkar                ! von Karman constant
       _RL T0_Kelvin             ! freezing T of fresh water (K)
       _RL K_Boltz      ! Boltzmann constant (eV/K)
       _RL rho_sw              ! density of salt water (g/cm^3)
       _RL epsC              ! small C concentration (mmol C/m^3)
       _RL epsTinv          ! small inverse time scale (1/year) (1/sec)
       _RL molw_Fe             ! molecular weight of iron (gFe / mol Fe)
       _RL R13C_std               ! actual 13C/12C PDB standard ratio (Craig, 1957) = 1123.72e-5_r8
       _RL R14C_std              ! actual 14C/12C NOSAMS standard ratio = 11.76e-13_r8


            
CC This from K. Suselj
C numerical parameters (names start with n_)      
C      _RL n_NumAlg
C      _RL n_AAlg

C temp
C      _RL test_day
C      _RL test_lat
C      _RL test_dl


C time constant
	_RL day_in_s
	_RL h_in_s
	
C number of seconds in day	
	  parameter(day_in_s=86400.)
C number of s in hour
      parameter(h_in_s=3600.)

CC-MM PARAMERTERS OF MACMODS WITH 
C ASSIGNED VALUES
C for the Giant Kelp (Macrocystis pyrifera)

       parameter(mp_spp_Vmax=12.8*24)         ! [umol N/g-DW/d], expecting umol N/m2/d, conversion below
       parameter(mp_spp_Ks_NO3= 10130)        ! [umol N/m3]  %
       parameter(mp_spp_kcap= 2000)           ! [g(dry)/m]
       parameter(mp_spp_Gmax_cap= 0.2)        ! [1/day]
       parameter(mp_spp_PARs= 212.4/4.57)     ! [W/m2]
       parameter(mp_spp_PARc= 20.45/4.57)     ! [W/m2]
       parameter(mp_spp_Q0= 32.0)             ! initial Q [mg N/g(dry]
       parameter(mp_spp_Qmin= 10.18)          ! [mg N/g(dry)]
       parameter(mp_spp_Qmax= 54.0)           ! [mg N/g(dry)]
       parameter(mp_spp_BtoSA= 1.0)           ! Hmm not used right now???
       parameter(mp_spp_line_sep= 1.0)        ! m
       parameter(mp_spp_kcap_rate= 0.05)      ! [1/day]
       parameter(mp_spp_Topt1= 13.0)          ! [deg C]
       parameter(mp_spp_K1= 0.04)             ! temp func slope 1
       parameter(mp_spp_Topt2= 18.0)          ! [deg C]
       parameter(mp_spp_K2= 0.05)             ! temp func slope 2
       parameter(mp_spp_CD=0.5)               ! drag coefficient (unitless)
       parameter(mp_spp_dry_sa= 58.0)         ! [g(dry)/m2]
       parameter(mp_spp_dry_wet= 0.094)       ! [g(dry)/g(wet)] % Not changed from the macrocystis values
       parameter(mp_spp_E= 0.01)              ! [d-1] % No info specific for Eucheuma
       parameter(mp_spp_seed= 50.0)           ! initial biomass [g(dry)/m]
       parameter(mp_spp_death= 0.01)          ! death rate [1/day]

CC - other MACMODS parameters with assigned values 

      parameter(break_Duarte_Ferr = 0)    ! fancy breakage
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

C-MM more MACMODS parameters

      parameter(mw_n      = 14.00672)   ! molecular weight N
      parameter(sphr      = 3600.0)     ! number of seconds in an hour
      parameter(spd       = 86400.0)    ! number of seconds in a day
      parameter(dpy       = 365.0)      ! number of days in a year
      parameter(spy       = dpy*spd)    ! number of seconds in a year
      parameter(hrps      = c1 / sphr)  ! number of hours in a second
      parameter(dps       = c1 / spd)   ! number of days in a second
      parameter(ypd       = c1 / dpy)   ! number of years in a day
      parameter(yps       = c1 / spy)   ! number of years in a second
      parameter(cmperm    = 100.0)      ! cm per meter
      parameter(mpercm    = 0.01)        ! meters per cm

C-MM Some more physical constants for MACMODS code


      parameter(vonkar    =   0.4)            ! ! von Karman constant
      parameter(T0_Kelvin = 273.15)           ! ! freezing T of fresh water (K)
      parameter( K_Boltz   =   8.617330350*1E-5)  ! ! Boltzmann constant (eV/K)
      parameter(rho_sw    =   1.026)           ! density of salt water (g/cm^3)
      parameter(epsC      =   1.0*1E-8)         ! ! small C concentration (mmol C/m^3)
      parameter(epsTinv   =   3.17*1E-8)        ! ! small inverse time scale (1/year) (1/sec)
      parameter(molw_Fe   =  55.845)          ! molecular weight of iron (gFe / mol Fe)
      parameter(R13C_std  =   1.0)            ! actual 13C/12C PDB standard ratio (Craig, 1957) = 1123.72e-5_r8
      parameter(R14C_std =    1.0)            ! actual 14C/12C NOSAMS standard ratio = 11.76e-13_

      parameter(tstep_mag =    1200)

C       write(*,*)'WAVE_HEIGHT = ',Hw
C       write(*,*)'WAVE_PERIOD = ',Tw

C End of all the variables and parameters
C needed by MACMODS code
CC-MM
      
CC- MM This code from Kay Suselj      
C !! NOW COMMENTED OUT
C  structural dry weight per unit area [g  SW dm-2] - wgt=c_ka*A
C      parameter(c_ka=0.6)     
C  minimal N reserve [(g N)/(g SW)]    
C      parameter(c_Nmin=0.01)
C  minimal C reserve  [(g C)/(g SW)]   
C      parameter(c_Cmin=0.01)
C  maximum N reserve [(g C)/(g SW)]   
C      parameter(c_Nmax=0.022)
C growth-rate adjustment parameter [-]
C      parameter(c_m1=0.1085)
C growth-rate adjustment parameter [-]
C      parameter(c_m2=0.03)
C growth-rate adjustment parameter [dm2]
C      parameter(c_A0=6.0)
C front erosion parameter [dm^-2]
C      parameter(c_eps=0.22)
C maximal nitrate uptake rate [g N dm-2 h-1]
C      parameter(c_Jmax=1.4e-4)
C current speed at which J = 0.65 Jmax [m s-1]
C      parameter(c_U065=0.03)
C irradiance for maximal photosynthesis [mu mol photons m-2 s-1 = mu E m-2 s-1]
C      parameter(c_Isat=200.0)
C photosynthetic efficiency [(g C dm-2 h-1) /(mu mol photons m-2 s-1)]
C      parameter(c_alph=3.75e-5)
C arrhenious temperature for respiration [K]
C      parameter(c_Tar=11.0)
C reference temperature for respiration [K]
C      parameter(c_Tr1=285.0)
C reference temperature for respiration [K]
C      parameter(c_Tr2=290.0)
C respiration rate at T=Tr [g C dm-2 h-1]
C      parameter(c_r1=2.785e-4)
C exudation parameter [gC/g]
C      parameter(c_gamm=0.5)
C nitrate uptake half-saturation constant [mu mol L-1]
C      parameter(c_kx=4.0)
C molar mass of C [g mol-1]
C      parameter(c_Cmol=12.0)      
C molar mass of N [g mol-1]
C      parameter(c_Nmol=14.0)      
C  mass of structural C per mass of structural mass [gC /(g sw)]
C      parameter(c_Cstr=0.20)     
C  mass of structural N per mass of structural mass [gC /(g sw)]
C      parameter(c_Nstr=0.01) 
C photoperiod parameter [-]
C      parameter(c_a1=0.85) 
C photoperiod parameter [-]
C      parameter(c_a2=0.3) 
C photosynthetic temperatures, all in [K]      
C       Tap=1./(1./c_Tp1-1./c_Tp2)*log(c_p2/c_p1)
C      parameter(c_Tap=1649.4)      
C      parameter(c_Tp1=285.)
C      parameter(c_Tapl=27774.)
C      parameter(c_Tpl=271.)
C      parameter(c_Taph=25924.)
C      parameter(c_Tph=296.)
C max photosynthetic rate at T=Tp1 [g C dm-2 h-1]
C      parameter(c_p1=1.22e-3)

C minimum amount of algae [#/km^3]
C       parameter(n_NumAlg=1.e-12)
C minimum frond area [dm^2]      
C       parameter(n_AAlg=1.e-7)

C  STARTING HERE MACMODS CODE
C  COMPUTATION 

C ASSINGING VALUE OF FROND BREAKAGE TYPE          

       mp_break_type = break_Duarte_Ferr

CC-MM
C Now code from MACMODS starts here

      NO3_mag=ptr(iNO3) ! mmol/m3

C Here exclude any value that was less than 0
      IF (NO3_mag .LT. 0 ) THEN 

         NO3_mag=0.0      

      ENDIF


      t_mag=t ! deg C
      PAR_mag=parl ! Watts/m2

C Lines used for debugging
C      NO3_mag=0.1458  ! mmol/m3
C      t_mag= 24       ! deg C
C      par_mag= 259.2042 ! Watts/m2


C Now fake value of chl from plakton for test 
      chl_plk=0.1 ! mg/m3      

C growth Limitation model
C if = 0 ==>  multiply all the limitaton factors
C if > 0 ==>  independent light and nutrient limitation growth model

      mp_growth_lim_type=0


C Assigning Q parameters

        Qmin=mp_spp_Qmin           ! [mg N/g(dry)]
        Qmax=mp_spp_Qmax


C Assiging thermal parameters

        Topt1=mp_spp_Topt1
        Topt2=mp_spp_Topt2

        K1=mp_spp_K1
        K2=mp_spp_K2


C Assigning drag coeff.
         CD=mp_spp_CD

         VmaxNO3=mp_spp_Vmax
         KsNO3=mp_spp_Ks_NO3

C Setting the kind of nuts uptake

          mp_N_uptake_type=1

C Assigning drag coeff  light parameters

          PARc = mp_spp_PARc;
          PARs = mp_spp_PARs;


C LAMBDA_NO3 FUNCTION

         
C         write(*,*),'Wave Height = ' WvHeight(i,j,k)


C         Tw = 10; ! Wave Period in seconds (typical number)
         n_length = 1 
         visc = 1E-6* 86400
         Dm = (18*3.65*1E-11 + 9.72*1E-10) * 86400

         magu = hor_vel ! m/s (from ECCO CCS model)

C         ! unit conversions, minima
         magu_m = max(1.0,magu*86400) ! m/day
         Tw_s = max(0.01,Tw)/86400.0  ! whaa...?

         DBL = 10 * (visc / (sqrt(CD) * abs(magu_m)));

C for now set it to 1
         n_length=1
C  Oscillatory Flow

        DO n = 1,n_length
        vval = (1-exp((-Dm * n**2 * pi**2 *Tw_s)/(2.0*DBL**2)))/(n**2 *
     &   pi**2)
        ENDDO
C        oscillatory = ((4.0*DBL)/Tw_s) * sum(vval)  ! original line
         oscillatory = ((4.0*DBL)/Tw_s) * vval

C       2. Uni-directional Flow
        flow = Dm / DBL

       beta = flow + oscillatory

       lambda_NO3 = c1 + (VmaxNO3 / (beta*KsNO3)) - (NO3_u/KsNO3)
       lambda=lambda_NO3

C   end of lambda

C Part on growth depending on temperature       
        
        IF (t_mag .GE. Topt1) THEN 
           IF (t_mag .LE. Topt2) THEN
            temp_lim = 1.0
           ELSE 
            temp_lim = exp(-K2*(t_mag-Topt2)**2)
           ENDIF 
        ELSE
        temp_lim = exp(-K1*(t_mag-Topt1)**2)
       ENDIF 
      temp_lim = max(0.0,temp_lim)
      temp_lim = min(1.0,temp_lim)
      gT=temp_lim

C end part on temp limitation 

C Computing now 

C      Bmag= mp_spp_seed
      Q = Qmin + NO3_mag*(Qmax-Qmin)/35
C here converting into micromol/m3 from MM
      Q = Qmin + (NO3_mag*1E3)*(Qmax-Qmin)/35

      write(*,*)'Q with formula = ',Q

       Bmag = ptr(iMagB)
       Nmag = ptr(iMagN)
       Q = Nmag/Bmag 
       

       write(*,*)'Bmag = ', Bmag
       write(*,*)'Nmag = ', Nmag
       write(*,*)'Nmag/Bmag= ', Q

      IF (mp_N_uptake_type .EQ.0) THEN

       vQ = (Qmax-Q)/(Qmax-Qmin)
       vQ = max(c0,vQ)
       vQ = min(c1,vQ)


      ELSEIF (mp_N_uptake_type .EQ.0) THEN

      vQ = 1.0 - 1.0/(1.0+(max(0.0,Qmax-Q))*50.0/(Qmax-Qmin))
      vQ = max(c0,vQ)
      vQ = min(c1,vQ)

      ENDIF

C Below is what we call "Uptake Factor." It varies betwen 0
C and 1 and includes kinetically limited uptake and
C mass-transfer-limited uptake (oscillatory + uni-directional flow)

        NO3_u = NO3_mag*1000.0
        vNuTw_NO3 = NO3_u / (KsNO3 * ((NO3_u/KsNO3)+ 0.5
     & *(lambda+sqrt(lambda**2 + c4 * (NO3_u/KsNO3)))))
        vNuTw_NO3 = max(c0,vNuTw_NO3)
        vNuTw_NO3 = min(c1,vNuTw_NO3)   

C Uptake Rate [mg N/g(dry)/d]
C Nutrient Uptake Rate = Max Uptake * v[Ci,u,Tw] * vQ
C converted from umol N/m2/d -> mg N/g(dry)/d by 14.0067 / 1e3

         Uptake_NO3 = VmaxNO3 * vNuTw_NO3 * vQ ! [umol/m2/d]
         Uptake_NO3 = Uptake_NO3 * mw_n / 1E3
         UptakeN = Uptake_NO3
         UptakeN = UptakeN / mp_spp_dry_sa

C Growth, nitrogen movement from Ns to Nf = umax*gQ*gT*gE*gH; [per day]
C Output:
C!   Growth, [h-1]
C   gQ, quota-limited growth
C!       from Wheeler and North 1980 Fig. 2
C   gT, temperature-limited growth
C       piecewise approach taken from Broch and Slagstad 2012 (for sugar
C       kelp) and optimized for Macrocystis pyrifera
C!   gE, light-limited growth
C!       from Dean and Jacobsen 1984
C   gH, carrying capacity-limited growth

CgQ = (Q(i,j) - params(mp_spp_Qmin)) / Q(i,j) % Droop equation
        gQ = (Q - Qmin) / Q * Qmax/(Qmax-Qmin) ! Droop scaled from 0-1
CgQ = (Q(i,j) - Qmin) / (Qmax-Qmin) ! Freider et al.

        write(*,*)'gQ = ',gQ
              
        gQ = max(c0,gQ)
        gQ = min(c1,gQ)

        par_watts = par_mag

        IF (par_watts .LT. PARc) THEN
          gE = c0
        ELSEIF (par_watts .GT. PARs) THEN
          gE = c1
         ELSE
        gE = (par_watts-PARc)/(PARs-PARc)*exp(-(par_watts-PARc)
     &  /(PARs-PARc)+c1)
         ENDIF

C         write(*,*)'gE (dar_mag) = ',gE         

C         consider daylength if timestep is > 1/2 a day
C        if (params(mp_dt_mag) > 0.51) then
C            day_h = daylength(lat(i,j),lon(i,j),c0,params(mp_dte)) ! daylength in h
C            gE = min(day_h,gE*day_h)
C        endif

C         ! Carrying capacity
C        ! ----------------------------------------------------------
C          ! gH -> density-limited growth (ranges from the max growth rate to the death rate, excluding wave mortality)
C          ! This expression follows Xiao et al (2019 and ignores wave mortality when
C        ! thinking about the death rate

C        Bmag = ptr(iMagB)

        Bnew = Bmag *mp_spp_line_sep ! converting from g/m2 to g/m


         A = mp_spp_kcap_rate/mp_spp_kcap**(-0.75)
         gH = A*Bnew**(-0.75)

C

C Now computing the growth
C Using all the factors computed above

         IF (mp_growth_lim_type .EQ. 0) THEN
                 ! interacting limitation terms growth model
             muMAG =  min(gH,mp_spp_Gmax_cap) * gT * gE * gQ
         ELSE 
                ! independent light and nutrient limitation growth model
                ! this is the preferred formula
                ! used by MACMODS group
             muMAG =  min(gH,mp_spp_Gmax_cap) * gT * min(gE,gQ)
         ENDIF

        
C             write(*,*)'temp = ',t_mag 
C             write(*,*)'PAR = ', par_watts
C             write(*,*),'NO3 = ',NO3_mag
C             write(*,*)'gT = ',gT  
C             write(*,*)'gH = ',gH
C             write(*,*)'gQ = ',gQ
C             write(*,*)'gE = ',gE
C             write(*,*)'muMAG = ',muMAG

C Now computing the mortality rate due to wave activity      
C Wave height is computed by the ERA wave model

C !------- 
C ! d_wave = frond loss due to waves; dependent on Hs, significant
C ! wave height [m]; Rodrigues et al. 2013 demonstrates linear relationship
C ! between Hs and frond loss rate in Macrocystis [d-1] (continuous)
C ! Duarte and Ferreira (1993) find a linear relationship between wave
C power and mortality in red seaweed        

C Assign here the value of wave height

       mwh  = Hw    ! Assign wave height 
       wp_tmp = Tw  ! Assign wave period

C       write(*,*)'mwh = ',mwh


       IF (mp_break_type.EQ.break_Duarte_Ferr ) THEN

       WP = (rho_sw * (9.8**2)) /(64.0*pi) * (mwh**2) *wp_tmp / 1E3     
       Mwave=(2.3*1E-4) * WP + (2.2*1E-3) * mp_wave_mort_factor
      
C       write(*,*)'mwh = ', mwh 
C       write(*,*)'wp_tmp =', wp_tmp

C       write(*,*)'WP_DFerr = ', WP 
C       write(*,*)'Mwave_DFerr = ',Mwave
       ELSEIF (mp_break_type.EQ.break_Rodrig) THEN

C death rate scales with Hw (Hw, surface wave height)

       Mwave=mp_spp_death * Hw * mp_wave_mort_factor 

C        write(*,*)'Mwave_Rod = ',Mwave

       ENDIF

C M_wave + general Mortality rate; [d-1]       
       mortMAG = mp_spp_death + Mwave
C        mortMAG = 50



C        write(*,*),'Wave_height_mag = ', Hw 
C        write(*,*),'Wave_per_mag = ', Tw
C        write(*,*),'mp_spp_death = ',mp_spp_death
C        write(*,*),'Mwave = ', Mwave
C        write(*,*)'mortMAG = ',mortMAG



CC-MM: End of part of mortality computation

       dNs = UptakeN * tstep_mag * Bmag   ! [mg-N/ m2]


       UptakeN = dNs/tstep_mag/Bmag ! [mg-N/g-dry/day] - Redefining Uptake

C Now compute exudation

       dE = mp_spp_E * (Q-Qmin) * tstep_mag

C now growth via solver

       IF (vQ .LE. c0) THEN
            Upc1 = c0
       ELSE 
            Upc1 = UptakeN/vQ
       ENDIF 

        Upc2 = 50.0/(Qmax-Qmin)

C Now compute Nitrogen in the macroalgae = Biomass * Quota

        Nmag = Q*Bmag

C Now intermaediate step to compute the tendencies : dBmagdt and dNmagdt 
C (1/(3600*24)) to convert from per day into per second

        dBmagdt = Bmag * (muMAG - mortMAG) * (1/(3600.*24.)) 

        dNmagdt= (Bmag*Upc1* (1.0- 1.0/(1.0+ (Qmax-Nmag/Bmag)* 
     &      Upc2))-(Nmag*mortMAG)- 
     &    (Nmag - Qmin*Bmag)*mp_spp_E)* (1/(3600.*24.))   


C        gtr(iMagB)=gtr(iMagB)+dBmagdt
C        gtr(iMagN)=gtr(iMagN)+ dNmagdt
C        gtr(iMagC)=gtr(iMagC)+ (dBmagdt * (30/100))

         gtr(iMagB)=dBmagdt
         gtr(iMagN)=dNmagdt
         gtr(iMagC)=(dBmagdt * (30./100.))

C           gtr(iMagB)=10
C           gtr(iMagN)=20
C           gtr(iMagC)=30


C        IF (muMAG.GT.0) THEN

C         write(*,*)'iMagB = ', iMagB 
C         write(*,*)'iMagN = ', iMagN
C         write(*,*)'iMagC = ', iMagC


C        write(*,*)'muMAG = ', muMAG
C        write(*,*)'dBmagdt = ', dBmagdt
C        write(*,*)'dNmagdt = ', dNmagdt
C        write(*,*)'Bmag = ',Bmag

C        IF (Upc1.GT.0) THEN
C        write(*,*)'Upc1 > 0 ', Upc1
C        ENDIF
        

C        write(*,*)'Upc2 = ', Upc2
C        write(*,*)'Q = ',Q
C        write(*,*)'Qmin = ',Qmin
C        write(*,*)'Qmax = ',Qmax
C        write(*,*)'mortMAG =', mortMAG
C        write(*,*),'mp_spp_E = ',mp_spp_E
C        write(*,*)'tstep_mag = ', tstep_mag

C        write(*,*)'gtr_MagB = ', gtr(iMagB)
C        write(*,*)'gtr_MagN = ', gtr(iMagN)
C        write(*,*)'gtr_MagC = ', gtr(iMagC)

C        ENDIF


C  30 % conversion factor from biomass to carbon 


C             d_Bmag = dBmagdt
C             Nitrogen of MAG
C             Q_new = (Nmag+dNmagdt) / (Bmag+dBmagdt);
C             d_Q = Q_new - Q;
C             Q = Q_new;


          
C End of MACMODS code
CC-MM


C
C   values for numerical stability (n_)
C
C minimum amount of algage [#/km^3]
C       parameter(n_NumAlg=1.e-12)
C minimum frond area [dm^2]      
C       parameter(n_AAlg=1.e-7)


C
C       print*,'XXXXXXXXXXXXXXXXXXX'
C       print*,'BEGINNING COMP PARAMETERS'


C     temperature in K
C       tk=t+273.16 ;


C       print*,'tk =', tk
C       print*,'Temp = ',t


C number of macroalgage per km^3       
C       NumAlg=ptr(maNum)

C       print*,'n_NumAlg = ' , n_NumAlg
C       print*,'NumAlg = ',  NumAlg




C       IF (NumAlg .gt. n_NumAlg) THEN
       
C frond area of a typical/average algae [dm^2]         
C        Aalg=ptr(MaFrA)/NumAlg

C        print*,'Aalg = ',Aalg
C        print*,'n_AAlg =',n_AAlg

C        IF (Aalg .gt. n_AAlg) THEN


C            print*,'Aalg = ',Aalg
C            print*,'n_AAlg =',n_AAlg

C
C do calculation only when there are
C any algage in the grid box and the frond area is nonzero 
C (strictly speaking when the values are greater than thresholds, 
C to avoid numerical issues dividing by small values)
C


C    N in the ocean that macroalge is uptaking [mmol l-1]
C    original units are in 1uM = 1 mmol l-1
C       Nenv=(ptr(iNH4)+ptr(iNO3))

C       N_NO3 = ptr(iNO3)

C        print*,'Nenv =',  Nenv 
C        print*,'parl = ', parl
C        print*,'NO3 = (mmol/m3)', N_NO3


C dry weight of a typical algage [g]          
C          WgtStr=c_ka*Aalg
C N reserve in macroalgage [g N/g dw]       
C
C          Nres=1.e6*c_Nmol*ptr(maN)/(NumAlg*WgtStr)

C          Nres=c_Nmol*ptr(maN)/(NumAlg*WgtStr)
C C reserve in macroalgage [g C/ g dw]
C          Cres=1.e6*c_Cmol*ptr(maC)/(NumAlg*WgtStr)

C          Cres=c_Cmol*ptr(maC)/(NumAlg*WgtStr)


C         print*,'Cres = ',Cres
c         print*,'c_Cmol = ',c_Cmol
c         print*,'ptr(maC) = ', ptr(maC)
c         print*,'NumAlg = ',NumAlg
C         print*,'WgtStr = ',WgtStr

C         print*,'Nres = ',Nres
C         print*,'ptr(maN) = ',ptr(maN)

C   Eq. 3
C       f_area=c_m1*exp(-(Aalg/c_A0)**2)+c_m2
C         print*,'f_area = ',f_area

C   Eq. 4     
C Orginal code i mplemented by Kaj Suselj
C From Slagstad  & Broch
C for Saccharina Latissima
C       IF (T .lt. -1.8) THEN
C       f_temp=0.0
C       ELSEIF (T .lt. 10.) THEN
C       f_temp=0.08*T + 0.2 
C       ELSEIF (T .le. 15.) THEN
C       f_temp=1.0
C       ELSEIF (T .le. 19.) THEN
C       f_temp=4.75-T/4.0
C       ELSE
C       f_temp=0.0
C       ENDIF

C Now implementing here the  
C code from MACMODS used the 
C the Giant Kelp of the Pacific Ocean

C       IF (T .lt. 14) THEN
C       f_temp = T/14

C       ELSEIF (T.ge.14 .AND. T.le.20) THEN 

C       f_temp = 1

C       ELSEIF (T.ge.20 .AND. T.lt. 23  ) THEN

C       f_temp = ((-1/3)*T)+7.7


C       ELSEIF  (T .ge. 23) THEN        
C       f_temp = 0.0        

C       ENDIF

C Safety Check if f_temp becomes 
C greater than 1
C       IF (f_temp .GT. 1) THEN
C       f_temp = 1.0         
C       ENDIF         

C       print*,'f_temp = ',f_temp
C       print*,'XXX' 
       
C Eq. 5    
C       call LAMBDA_DAY(dayyr,lat,lambda)
        
C       IF (lambda .ge. 0.) THEN
C       	f_photo=c_a1*(1.+sqrt(lambda))+c_a2
C       ELSE
C        f_photo=c_a1*(1.-sqrt(-lambda))+c_a2
C       ENDIF
      
       
C growth function  - mu, Eq. 2 in B&S   [day-1]  
C - it seems like we just pick-up units of 1/day -       
C        mu=f_area*f_photo*f_temp*
C     &  min(1.0-c_Nmin/NRes,1.0-c_Cmin/CRes)
       
       
C frond loss - nu, Eq. 6 [day-1]
C       exp_f=exp(c_eps*AAlg)   
C       nu=(1.0e-6*exp_f)/(1.0+1.0e-6*(exp_f-1.0)) 
       

C  Eq. 8   
C      J=c_Jmax*(1.0-exp(-hor_vel/c_U065))
C     & *(c_Nmax-Nres)/(c_Nmax-c_Nmin)
C     & * Nenv/(c_Kx+Nenv)  
      
C  N uptake, Eq 7 [gN/gDW h-1]     
C      Nupt=J/c_ka 

C     Eq. 13
C      pmax=c_p1*exp(c_Tap/c_Tp1-c_Tap/tk)/
C     & (1.+exp(c_Tapl/tk-c_Tapl/c_Tpl)+exp(c_Taph/c_Tph-c_Taph/tk))
C first approx. for beta
C      beta=1.e-9
C iteratively solve for beta       
C      call SOLVE_BETA(c_alph,pmax,c_Isat,beta)

C     Eq 11     
C      ps = c_alph*c_Isat/(log(1.0+c_alph/beta))
C     Eq 10
C      p=ps*(1.0-exp(-c_alph*parl/ps))*exp(beta*parl/ps)  
C     Eq. 15 
C      e =1.0-exp(c_gamm*(c_Cmin-Cres))
C     Eq. 14
C      R = c_r1 * exp(c_Tar / c_Tr1 - c_Tar / (t+273.16))

C C uptake - Eq. 9     [gC/gDW h-1]

C   	  Cupt=(p*(1.0-e)-R)/c_ka

C          print*,'Cupt = ',Cupt
C          print*,'Cres = ',Cres
C          print*,'c_Cmin = ',c_Cmin
C          print*,'c_gamm = ',c_gamm
C          print*,'p = ',p
C          print*,'e = ',e
C          print*,'R = ',R
C          print*,'c_ka = ',c_ka          
 
C   
C extreme carbon limitations 
C (if expected C<Cmin then increase C, 
C so that C=Cmin & compute the loss of C and N from structural mass)
C IMPLEMENT THIS

C

C
C source terms for macroalgae variables  [their units/s]
C

C The number of macroalgage stays constant
C d/dt(ptr(MaNum)) = 0 
C don't need to do anything


C  ptr(MaFrA)=Aalg*NumAlg
C  d/dt ptr(MaFrA)=NumAlg*d/dt(A)=NumAlg*Aalg*(mu-nu)/(86400)
C        mu,nu [day-1], (mu,nu)/86400 [s-1] 
C        NumAlg [# macro/km^3]
C        Aalg  [dm^2]
C  t_Fra in dm^2 km^-3 s^-1
      
C      t_FrA=NumAlg*Aalg*(mu-nu)/day_in_s

C      print*,'t_FrA =', t_FrA

      
C d/dt ptr(MaN)=(NumAlg*c_ka)/(1.e6 c_Nmol) d/dt(Aalg*N_res)
C d/dt(Aalg*N_res)=Aalg(N_upt-N_res*nu-mu*Nstr)
C mu,nu = day-1
c Nupt = h-1

C      t_N=1.0e-6*NumAlg*c_ka/c_Nmol*Aalg*
C     & (Nupt/h_in_s-Nres*nu/day_in_s
C     & -mu*c_Nstr/day_in_s) 

C       print*,'t_N = ',t_N
C       print*,' c_Nmol =',c_Nmol
C       print*,' Nupt=',Nupt
C       print*,'Nres =', Nres



C d/dt ptr(MaN)
C same as for C
C       t_C=1.0e-6*NumAlg*c_ka/c_Cmol*Aalg*
C     & (Cupt/h_in_s-Cres*nu/day_in_s
C     & -mu*c_Cstr/day_in_s) 

C        print*,'t_C = ',t_C
C        print*,' c_Cmol =',c_Cmol
C        print*,'Cupt =', Cupt
C        print*,'Cres =', Cres

C        print*,'XXXXXXXXXXXXXX'
C        print*,'END COMPUTATION PARAMS'

C
C apply tendencies to the macroalgage tracers
C 
C      gtr(MaFrA)=gtr(MaFrA)+t_FrA
C      gtr(MaN)=gtr(MaN)+t_N
C      gtr(MaC)=gtr(MaC)+t_C 
 

C
C source terms for ocean nutrients [their units/s ]
C IMPLEMENT THIS



C
C debug
C     
     
C      print *,'inputs:'
C      print *,'tendency, macroalgage mass',gTR(maC),t,hor_vel
      
C      print *,'macroalge properties:'
C      print *,'inputs, maNum,MaFRA,MaC,MaN',
C     & ptr(maNum),ptr(MaFrA),ptr(MaC),ptr(MaN) 
      
C      print *,'other inputs',t,hor_vel,parl,lat,dayyr
      
C
C derived quantities
C      
      
C      print *,'num_alg [km-3],frond A [dm2], C and N reserves:',
C     &  NumAlg,AAlg,Cres,Nres
 
C
C growth rates
C 
C       print *,'maroalgage growth parameters: mu,nu [day-1]:',
C     & mu,nu
     
      
C      print *,'days to dep/doubl frondA, Nres and Cres:',
C     & ptr(MaFrA)/(t_FrA*86400.),
C     & ptr(MaN)/(t_N*86400.),ptr(MaC)/(t_C*86400.)
           
C       ENDIF
C       ENDIF
             
      RETURN
      END SUBROUTINE
      
C-MM
C - !!! CODE BELOW ALL COMMENTED OUT
C - NOT IN USE WITH MACMODS

C      SUBROUTINE SOLVE_BETA(a,pmax,Isat,beta)
C solve implicit equation for beta - Eqs. 12 & 13
C
C find zeros using Newton-Rapson's method
C x_{n+1}=x_n+f(x_n)/f'(x_n)
C where:
C f= Eq. 12 [r.h.s]- Eq. 13 [r.h.s.] 
      
C C     IMPLICIT NONE
      
C      _RL a
C      _RL pmax
C      _RL Isat
C in/out: in-first guess, out-final guess      
C      _RL beta
      
C internals      
C      integer max_iter,k
C      _RL ax,axx,lax,laxx,a2
C      _RL f1,f2,f3
C      _RL x
      
C      max_iter=10
      
C      x=beta
      
C      DO k=1,max_iter

C       ax=a+x
C       axx=(a+x)/x     
C       lax=log(a+x)
C       laxx=log((a+x)/x)       
C       a2=a**2      
                  
C       f1=(x**(1.-x/a))*ax*laxx
C       f2=a2*isat*(x**(x/a))-pmax*(ax**(ax/a))*laxx
C       f3=a*isat*(a2+x*ax*(log(x)-lax)*laxx)
C
C f1*f2/f3=f(x_n)/f'(x_n)
C (help with wolframalpha.com) 
C
C       IF (abs(f3) .gt. 1.e-12) THEN           
C       	x=x-f1*f2/f3
C       ENDIF
       
       
C      ENDDO
      
C       beta=x      
    
C  	  RETURN
C	  END    
      
C      SUBROUTINE LAMBDA_DAY(dayyr,lat,lambda) 
C
C using Broch model for length of day in which the longest 
C days are at day 173 (N hemisphere) and at 365 (S hemisphere)
C and days are equally long at equator	
C

C      IMPLICIT NONE
C in-s
C      _RL lat
C      _RL dayyr 

C out
C      _RL lambda
C internals
C      _RL day_l
C      _RL day_lm
C      _RL day_max
C      _RL longest_day

C length of the longest day
C (note that at 0. all days have the same length)
C	  IF (LAT .gt. 0.) THEN
C	  longest_day=173.
C	  ELSE 
C	  longest_day=365.
C	  ENDIF

C length of day    
C      call length_day(dayyr,lat,day_l)  
C length of day -1      
C      call length_day(dayyr-1.,lat,day_lm)
C length of the longest day       
C      call length_day(longest_day,lat,day_max) 
           
C	  lambda=(day_l-day_lm)/day_max
C	  
C	  RETURN
C	  END
C	     
C      SUBROUTINE LENGTH_DAY(dayyr,lat,lenday)
C Broch 1981 model for length of day 
C see Forthsyle et al. 1995 https://www.sciencedirect.com/science/article/abs/pii/030438009400034F            
C Eqs. 4-6         
       
C       IMPLICIT NONE
       
C       _RL lenday
C       _RL dayyr
C       _RL lat
C       _RL phi
C       _RL fac
         
C 360/365  =   0.98630 ...   
C       phi=23.45 * sind(0.98630136986 *(283.+dayyr)) 
C       fac=tand(phi)*tand(lat)
C       fac=max(min(fac,1.),-1.)
C need to convert radians to degrees
C 2/15=0.133 
C 180/pi = 57.2957      
C       lenday=0.1333333* 57.2957795131*acos(-fac)    
C      RETURN 
C      END      
