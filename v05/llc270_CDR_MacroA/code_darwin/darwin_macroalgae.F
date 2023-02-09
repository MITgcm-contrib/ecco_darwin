#include "DARWIN_OPTIONS.h"

CBOP


        SUBROUTINE DARWIN_MACROALGAE(
     &  t,hor_vel,parl,ptr,gtr)

C   t:         temperature  [C]
C   hor_vel:   horizontal velocity [m s-1]
C   parl: photosynthetic radiation [mu Ein s-1 m-2]
C   ptr(ndarwin): tracer concentrations  
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

      _RL ptr(nDarwin)
      _RL gTr(nDarwin)
      _RL t,hor_vel,parl

C   internal variables
	  _RL NumAlg
	  _RL Aalg 
	  _RL WgtStr
      _RL NRes
      _RL CRes   		
      _RL f_area
      _RL f_photo
      _RL f_temp
      _RL mu
      _RL nu
      _RL s_n
      _RL exp_f
      _RL ps
      _RL p
      _RL e
      _RL r
      _RL Cupt
      _RL Nupt
      _RL Nenv
      _RL beta
      _RL J

C tendencies
      _RL t_FrA
      _RL t_N
      _RL t_C

      
C constants (names start with c_)
      _RL c_ka
      _RL c_m1
      _RL c_m2
      _RL c_A0
      _RL c_eps
      _RL c_Jmax
      _RL c_U065
      _RL c_Nmax
      _RL c_Nmin
      _RL c_J
      _RL c_Nupt
      _RL c_Isat
      _RL c_alph
      _RL c_Cmin
      _RL c_Tar
      _RL c_Tr1
      _RL c_Tr2
      _RL c_r1
      _RL c_gamm
      _RL c_kx
      _RL c_Cmol
      _RL c_Nmol
      _RL c_Cstr
      _RL c_Nstr
C numerical parameters (names start with n_)      
      _RL n_NumAlg
      _RL n_AAlg


C time constant
	_RL day_in_s
	_RL h_in_s
	
C number of seconds in day	
	  parameter(day_in_s=86400.)
C number of s in hour
      parameter(h_in_s=3600.)

      
      
C  structural dry weight per unit area [g  SW dm-2] - wgt=c_ka*A
      parameter(c_ka=0.6)     
C  minimal N reserve [(g N)/(g SW)]    
      parameter(c_Nmin=0.01)
C  minimal C reserve  [(g C)/(g SW)]   
      parameter(c_Cmin=0.01)
C  maximum N reserve [(g C)/(g SW)]   
      parameter(c_Nmax=0.022)
C growth-rate adjustment parameter [-]
      parameter(c_m1=0.1085)
C growth-rate adjustment parameter [-]
      parameter(c_m2=0.03)
C growth-rate adjustment parameter [dm2]
      parameter(c_A0=6.0)
C front erosion parameter [dm^-2]
      parameter(c_eps=0.22)
C maximal nitrate uptake rate [g N dm-2 h-1]
      parameter(c_Jmax=1.4e-4)
C current speed at which J = 0.65 Jmax [m s-1]
      parameter(c_U065=0.03)
C irradiance for maximal photosynthesis [mu mol photons m-2 s-1 = mu E m-2 s-1]
      parameter(c_Isat=200.0)
C photosynthetic efficiency [(g C dm-2 h-1) /(mu mol photons m-2 s-1)]
      parameter(c_alph=3.75e-5)
C arrhenious temperature for respiration [K]
      parameter(c_Tar=11.0)
C reference temperature for respiration [K]
      parameter(c_Tr1=285.0)
C reference temperature for respiration [K]
      parameter(c_Tr2=290.0)
C respiration rate at T=Tr [g C dm-2 h-1]
      parameter(c_r1=2.785e-4)
C exudation parameter [gC/g]
      parameter(c_gamm=0.5)
C nitrate uptake half-saturation constant [mu mol L-1]
      parameter(c_kx=4.0)
C molar mass of C [g mol-1]
      parameter(c_Cmol=12.0)      
C molar mass of N [g mol-1]
      parameter(c_Nmol=14.0)      
C  mass of structural C per mass of structural mass [gC /(g sw)]
      parameter(c_Cstr=0.20)     
C  mass of structural N per mass of structural mass [gC /(g sw)]
      parameter(c_Nstr=0.01) 

C
C   values for numerical stability (n_)
C
C minimum amount of algage [#/km^3]
      parameter(n_NumAlg=1.e-12)
C minimum frond area [dm^2]      
	  parameter(n_AAlg=1.e-7)	


C number of macroalgage per km^3       
       NumAlg=ptr(maNum)
       
       IF (NumAlg .gt. n_NumAlg) THEN
       
C frond area of a typical/average algae [dm^2]         
        Aalg=ptr(MaFrA)/NumAlg

        IF (Aalg .gt. n_AAlg) THEN

C
C do calculation only when there are
C any algage in the grid box and the frond area is nonzero 
C (strictly speaking when the values are greater than thresholds, 
C to avoid numerical issues dividing by small values)
C


C    N in the ocean that macroalge is uptaking [mmol l-1]
C    original units are in 1uM = 1 mmol l-1
       Nenv=(ptr(iNH4)+ptr(iNO3))
C dry weight of a typical algage [g]          
          WgtStr=c_ka*Aalg
C N reserve in macroalgage [g N/g dw]       
          Nres=1.e6*c_Nmol*ptr(maN)/(NumAlg*WgtStr)
C C reserve in macroalgage [g C/ g dw]
         Cres=1.e6*c_Cmol*ptr(maC)/(NumAlg*WgtStr)       
       
C   Eq. 3
       f_area=c_m1*exp(-(Aalg/c_A0)**2)+c_m2

 
C   Eq. 4     
       IF (T .lt. -1.8) THEN
       f_temp=0.0
       ELSEIF (T .lt. 10.) THEN
       f_temp=0.08*T + 0.2 
       ELSEIF (T .le. 15.) THEN
       f_temp=1.0
       ELSEIF (T .le. 19.) THEN
       f_temp=4.75-T/4.0
       ELSE
       f_temp=0.0
       ENDIF



C Eq. 5 (IMPLEMENT THAT)      
       f_photo=1.
       
C growth function  - mu, Eq. 2 in B&S   [day-1]  
C - it seems like we just pick-up units of 1/day -       
        mu=f_area*f_photo*f_temp*
     &  min(1.0-c_Nmin/NRes,1.0-c_Cmin/CRes)
       
       
C frond loss - nu, Eq. 6 [day-1]
       exp_f=exp(c_eps*AAlg)   
       nu=(1.0e-6*exp_f)/(1.0+1.0e-6*(exp_f-1.0)) 
       

C  Eq. 8   
      J=c_Jmax*(1.0-exp(-hor_vel/c_U065))
     & *(c_Nmax-Nres)/(c_Nmax-c_Nmin)
     & * Nenv/(c_Kx+Nenv)  
      
C  N uptake, Eq 7 [gN/gDW h-1]     
      Nupt=J/c_ka 
       
     
C solve Eqs. 12 and 13 (this is kind of annoying, also T in K)   (FIX THAT)   
      beta=1.e-9

C     Eq 11     
      ps = c_alph*c_Isat/(log(1.0+c_alph/beta))
C     Eq 10
      p=ps*(1.0-exp(-c_alph*parl/ps))*exp(beta*parl/ps)  
C     Eq. 15 
      e =1.0-exp(c_gamm*(Cres-c_Cmin))
C     Eq. 14
      R = c_r1 * exp(c_Tar / c_Tr1 - c_Tar / (t+273.16))

C C uptake - Eq. 9     [gC/gDW h-1]
   	  Cupt=(p*(1.0-e)-R)/c_ka
   

C extreme carbon limitations (IMPLEMENT THAT)



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
      
      t_FrA=NumAlg*Aalg*(mu-nu)/day_in_s

      
C d/dt ptr(MaN)=(NumAlg*c_ka)/(1.e6 c_Nmol) d/dt(Aalg*N_res)
C d/dt(Aalg*N_res)=Aalg(N_upt-N_res*nu-mu*Nstr)
C mu,nu = day-1
c Nupt = h-1

C      t_N=1.0e-6*NumAlg*c_ka/c_Nmol*Aalg*
C     & (Nupt/h_in_s-Nres*nu/day_in_s
C     & -mu*c_Nstr/day_in_s) 



C d/dt ptr(MaN)
C same as for C
C       t_C=1.0e-6*NumAlg*c_ka/c_Cmol*Aalg*
C     & (Cupt/h_in_s-Cres*nu/day_in_s
C     & -mu*c_Cstr/day_in_s) 


C apply tendencies 
C      gtr(MaFrA)=gtr(MaFrA)+t_FrA
C      gtr(MaN)=gtr(MaN)+t_N
C      gtr(MaC)=gtr(MaC)+t_C 
      print *,'test'  






C
C source terms for ocean nutrients [their units/s ]
C





C
C debug
C     
     
      print *,'inputs:'
      print *,'tendency, macroalgage mass',gTR(maC),t,hor_vel
      
      print *,'macroalge properties:'
      print *,'inputs, maNum,MaFRA,MaC,MaN',
     & ptr(maNum),ptr(MaFrA),ptr(MaC),ptr(MaN) 
      
      print *,'num_alg [km-3],frond A [dm2], C and N reserves:',
     &  NumAlg,AAlg,Cres,Nres
 
C
C growth rates
C 
       print *,'maroalgage growth parameters: mu,nu [day-1]:',
     & mu,nu
     
      
      print *,'days to dep/doubl frondA, Nres and Cres:',
     & ptr(MaFrA)/(t_FrA*86400),
     & ptr(MaN)/(t_N*86400),ptr(MaC)/(t_C*86400)
           
       ENDIF
       ENDIF
       
       


      
      RETURN
      END SUBROUTINE
