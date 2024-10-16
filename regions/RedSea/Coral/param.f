	  module param
	
	  implicit none

	
	  REAL, PARAMETER:: dt=0.1 		!time step 
	  REAL, PARAMETER:: nyear=365.	!days of one year 
	  REAL, PARAMETER:: years=5.		!years to run the model
	  REAL :: tr						!running time
	  REAL :: tt		        		!total time in days
	  REAL :: init_day_r				
          INTEGER :: ts 				!total steps
	  INTEGER :: n,init_day
	  INTEGER,PARAMETER :: nr=2 		!number of ODE to be solved in rk4 routine
	  INTEGER :: i,trr,nrun
	 
	  REAL:: H,Hf,dH    ! Host Biomass
	  REAL:: S,Sf,dS    ! Symbiont Biomass
	  REAL:: L    	    ! Light downwelling irradiance
	  REAL:: PARdata(1825), PARdata2(365)
	  REAL:: SSTdata(1825),  SSTdata2(365), DINdata(365)
	  REAL:: L2,L1!Maximu-Minimum light for sinus
	  REAL:: X,X2,X1    ! Prey avialiability - C mol/L
	  REAL:: DIN,DIN2,DIN1  ! External Nitrogen 
          REAL:: T    ! Temperature
          REAL:: kHT  ! Temperature Dependance Function for host
          REAL:: kST  ! Temperature Dependance Function for Symbiont
	  REAL:: JHG  ! Host Biomass Formation Rate
	  REAL:: pc   ! Fixed carbon produced by phytosynthesis 
	  REAL:: JL   ! Total light absorption
	  REAL:: JCO2 ! CO2 input to photosynthesis
	  REAL:: Al   ! Ligth amplification factor
	  REAL:: JEC  ! Excess carbon used to activate host CCMs
	  REAL:: Jx   ! Functional Response Function for prey
	  REAL:: rCH  ! metabolic CO2 generated by the host availiable for photosynthesis
	  REAL:: rCs  ! CO2 generated by symbiont
	  REAL:: JHT  ! Host biomass turnover rate
	  REAL:: JSG  ! Symbiont biomass formation rate
	  REAL:: pN   ! Nitrogen shared with the symbiont
	  REAL:: JN   ! Functional Response Function for nitrogen
	  REAL:: rNH  ! Recycled nitrogen from host turnover
	  REAL:: CROS ! ROS production proportional to baseline
	  REAL:: JEL  ! Light energy in excess of photochemistry
	  REAL:: JNPQ ! Total capacity of NPQ
	  REAL:: rNS  ! Recycled nitrogen from symbiont turnover
	  REAL:: JST  ! Symbbiont Biomass Turnover
          REAL:: y(nr),dydx(nr),yout(nr)
          REAL ::JCP
          REAL :: aT,bT,cT,aT2,bT2,cT2 ! for calculation of kT
          REAL :: Length,Vol,Length_f,Mv,Volf,Vinit,Area,SGR
          REAL :: S_g,S_cells,S_cells_cm2
  	  REAL :: V,dV,SA,dSA,Vf,SAf
	  REAL :: JXm_T,JNm_T,JHGm_T, JHT0_T,JST0_T,JCPm_T,JSGm_T,kNPQ_T,kROS_T

!################## Generic ####################
	  REAL, PARAMETER:: yc=0.8  		!Yield of biomass formation from carbon | C-mol mol C −1
	  REAL, PARAMETER:: Ti=293.15    	! Reference Temperature
!################## HOST #######################
	  REAL, PARAMETER:: JHGm=0.0244	!Maximum specific growth rate of host | C-mol H C-mol H −1 d −1
	  REAL, PARAMETER:: Jxm=2.0916		!Maximum prey assimilation rate from host feeding
 	  REAL, PARAMETER:: Kx=4.766e-7 	!Half-saturation constant for prey assimilation
	  REAL, PARAMETER:: JHT0=0.0209	!Maintenance rate of host biomass
	  REAL, PARAMETER:: Kn=1.50e-6	!Half-saturation constant for host DIN uptake
	  REAL, PARAMETER:: kv=16.9		!L c-mol H-1
	  REAL, PARAMETER:: kSA=3.4551e3!		
!################## Symbiont ###################
 	  REAL, PARAMETER:: JSGm=0.0210 	!Maximum specific growth rate of symbiont
	  REAL, PARAMETER:: JST0=0.0094	    !Maintenance rate of symbiont biomass
	  REAL, PARAMETER:: JCPm=13.3937	!Maximum specific photosynthesis rate of symbiont | mol C C-mol S−1 d−1
	  REAL, PARAMETER:: ycl=0.1  		!Quantum yield of photosynthesis | mol C mol photons −1
	  REAL, PARAMETER:: a=1.34		   	!Effective light-absorbing cross-section of symbiont | m 2 C-mol S −1
	  REAL, PARAMETER:: kCO2=3.0383 	!Efficacy of CO 2 delivery to photosynthesis by host CCMs
	  REAL, PARAMETER:: sigma_CH=0.1	!Proportion host metabolic CO2 recycled to photosynthesis
	  REAL, PARAMETER:: sigma_CS=0.9 	!Proportion symbiont metabolic CO 2 recycled to photosynthesis
	  REAL, PARAMETER:: JNm=0.0353	    !Maximum host DIN uptake rate
	  REAL, PARAMETER:: sigma_NH=0.9  	!Proportion N turnover recycled in host
	  REAL, PARAMETER:: nNH=0.18		!N:C molar ratio in host biomass
	  REAL, PARAMETER:: nNS=0.13		!N:C molar ratio in symbiont biomass
	  REAL, PARAMETER:: nNX=0.2		    !N:C molar ratio in prey biomass
	  REAL, PARAMETER:: kNPQ=113.7425   !NPQ capacity of symbiont
	  REAL, PARAMETER:: kROS=83.6141	!Excess photon energy that doubles ROS production, relative to baseline levels
 	  REAL, PARAMETER:: sigma_NS=0.9  	!Proportion N turnover recycled in symbiont
	  REAL, PARAMETER:: b=6.0285		!Scaling parameter for bleaching response
	  REAL, PARAMETER:: Ccell=85e-12 	!Carbon content of 1 cell - g// 85pg/cell - Muscatine 1984


  	  REAL, PARAMETER:: TSa=5800.  		! Arrhenius Temperature for Symbiont
      REAL, PARAMETER:: TSal=2.0096e4	! Rate of decrease of lower boundary
      REAL, PARAMETER:: TSah=5.0710e4	! Rate of decrease of upper boundary
      REAL, PARAMETER:: TSl=280.6		! Lower boundary of tolerance range for Symbiont
      REAL, PARAMETER:: TSh=302.1		! Upper boundary of tolerance range for Symbiont
      
	  REAL, PARAMETER:: THa=5800.  		! Arrhenius Temperature for Host
      REAL, PARAMETER:: THal=2.0108e4	! Rate of decrease of lower boundary
      REAL, PARAMETER:: THah=5.0629e4	! Rate of decrease of upper boundary
      REAL, PARAMETER:: THl=280.6		! Lower boundary of tolerance range for Host
      REAL, PARAMETER:: THh=302.08		! Upper boundary of tolerance range for Host
	

	  end module param



