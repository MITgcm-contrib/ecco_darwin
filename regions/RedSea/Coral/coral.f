	  Program Coral

	  USE param
	  USE sub

	  implicit none
!####################################################

		init_day_r=1*dt

  	  tt= nyear*years  ! total time in days
         ts=int(tt*(1./dt))  ! total steps
         init_day=int(init_day_r/dt)
     
	  OPEN(100,file='Host.dat')
  	  OPEN(200,file='Symbiont.dat')
  	  OPEN(300,file='Symbiosis.dat')
  	  OPEN(400, file='Forcing.dat')


	  H=0.0031
	  S=0.00013
          SA=10.75
         
	
!#### External Forcning ############################
	


	    do i=1,1825
        open(150, file = 'Forcing/PAR_5m_2008_2012_Th.txt',
     & form='formatted')
        read(150,'(e12.5)'),PARdata
        close(150)
      end do

	    do i=1,1825
        open(150, file = 'Forcing/SST_2008_2012_Th.txt',
     & form='formatted')
        read(150,'(e12.5)'),SSTdata
        close(150)
      end do


!####################################################
       do n= 15,ts !

    	  tr=n*dt 
	   	  

!       trr=mod(int(tr),trr)

		trr=int(tr)+1
	    		write(*,*), trr

!	    if (trr .eq. 0) then
!
!	   	    trr=365

!		end if

         L=PARdata(trr)
	   DIN=1e-7
	   X=2e-7
	
!#######  Temperature Effect #######    

      T=SSTdata(trr)
      T=SSTdata(trr)
      T=T+273.15
      aT=exp((TSa/Ti)-(TSa/T))
      bT=exp((TSal/T)-(TSal/TSl))
      cT=exp((TSah/TSh)-(TSah/T))
      kST= aT /(1. + bT  +cT)


	  aT2=exp((THa/Ti)-(THa/T))
      bT2=exp((THal/T)-(THal/THl))
      cT2=exp((THah/THh)-(THah/T))
      kHT= aT2 /(1. + bT2  +cT2);


!      kHT=exp(THa/Ti - THa/T)


	  JXm_T=JXm*kHT
	  JNm_T=JNm*kHT
	  JHGm_T=JHGm*kHT
	  JHT0_T=JHT0*kHT	
	  JST0_T=JST0*kST
	  JCPm_T=JCPm*kST
  	  JSGm_T=JSGm*kST
  	  kNPQ_T=kNPQ!*kST
	  kROS_T=kROS*kST

      JX=(Jxm_T*X)/(X+Kx)
	  JN=(JNm_T*DIN)/(DIN+Kn)
	  rNH=sigma_NH*nNH*JHT0_T
	  rNS=sigma_NS*nNS*JST0_T
	     	  
	  if (n .eq.1) then 
	
	 	  
	  JL= L*a
	  JHT=JHT0_T
  	  JEC=10.
	  rCH=sigma_CH*JHT0_T
	  JCO2=kCO2*JEC
      
      call synth(JCP,JL*ycl,(JCO2+rCH)*H/S+rCS,JCPm_T)
      JCP=JCP/CROS 

	  JEL=(JL-JCP*1./ycl)
	  JNPQ=kNPQ_T
	  CROS=1.
	  JHG=0.25
	  pN=JN
!	  dH=JHGm_T
!	  dS=JSGm_T  
	  JSG=JSGm_T
	  pC=JCP
	  rCS=JST0_T*sigma_CS
	  JST=JST0_T

	  elseif (n .gt. 1) then
	

!###### Symbiont ##################################

	  Al=1.26 +1.39*exp(-6.48*(S/H))
	  JL= Al*L*a
  	  rCS=sigma_CS*JST0_T+(1-yc)*JSG/yc
      call synth(JCP,JL*ycl,(JCO2+rCH)*H/S+rCS,JCPm_T)
      JCP=JCP/CROS 
	  JEL=JL-(JCP/ycl)

	  if (JEL .lt.0) then
	   JEL=0.
	  end if

      JNPQ=1./((1./kNPQ)+(1./JEL))
!      JNPQ=(((1./kNPQ_T)+(1./JEL)))
	  
	  if ((JEL-JNPQ) .lt.0) then
	  CROS=1.

	  else 	

	  CROS=1+(JEL-JNPQ)/kROS_T

	  end if	 

      call synth(JSG, yC*JCP,(pN*H/S+rNS)/nNS,jSGm_T)

      
      pc= JCP-JSG/yc

	  if (pC .lt.0) then
	  pC=0.
	  end if

 	  JST=JST0_T*(1+b*(CROS-1))

!###### Host ###############################

      call synth(JHG,yC*(pC*S/H + jX),(jN+nNX*jX + rNH)/nNH,jHGm_T)

	  pN=JN+(nNX*Jx)+rNH-(nNH*JHG)!*(1./yc))

	  if (pN .lt.0) then
	  pN=0. 
	  !write(*,*) 'pN=0'
	  end if

  	  JEC=Jx+pc*(S/H)-JHG/yc


  	  if (JEC .lt.0) then
	  JEC=0.
	  end if

  	  rCH=sigma_CH*(JHT0_T+(1-yc)*JHG/yc)
	  JCO2=kCO2*JEC
	  JHT=JHT0_T

	  endif


       dS=(JSG-JST)*S
 	   dH=(JHG-JHT)*H

	   Sf=S+dS*dt
	   Hf=H+dH*dt


  	   if (dH .gt.0) then     
           dV=kV*dH;
           else
           dV=0;
           endif


  	   if (dH .gt.0) then     
           dSA=kSA*dH;
           else
           dSA=0;
           endif

	   Vf=V+dV*dt
	   SAf=SA+dSA*dt


	   S_g=S*12.01
	   S_cells=S_g/Ccell
	   S_cells_cm2=S_cells/SA

	   SGR=(log(Hf)-log(H))*100/dt

        write(100,1000)n,H,Jx,Jn,JHG,JHT,dH,dt,Vol,SA,SGR
1000    format(i7,10(2x,e12.5)) 

        write(200,2000)n,S,JL,JCO2,JeL,JNPQ,JCP,CROS,JST,JSG,dS,
     & S_cells_cm2
2000    format(i7,11(2x,e12.5)) 

 		write(300,3000)n,S/H,pC,pN,JeC,rNH,rCH,rNS,rCS
3000    format(i7,8(2x,e12.5))  

 		write(400,4000)n,X,DIN,L,T,kHT,kST
4000    format(i7,6(2x,e12.5))  



	  S=Sf
	  H=Hf
	  V=Vf
	  SA=SAf

      enddo

      stop

      end program coral


