#ifdef ALLOW_DARWIN

CBOP
C    !ROUTINE: DARWIN_DIAGS.h
C    !INTERFACE:
C #include DARWIN_DIAGS.h

C    !DESCRIPTION:
C Contains indices into diagnostics array

      integer iPP
      integer iNfix
      integer iDenit
      integer iDenitN
      integer iConsPO4
      integer iConsSi
      integer iConsFe 
      integer iConsDIN
      integer iPPplank
      integer iConsNO3
      integer iConsNO2
      integer iConsNH4
      integer iConsumDIC
      integer iConsumDIC_PIC
      integer iRespirDIC
	  integer iReminDIC_DOC
	  integer iReminDIC_POC
	  integer iDisscDIC_PIC
      integer iGRplank
      integer iGrGn
      integer darwin_nDiag
      
      PARAMETER(iPP=     1)
      PARAMETER(iNfix=   2)
      PARAMETER(iDenit=  3)
      PARAMETER(iDenitN= 4)
      PARAMETER(iConsPO4=5)
      PARAMETER(iConsSi= 6)
      PARAMETER(iConsFe= 7)
      PARAMETER(iConsDIN=8)
      PARAMETER(iPPplank=9)
      PARAMETER(iConsNO3=10)
      PARAMETER(iConsNO2=11)
      PARAMETER(iConsNH4=12)
      PARAMETER(iConsumDIC=13)
      PARAMETER(iConsumDIC_PIC=14)
      PARAMETER(iRespirDIC=15)
	  PARAMETER(iReminDIC_DOC=16)
	  PARAMETER(iReminDIC_POC=17)
	  PARAMETER(iDisscDIC_PIC=18)
      
#ifdef DARWIN_DIAG_PERTYPE
      PARAMETER(iGRplank=iPPplank+nplank)
      PARAMETER(iGrGn=iGRplank+nplank)
      PARAMETER(darwin_nDiag=iGrGn+nplank-1)
#else
      PARAMETER(iGRplank=iDisscDIC_PIC)
      PARAMETER(iGrGn=iGRplank)
      PARAMETER(darwin_nDiag=iGrGn-1)
#endif

CEOP
#endif /* ALLOW_DARWIN */
