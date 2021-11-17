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
      integer iPPplank
      integer iPCplank
      integer iConsNO3
      integer iConsNO2
      integer iConsNH4
      integer iConsPO4
      integer iConsSi
      integer iConsFe
      integer iConsDIC
      integer iConsALK
      integer iConsDIC_PIC
      integer iRespirDIC
      integer iReminDIC_DOC
      integer iReminDIC_POC
      integer iDisscDIC_PIC
      integer iGRplank
      integer iGrGn
      integer darwin_nDiag

      PARAMETER(iPP=      1)
      PARAMETER(iNfix=    2)
      PARAMETER(iDenit=   3)
      PARAMETER(iDenitN=  4)
      PARAMETER(iConsNO3= 5)
      PARAMETER(iConsNO2= 6)
      PARAMETER(iConsNH4= 7)
      PARAMETER(iConsPO4= 8)
      PARAMETER(iConsSi=  9)
      PARAMETER(iConsFe= 10)
      PARAMETER(iConsDIC=11)
      PARAMETER(iConsALK=12)
      PARAMETER(iConsDIC_PIC=13)
      PARAMETER(iRespirDIC=14)
      PARAMETER(iReminDIC_DOC=15)
      PARAMETER(iReminDIC_POC=16)
      PARAMETER(iDisscDIC_PIC=17)
      PARAMETER(iPPplank=18)
#ifdef DARWIN_DIAG_PERTYPE
      PARAMETER(iPCplank=iPPplank+nplank)
      PARAMETER(iGRplank=iPCplank+nplank)
      PARAMETER(iGrGn=iGRplank+nplank)
      PARAMETER(darwin_nDiag=iGrGn+nplank-1)
#else
      PARAMETER(iPCplank=iPPplank)
      PARAMETER(iGRplank=iPCplank)
      PARAMETER(iGrGn=iGRplank)
      PARAMETER(darwin_nDiag=iGrGn-1)
#endif

CEOP
#endif /* ALLOW_DARWIN */
