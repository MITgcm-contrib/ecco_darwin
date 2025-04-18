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
      integer iConsDIC
      integer iConsDIC_PIC
      integer iRespirDIC
      integer iReminDIC_DOC
      integer iReminDIC_POC
      integer iDisscDIC_PIC 
      integer iConsALK
      integer iConsO2
C-MM      
      integer igrateMAG
      integer igTMAG
      integer igQMAG
      integer igEMAG
      integer igHMAG  
C-MM
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
      PARAMETER(iConsNO3=9)
      PARAMETER(iConsNO2=10)
      PARAMETER(iConsNH4=11)
      PARAMETER(iConsDIC=12)
      PARAMETER(iConsDIC_PIC=13)
      PARAMETER(iRespirDIC=14)
      PARAMETER(iReminDIC_DOC=15)
      PARAMETER(iReminDIC_POC=16)
      PARAMETER(iDisscDIC_PIC=17)
      PARAMETER(iConsALK=18)
      PARAMETER(iConsO2=19)
      PARAMETER(iPPplank=20)
C#ifdef DARWIN_ALLOW_MACROALGAE
C
C-MM 
C Adding new indices for Kelp diagnostics
C      PARAMETER(igrateMAG=21)
C      PARAMETER(igTMAG=22)	
C      PARAMETER(igQMAG=23)
C      PARAMETER(igEMAG=24)
C      PARAMETER(igHMAG=25)
C#endif

#ifdef DARWIN_DIAG_PERTYPE
      PARAMETER(iGRplank=iPPplank+nplank)
C      PARAMETER(iGRplank=igHMAG+nplank)
      PARAMETER(iGrGn=iGRplank+nplank)
      PARAMETER(darwin_nDiag=iGrGn+nplank-1)
#else
      PARAMETER(iGRplank=iPPplank)
C      PARAMETER(iGRplank=igHMAG)        
      PARAMETER(iGrGn=iGRplank)
      PARAMETER(darwin_nDiag=iGrGn-1)
#endif

CEOP
#endif /* ALLOW_DARWIN */
