#ifdef ALLOW_DARWIN

CBOP
C     !ROUTINE: DARWIN_EXF_FIELDS_FIELDS.h
C     !INTERFACE:
C #include DARWIN_EXF_FIELDS_FIELDS.h

C     !DESCRIPTION:
C Contains fields for darwin package read through exf
C
C Requires: SIZE.h

      COMMON /darwin_fields_c/
     &    PARfile,
     &    ironfile,
     &    icefile,
     &    windfile,
     &    pCO2file,
     &    DOCrunofffile,
     &    tDOCrunofffile,
     &    DONrunofffile,
     &    DOPrunofffile,
     &    DSirunofffile,
     &    DICrunofffile,
     &    ALKrunofffile

      CHARACTER*128 PARfile
      CHARACTER*128 ironfile
      CHARACTER*128 icefile
      CHARACTER*128 windfile
      CHARACTER*128 pCO2file
	    CHARACTER*128 DOCrunofffile
	    CHARACTER*128 tDOCrunofffile
      CHARACTER*128 DONrunofffile
      CHARACTER*128 DOPrunofffile
      CHARACTER*128 DSirunofffile
      CHARACTER*128 DICrunofffile
      CHARACTER*128 ALKrunofffile

      COMMON /darwin_fields_r/
     &    PAR0, PAR1, surfPAR,
     &    iron0, iron1, inputFe,
     &    ice0, ice1, iceFrac,
     &    wind0, wind1, windSpeed,
     &    pCO20, pCO21, atmospCO2,
     &    DOCrunoff0,DOCrunoff1,DOCrunoff,
     &    tDOCrunoff0,tDOCrunoff1,tDOCrunoff,
     &    DONrunoff0,DONrunoff1,DONrunoff,
     &    DOPrunoff0,DOPrunoff1,DOPrunoff,
     &    DSirunoff0,DSirunoff1,DSirunoff,
     &    DICrunoff0,DICrunoff1,DICrunoff,
     &    ALKrunoff0,ALKrunoff1,ALKrunoff

      _RL PAR0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL PAR1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL surfPAR(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL iron0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL iron1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL inputFe(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL ice0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ice1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL iceFrac(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL wind0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL wind1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL windSpeed(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL pCO20(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL pCO21(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL atmospCO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL DOCrunoff0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DOCrunoff1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DOCrunoff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL tDOCrunoff0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL tDOCrunoff1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL tDOCrunoff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL DONrunoff0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DONrunoff1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DONrunoff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL DOPrunoff0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DOPrunoff1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DOPrunoff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL DSirunoff0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DSirunoff1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DSirunoff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL DICrunoff0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DICrunoff1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DICrunoff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL ALKrunoff0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ALKrunoff1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ALKrunoff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

CEOP

#endif /* ALLOW_DARWIN */
