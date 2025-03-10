#ifdef ALLOW_DARWIN

CBOP
C     !ROUTINE: DARWIN_EXF_FIELDS_FIELDS.h
C     !INTERFACE:
C #include DARWIN_EXF_FIELDS_FIELDS.h

C     !DESCRIPTION:
C Contains fields for darwin package read through exf
C
C Requires: SIZE.h

C--   COMMON /DARWIN_FIELDS_C/
C     ventHe3file   :: file with He3 flux from hydrothermal vents (mmol He/m2/s)
      COMMON /DARWIN_FIELDS_C/
     &    PARfile,
     &    ironfile,
     &    icefile,
     &    windfile,
     &    pCO2file,
     &    ventHe3file,
     &    DOCrunofffile,
     &    rDOCrunofffile,
     &    CDOMrunofffile,
     &    DONrunofffile,
     &    DOPrunofffile,
     &    NO3runofffile,
     &    NO2runofffile,
     &    NH4runofffile,
     &    PO4runofffile,
     &    DSirunofffile,
     &    POCrunofffile,
     &    POPrunofffile,
     &    PONrunofffile,
     &    DICrunofffile,
     &    ALKrunofffile

      CHARACTER*128 PARfile
      CHARACTER*128 ironfile
      CHARACTER*128 icefile
      CHARACTER*128 windfile
      CHARACTER*128 pCO2file
      CHARACTER*128 ventHe3file
      CHARACTER*128 DOCrunofffile
      CHARACTER*128 rDOCrunofffile
      CHARACTER*128 CDOMrunofffile
      CHARACTER*128 DONrunofffile
      CHARACTER*128 DOPrunofffile
      CHARACTER*128 NO3runofffile
      CHARACTER*128 NO2runofffile
      CHARACTER*128 NH4runofffile
      CHARACTER*128 PO4runofffile
      CHARACTER*128 DSirunofffile
      CHARACTER*128 POCrunofffile
      CHARACTER*128 POPrunofffile
      CHARACTER*128 PONrunofffile
      CHARACTER*128 DICrunofffile
      CHARACTER*128 ALKrunofffile

C--   COMMON /DARWIN_FIELDS_R/
C     ventHe3   :: He3 flux from hydrothermal vents (mmol He/m2/s)
      COMMON /DARWIN_FIELDS_R/
     &    PAR0, PAR1, surfPAR,
     &    iron0, iron1, inputFe,
     &    ice0, ice1, iceFrac,
     &    wind0, wind1, windSpeed,
     &    pCO20, pCO21, atmospCO2,
     &    ventHe30, ventHe31, ventHe3,
     &    DOCrunoff0, DOCrunoff1, DOCrunoff,
     &    rDOCrunoff0, rDOCrunoff1, rDOCrunoff,
     &    CDOMrunoff0, CDOMrunoff1, CDOMrunoff,
     &    DONrunoff0, DONrunoff1, DONrunoff,
     &    DOPrunoff0, DOPrunoff1, DOPrunoff,
     &    NO3runoff0, NO3runoff1, NO3runoff,
     &    NO2runoff0, NO2runoff1, NO2runoff,
     &    NH4runoff0, NH4runoff1, NH4runoff,
     &    PO4runoff0, PO4runoff1, PO4runoff,
     &    DSirunoff0, DSirunoff1, DSirunoff,
     &    POCrunoff0, POCrunoff1, POCrunoff,
     &    POPrunoff0, POPrunoff1, POPrunoff,
     &    PONrunoff0, PONrunoff1, PONrunoff,
     &    DICrunoff0, DICrunoff1, DICrunoff,
     &    ALKrunoff0, ALKrunoff1, ALKrunoff

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

      _RL ventHe30(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ventHe31(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ventHe3(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL DOCrunoff0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DOCrunoff1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DOCrunoff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL rDOCrunoff0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL rDOCrunoff1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL rDOCrunoff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL CDOMrunoff0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL CDOMrunoff1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL CDOMrunoff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL DONrunoff0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DONrunoff1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DONrunoff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL DOPrunoff0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DOPrunoff1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DOPrunoff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL NO3runoff0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL NO3runoff1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL NO3runoff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL NO2runoff0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL NO2runoff1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL NO2runoff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL NH4runoff0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL NH4runoff1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL NH4runoff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL PO4runoff0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL PO4runoff1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL PO4runoff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL DSirunoff0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DSirunoff1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DSirunoff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL POCrunoff0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL POCrunoff1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL POCrunoff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL POPrunoff0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL POPrunoff1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL POPrunoff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL PONrunoff0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL PONrunoff1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL PONrunoff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL DICrunoff0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DICrunoff1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DICrunoff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL ALKrunoff0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ALKrunoff1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ALKrunoff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

CEOP

#endif /* ALLOW_DARWIN */
