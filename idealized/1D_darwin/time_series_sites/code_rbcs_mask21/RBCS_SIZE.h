#ifdef ALLOW_RBCS

CBOP
C    !ROUTINE: RBCS_SIZE.h
C    !INTERFACE:

C    !DESCRIPTION:
C Contains RBCS array size (number of tracer mask)
C 1-D Darwin (2026-10-09): maskLEN = 2 + 19 ptracers, so each BGC tracer can have its own
C relaxation mask (relaxMaskFile(2+iTr)); used for the Papa SiO2-only supply test.
CEOP

C---  RBCS Parameters:
C     maskLEN :: number of mask to read
      INTEGER maskLEN
      PARAMETER( maskLEN = 21 )

#endif /* ALLOW_RBCS */
