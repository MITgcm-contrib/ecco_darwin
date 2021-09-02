#ifdef ALLOW_RADTRANS

CBOP
C    !ROUTINE: RADTRANS_SIZE.h
C    !INTERFACE:
C #include RADTRANS_SIZE.h

C    !DESCRIPTION:
C Contains dimensions and index ranges for radtrans model.

      integer nlam
      parameter(nlam=13)

CEOP
#endif /* ALLOW_RADTRANS */
