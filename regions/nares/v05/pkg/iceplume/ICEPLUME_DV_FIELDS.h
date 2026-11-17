#ifdef ALLOW_ICEPLUME

C This header file will keep track of the iceplume fields
C so that they can be accessed by the output routine of
C the diagnostics_vec package.
C This is a bit of a work-around since iceplume is not
C currently set up to keep track of these variables outside
C of the iceplume_calc routine

C Header file pkg/ICEPLUME_DV_FIELDS

      COMMON /ICEPLUME_DV_FIELDS/
     &     icefrntA_dv, icefrntM_dv, icefrntR_dv,
     &     icefrntS_dv, icefrntT_dv, icefrntW_dv
      _RL icefrntA_dv  (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RL icefrntM_dv  (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RL icefrntR_dv  (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RL icefrntS_dv  (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RL icefrntT_dv  (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RL icefrntW_dv  (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)

#endif /* ALLOW_ICEPLUME */
