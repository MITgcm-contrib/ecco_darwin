C     !ROUTINE: OFFLINE_SWITCH.h
C -------------------------------
C   OFFLINE_SWITCH.h
C  variable for switching on/off some calculations
C -------------------------------
C
C  Patched (local override) to add offlineLoadGGL90, mirroring
C  offlineLoadKPP.

C     offlineLoadGMRedi :: load from file GMRedi tensor (do not compute it)
C     offlineLoadKPP    :: load from file KPP mixing coeff (do not compute it)
C     offlineLoadGGL90  :: load from file GGL90 diffusivity (do not compute it)
C     offlineLoadConvec :: load from file Convective mixing (do not compute it)
      COMMON /OFFLINE_SWITCH_L/
     &       offlineLoadGMRedi, offlineLoadKPP, offlineLoadGGL90,
     &       offlineLoadConvec
      LOGICAL offlineLoadGMRedi
      LOGICAL offlineLoadKPP
      LOGICAL offlineLoadGGL90
      LOGICAL offlineLoadConvec

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
