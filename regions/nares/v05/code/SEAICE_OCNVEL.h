C     SEAICE_OCNVEL.h (Nares 2026-11-15, Dustin: fix (a)). Ocean velocity used in the
C     ice-ocean drag (drag coefficients, LSR solver, ocean stress) = thickness-weighted
C     mean of uVel/vVel over the top SIdragDepth metres, instead of level kSrf = 1
C     (drF 1 m). Filled each time step by S/R SEAICE_OCNVEL_AVG (z-coordinates only).
      _RL SIdragDepth
      PARAMETER ( SIdragDepth = 10. _d 0 )
      _RL uOcnSI(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vOcnSI(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /SEAICE_OCNVEL_R/ uOcnSI, vOcnSI
