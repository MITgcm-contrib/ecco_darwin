#ifdef ALLOW_BBL

C     bbl_wvel    :: vertical detrainment velocity of a dense BBL (m/s)
C     bbl_hvel    :: default horizontal velocity of BBL (m/s); upper
C                    bound of the speed when bbl_useNofSpeed=T
C     bbl_initEta :: default initial thickness of BBL (m)
C     bbl_useNofSpeed :: compute the lateral BBL speed from the density
C                    contrast and slope, u = bbl_nofCoeff*gPrime*slope/f
C     bbl_nofCoeff:: fraction of the Nof speed gPrime*slope/f that drains
C                    downslope (default 1; Campin & Goosse 1999 / MOM
C                    use 1/3, which gave too little export in tests)
C     bbl_fMin    :: lower bound of |f| in the speed (1/s)
C     bbl_selectEntrain :: entrainment of ambient water into the BBL
C                    0 = none; 1 = Turner (1986) E(Ri), w_E = E*U, with
C                    Ri = gPrime*bbl_eta/U^2 and U the BBL speed
C     bbl_selectDetrain :: detrainment of a dense BBL into its cell
C                    0 = at bbl_wvel everywhere (original scheme);
C                    1 = plume-like: at bbl_wvel*max(0,1-uOut/bbl_uDet),
C                    uOut the fastest outflow speed of the cell, so a
C                    descending plume keeps its water and loses it only
C                    downslope and when it reaches neutral density
C     bbl_uDet    :: outflow speed (m/s) above which a dense BBL no
C                    longer detrains (bbl_selectDetrain = 1)
C     bbl_gpTaper :: reduced-gravity scale (m/s^2) over which the BBL
C                    fades out as it approaches neutral density: for
C                    0 < gPrime < bbl_gpTaper entrainment is scaled by
C                    gPrime/bbl_gpTaper and the BBL dissolves into its
C                    cell at the complementary rate. 0 (default) = the
C                    original switch (instant mix when gPrime <= 0),
C                    which is discontinuous at gPrime = 0
C     bbl_tauMix  :: timescale (s) on which a neutral or light BBL
C                    dissolves into its cell; 0 (default) = within one
C                    time step
C     bbl_coldStart :: on a restart, do not read pickup_bbl but start the
C                    BBL as at iteration 0 (BBL T/S = bottom cell, zero
C                    thickness, or bbl_*File if set)
C     bbl_doPtracers :: the BBL also carries the passive tracers of
C                    pkg/ptracers (default: .TRUE. when usePTRACERS)
C     bbl_skipIceCavity :: no BBL in ice-shelf cavity columns (kSurfC > 1);
C                    default .FALSE. In thin LLC90 cavity columns the
C                    BBL couples to the ice-shelf melt and the adjoint
C                    grows without bound (Getz Ice Shelf, Oct 2026)
      LOGICAL             bbl_useNofSpeed, bbl_coldStart, bbl_doPtracers
      LOGICAL             bbl_skipIceCavity
      COMMON /BBL_PARM_L/ bbl_useNofSpeed, bbl_coldStart,
     &                    bbl_doPtracers, bbl_skipIceCavity
      INTEGER             bbl_selectEntrain, bbl_selectDetrain
      COMMON /BBL_PARM_I/ bbl_selectEntrain, bbl_selectDetrain
      _RL                 bbl_wvel, bbl_hvel, bbl_initEta
      _RL                 bbl_nofCoeff, bbl_fMin
      _RL                 bbl_gpTaper, bbl_tauMix, bbl_uDet
      COMMON /BBL_PARM_R/ bbl_wvel, bbl_hvel, bbl_initEta,
     &                    bbl_nofCoeff, bbl_fMin,
     &                    bbl_gpTaper, bbl_tauMix, bbl_uDet

C     file names for initial conditions:
      CHARACTER*(MAX_LEN_FNAM) bbl_thetaFile
      CHARACTER*(MAX_LEN_FNAM) bbl_saltFile
      CHARACTER*(MAX_LEN_FNAM) bbl_etaFile
      COMMON /BBL_PARM_C/ bbl_thetaFile, bbl_saltFile, bbl_etaFile

C     BBL 2-dim. fields
      _RL bbl_theta (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL bbl_salt  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL bbl_eta   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL bbl_rho_nr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /BBL_STATE_2D/ bbl_theta, bbl_salt, bbl_eta, bbl_rho_nr

C     bbl_kBot :: bottom level of the BBL of each column: kLowC, or 0
C                 (no BBL) in ice-shelf cavities if bbl_skipIceCavity
      INTEGER bbl_kBot  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /BBL_KBOT/ bbl_kBot

C     BBL lateral-exchange limiter factors (scratch, valid in halos
C     after EXCH): fraction of the face volumes a cell may give / take
      _RL bbl_sclOut(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL bbl_sclIn (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /BBL_LIMITER/ bbl_sclOut, bbl_sclIn

C     BBL tendency terms
      _RL bbl_TendTheta(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL bbl_TendSalt(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /BBL_TENDENCY/ bbl_TendTheta, bbl_TendSalt

#endif /* ALLOW_BBL */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
