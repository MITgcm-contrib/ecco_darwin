#ifdef ALLOW_BBL
#ifdef ALLOW_PTRACERS
C     Passive tracers (pkg/ptracers) carried by the bottom boundary layer.
C     Requires PTRACERS_SIZE.h to be included before this file.

C     bbl_ptr         :: BBL concentration of each passive tracer
C     bbl_TendPtr     :: BBL tendency of each passive tracer in the
C                        bottommost wet cell (tracer units / s)
      _RL bbl_ptr    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,
     &                PTRACERS_num)
      _RL bbl_TendPtr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,
     &                PTRACERS_num)
      COMMON /BBL_PTR_FIELDS/ bbl_ptr, bbl_TendPtr

#endif /* ALLOW_PTRACERS */
#endif /* ALLOW_BBL */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
