C CPP options file for BBL
C Use this file for selecting options within package "BBL"

#ifndef BBL_OPTIONS_H
#define BBL_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

#ifdef ALLOW_BBL
C Place CPP define/undef flag here

C Print global sums of BBL T/S tendency x volume x tracer time step
C (should be zero for a conservative lateral exchange)
#undef BBL_CHECK_BUDGET

C Adjoint only: no sensitivity through the BBL speeds and the
C entrainment speed (passive in the adjoint, as with
C GMREDI_WITH_STABLE_ADJOINT); the full adjoint can grow without bound
C on steep slopes (ECCO LLC90, East Greenland shelf, Oct 2026)
#undef BBL_WITH_STABLE_ADJOINT

#endif /* ALLOW_BBL */
#endif /* BBL_OPTIONS_H */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
