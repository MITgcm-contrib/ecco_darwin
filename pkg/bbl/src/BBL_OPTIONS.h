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

#endif /* ALLOW_BBL */
#endif /* BBL_OPTIONS_H */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
