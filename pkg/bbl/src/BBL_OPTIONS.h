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

C Adjoint only: the BBL flow is passive in the adjoint (speeds,
C entrainment speed, face volumes incl. limiter, fade-out ramp), as with
C GMREDI_WITH_STABLE_ADJOINT, so the adjoint only carries tracers along
C the BBL flow of the forward run. The full adjoint can grow without
C bound on steep slopes (ECCO LLC90, Oct 2026)
#undef BBL_WITH_STABLE_ADJOINT

#endif /* ALLOW_BBL */
#endif /* BBL_OPTIONS_H */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
