#ifndef SEDIMENT_OPTIONS_H
#define SEDIMENT_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

CBOP
C !ROUTINE: SEDIMENT_OPTIONS.h
C !INTERFACE:
C #include "SEDIMENT_OPTIONS.h"

C !DESCRIPTION:
C *==================================================================*
C | CPP options file for pkg "sediment" (suspended sediment with
C | settling, deposition, erosion and a bed store; for use with pkg
C | wetdry on tidal flats, but works without it)
C *==================================================================*
CEOP

#ifdef ALLOW_SEDIMENT
C (no optional features yet)
#endif /* ALLOW_SEDIMENT */
#endif /* SEDIMENT_OPTIONS_H */
