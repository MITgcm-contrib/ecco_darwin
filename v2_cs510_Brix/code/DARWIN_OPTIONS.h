#ifndef DARWIN_OPTIONS_H
#define DARWIN_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#ifdef ALLOW_DARWIN

#include "CPP_OPTIONS.h"

CBOP
C    !ROUTINE: DARWIN_OPTIONS.h
C    !INTERFACE:

C    !DESCRIPTION:
C options for darwin package
CEOP

cSD#define READ_PAR    
cSD#undef  USE_QSW
cHB ----
#undef  READ_PAR
#define USE_QSW
#define USE_EXFWIND
#define USE_EXFCO2
cHB -----
#define MINFE
#undef  NUT_SUPPLY
#undef  CONS_SUPP
#undef  OLD_GRAZE
#undef  OLD_NSCHEME
#undef  ALLOW_MUTANTS
#define PORT_RAND
#undef  OLDSEED

#define TEMP_VERSION 2
#undef TEMP_RANGE 

#undef TWO_SPECIES_SETUP
#define NINE_SPECIES_SETUP

#undef  RELAX_NUTS
#undef  FLUX_NUTS

#define CALC_RATE_TOTALS

#define IRON_SED_SOURCE
#define IRON_SED_SOURCE_VARIABLE
#define PART_SCAV

#undef  ALLOW_DIAZ
#undef  ALLOW_DENIT
#undef  DENIT_RELAX

#define ALLOW_CARBON
cSD#define ALLOW_OLD_VIRTUALFLUX
#undef ALLOW_OLD_VIRTUALFLUX


#undef  GEIDER
#undef  OASIM
#undef  WAVEBANDS
#undef  WAVEBAND_4TYPE
#undef  WAVEBAND_10TYPE

#undef  DYNAMIC_CHL
#undef  DAR_CALC_ACDOM
#undef  DAR_RADTRANS
#undef  DAR_RADTRANS_USE_MODEL_CALENDAR
C truncation to 2 downward decreasing modes a la Aas
#undef  DAR_RADTRANS_DECREASING
C iterative solution
#undef  DAR_RADTRANS_ITERATIVE
C use rmus for all components to convert to scalar irradiance
C (not recommended)
#undef  DAR_RADTRANS_RMUS_PAR
#undef  DAR_NONSPECTRAL_BACKSCATTERING_RATIO


#undef  CHECK_CONS
#undef  DAR_DIAG_RSTAR
#undef  DAR_DIAG_DIVER
#undef  DAR_DIAG_GROW
#undef  DAR_DIAG_ACDOM
#undef  DAR_DIAG_ABSORP
#undef  DAR_DIAG_SCATTER
#undef  DAR_DIAG_IRR

C diagnostic chlorophyll
#undef  DAR_DIAG_CHL

C average PAR daily and store previous day
#undef  ALLOW_PAR_DAY

C dependencies
#ifdef DAR_DIAG_CHL
#define ALLOW_PAR_DAY
#endif

#endif /* ALLOW_DARWIN */
#endif /* DARWIN_OPTIONS_H */
