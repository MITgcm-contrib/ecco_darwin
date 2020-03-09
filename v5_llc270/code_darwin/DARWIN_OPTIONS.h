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

cHB ----
#undef  READ_PAR
#define USE_QSW
#define USE_EXFWIND
#define ALLOW_CARBON
#define USE_EXFCO2
cHB -----
#define MINFE
#undef  NUT_SUPPLY
#undef  CONS_SUPP
#undef  OLD_GRAZE
#undef  ALLOW_DIAZ
#undef  ALLOW_DENIT
#undef  DENIT_RELAX
#undef  OLD_NSCHEME
#undef  ALLOW_MUTANTS
#define PORT_RAND
#undef  OLDSEED
#undef  CALC_RATE_TOTALS

#undef NOTEMP
#define TEMP_VERSION 2
#undef TEMP_RANGE 

#undef TWO_SPECIES_SETUP
#define NINE_SPECIES_SETUP

#define CALC_RATE_TOTALS
#define IRON_SED_SOURCE
#define IRON_SED_SOURCE_VARIABLE
#define PART_SCAV

#define  GEIDER
#undef  OASIM
#undef  WAVEBANDS
#define  DYNAMIC_CHL
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
C define this to turn of reading of phyto backscattering spectra
C and revert to fixed backscat ratios darwin_bbphy(nabp) set in data.darwin
#undef  DAR_NONSPECTRAL_BACKSCATTERING_RATIO

C Define omegaC-dependent dissolution rate
C following Naviaux et al. 2019 (Marine Chemistry)
#undef NAVIAUX_DISSOLUTION

C allow DIC and ALK flux from bottom sediments
C following Sulpis et al. 2018 (PNAS)
#undef ALLOW_SED_DISS_FLUX

c allow DIC, ALK, and Ca surface flux
#undef ALLOW_ADKINS_SURF_FLUX

C compute CO2 flux budget terms
#define CO2_FLUX_BUDGET

#undef  RELAX_NUTS
#undef  FLUX_NUTS

#undef  CHECK_CONS
#undef  DAR_DIAG_RSTAR
#undef  DAR_DIAG_DIVER
#undef  DAR_DIAG_GROW
#undef  DAR_DIAG_ACDOM
#undef  DAR_DIAG_ABSORP
#undef  DAR_DIAG_SCATTER
#undef  DAR_DIAG_PART_SCATTER
#undef  DAR_DIAG_IRR

C diagnostic chlorophyll
#undef  DAR_DIAG_CHL

C average PAR daily and store previous day
#undef  ALLOW_PAR_DAY

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C dependencies
c if two or nine species setup we don't want specific temperature ranges
#ifdef  TWO_SPECIES_SETUP
#undef TEMP_RANGE
#endif
#ifdef  NINE_SPECIES_SETUP
#undef TEMP_RANGE
#endif

c can use either denit_relax or allow_denit but not both
#ifdef ALLOW_DENIT
#undef DENIT_RELAX
#endif

#ifdef DAR_DIAG_CHL
#define ALLOW_PAR_DAY
#endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C overrides
C if you want to override dependent options, do it here:
C
C #define TEMP_RANGE

#endif /* ALLOW_DARWIN */
#endif /* DARWIN_OPTIONS_H */
