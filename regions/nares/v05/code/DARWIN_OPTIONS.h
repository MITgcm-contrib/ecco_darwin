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

C tracer selection

C enable nitrogen quotas for all plankton
#undef  DARWIN_ALLOW_NQUOTA

C enable phosphorus quotas for all plankton
#undef  DARWIN_ALLOW_PQUOTA

C enable iron quotas for all plankton
#undef  DARWIN_ALLOW_FEQUOTA

C enable silica quotas for all plankton
#undef  DARWIN_ALLOW_SIQUOTA

C enable chlorophyll quotas for all phototrophs
#define  DARWIN_ALLOW_CHLQUOTA

C enable internal c store for all phototrophs; also enables exudation
#undef  DARWIN_ALLOW_CSTORE
#undef  DARWIN_ALLOW_CSTORE_DIAGS

C enable a dynamic CDOM tracer
#define  DARWIN_ALLOW_CDOM

C enable air-sea carbon exchange and Alk and O2 tracers
#define  DARWIN_ALLOW_CARBON

C consistently use the total pH scale for carbon chemistry coefficients
#define  DARWIN_TOTALPHSCALE

C Compile Munhoven (2013) "Solvesaphe" package for pH/pCO2
C  can still select Follows et al (2006) solver in data.darwin,
C  but will use solvesaphe dissociation coefficient options.
#define  DARWIN_SOLVESAPHE

C this needs to be defined for coupling to atmospheric model:
#undef  DARWIN_USE_PLOAD

C enable RADI sediment model
C DARWIN_ALLOW_RADI itself is RETIRED -- darwin_check.F:330-337 aborts if it is
C defined, and tells you to pick v1 or v2 instead. Leave it undef.
#undef  DARWIN_ALLOW_RADI

C RADI sediment diagenesis, version 2. Enabled 2026-09-30 (user request),
C matching ecco_darwin v06/1deg code_darwin_v4r5_v2_AO.
C
C Despite the name this is NOT radiative transfer -- RADI is the sediment
C model (benthic remineralisation / burial fluxes at the sea floor). It is
C independent of pkg radtrans, which was enabled in the same pass; the two
C were kept as separate decisions for that reason.
C
C Prerequisites checked against darwin_check.F:339-362 before enabling:
C   - requires DARWIN_BOTTOM_SINK            -> already #define'd (line ~76)
C   - mutually exclusive with DARWIN_ALLOW_RADIv1 -> not defined
C   - DARWIN_ALLOW_RADI (retired) must stay undef -> it is, above
C No input files are read, and the sed_global* coefficients all have built-in
C defaults in darwin_read_params.F:396+, so no data.darwin entries are
C required. Verified the reference deployment
C (/nobackupp19/dcarrol2/v06_1deg_latest_AO/darwin3/run/data.darwin) likewise
C carries no sed_* overrides, i.e. it runs on those same defaults.
C
C NOTE: this changes the benthic boundary condition and therefore the
C nutrient/carbon budget -- it is a solution-changing switch, not diagnostics.
C
C UNPARKED AND ENABLED 2026-10-01. The 2026-09-30 blocker was real; it was
C removed by an actual 3-way merge of the six stale overrides, not by a
C hand-patch. History kept below because the failure mode is instructive.
C
C WHAT WAS BLOCKING (2026-09-30): enabling it failed to compile,
C   darwin_cons_reset.f(1163): error #6410 ... [RADIFLUXC]
C because this directory's code_darwin OVERRIDES shadow pkg/darwin, and they
C are an AO_ecosystem setup_files_v0 snapshot predating RADIv2 (added
C upstream 2025-03-03). Six overrides lacked the feature entirely, while
C darwin_cons_reset.F -- NOT overridden -- came from ckpt69k and did use
C radiFluxC. The build mixed a 2024 darwin with a 2026 one.
C
C WHY IT LOOKED UNMERGEABLE, AND WHY THAT WAS WRONG: the 09-30 note measured
C divergence against pkg HEAD (darwin_forcing.F: "462 lines absent, 43
C added") and concluded the overrides carried too much AO science to merge.
C Most of that 462 was simply a YEAR OF UPSTREAM EVOLUTION, not AO science.
C Against the actual fork point the AO delta is small:
C   file                       vs HEAD   vs fork point (= real AO delta)
C   darwin_forcing.F             505            51
C   darwin_check.F               398           201
C   darwin_read_params.F          94            10
C   DARWIN_PARAMS.h              118             4
C   darwin_diagnostics_init.F    323             0  (AO never touched it)
C   darwin_cons.F                 72            32
C The fork point is darwin3 fd20fa810 (2024-09-25), identified because
C darwin_diagnostics_init.F is byte-identical to pkg at that commit.
C *** Measure divergence from the COMMON ANCESTOR, not from HEAD. Diffing
C *** against HEAD makes every stale fork look unmergeable.
C
C THE MERGE (2026-10-01): git merge-file per file,
C   base = fd20fa810:pkg/darwin/<f>, ours = this override, theirs = ckpt69k.
C Three merged clean; three conflicted once each. All three conflicts were
C the SAME upstream feature -- "separate NO3/NO2/NH4 runoff" (11d365acf),
C which introduces IPrunoff -- plus the scavenging-default revert. pan_AO's
C AO_ecosystem runoff scheme uses PO4runoff and has NO IPrunoff at all, so
C all three were resolved to the AO side. Taking pkg there would not even
C have compiled.
C
C TWO THINGS THE MERGE GOT RIGHT ONLY AFTER CHECKING:
C  1. hillnumUptake was dropped, correctly per upstream (removed by
C     a092808e6 when macromolecular growth replaced it) -- but AO's
C     darwin_plankton.F is NOT merged and still uses it at 4 sites. It was
C     re-added by hand. *** A merge will faithfully follow upstream's
C     deletion of a symbol your UNMERGED files still depend on. ***
C  2. DARWIN_MACROMOLECULAR_GROWTH and the other features the merge pulled
C     forward are CPP-guarded and #undef, so they are inert here.
C A sweep for identifiers dropped by the merge but still referenced by the
C non-merged overrides found exactly one more, DARWIN_ALLOW_RADI, and it
C appears only in comments and its own #undef.
C
C Backups: *.bak_preRADIv2merge_20261001 for all six.
C DARWIN_FIELDS.h was re-synced separately on 09-30 (*.bak_preRADIv2_20260930).
C
C Note the trap this exposes generally: an overriding header/source FREEZES
C whatever the package later moves on from, and nothing warns you.
#define  DARWIN_ALLOW_RADIv2


C optional bits

C enable denitrification code
#undef  DARWIN_ALLOW_DENIT

C enable separate exudation of individual elements
#undef  DARWIN_ALLOW_EXUDE

C enable old virtualflux code for DIC and Alk
#undef  ALLOW_OLD_VIRTUALFLUX

C reduce nitrate uptake by iron limitation factor
C only effective with both nitrogen and iron quota
C this option was used pre-2015 in the quota model, together with Geider,
C so the iron limitation factor is also applied to alpha_I
#undef DARWIN_NITRATE_FELIMIT

C allow organic matter to sink into bottom (sedimentize)
#define  DARWIN_BOTTOM_SINK

C CDOM is in carbon units and follows POC
#define DARWIN_CDOM_UNITS_CARBON

C include code for reading nutrient runoff from files
#define  DARWIN_NUTRIENT_RUNOFF


C light

C compute average PAR in layer, assuming exponential decay
C (ignored when radtrans package is used)
#undef  DARWIN_AVPAR

C enable GEIDER light code
#define  DARWIN_ALLOW_GEIDER

C use rho instead of acclimated Chl:C for chlorophyll synthesis
#undef  DARWIN_GEIDER_RHO_SYNTH

C initialize chl as in darwin2 (with radtrans package)
#undef  DARWIN_CHL_INIT_LEGACY

C scattering coefficients are per Chlorophyll (with radtrans package)
#undef  DARWIN_SCATTER_CHL

C make diagnostics for instrinsic optical properties available
C Enabled 2026-09-30 alongside pkg oasim/sun, matching ecco_darwin v06/1deg
C code_darwin_v4r5_v2_AO/DARWIN_OPTIONS.h. Output-only (it adds absorption/
C scattering/backscatter diagnostics in darwin_light_radtrans.F), so it
C changes no prognostic field -- the point is to be able to SEE whether CDOM
C and the Arctic optics are behaving, which is the whole reason radtrans and
C the CDOM tracer are in this configuration.
#define DARWIN_DIAG_IOP


C grazing

C for quadratic grazing as in darwin2+quota
#undef  DARWIN_GRAZING_SWITCH

C compute palat from size ratios
#undef  DARWIN_ALLOMETRIC_PALAT

C turn off grazing temperature dependence
#undef  DARWIN_NOZOOTEMP


C temperature

C turn off all temperature dependence
#undef  DARWIN_NOTEMP

C select temperature version: 1, 2 or 3
#define DARWIN_TEMP_VERSION 4

C restrict phytoplankton growth to a temperature range
#undef  DARWIN_TEMP_RANGE


C iron

C restrict maximum free iron
#define DARWIN_MINFE

C enable particle scavenging code
#define  DARWIN_PART_SCAV

C enable variable iron sediment source
#define  DARWIN_IRON_SED_SOURCE_VARIABLE

C revert to old variable iron sediment source in terms of POP
#undef  DARWIN_IRON_SED_SOURCE_POP

C add iron source from hydrothermal vents
#undef  DARWIN_ALLOW_HYDROTHERMAL_VENTS


C diagnostics

C include code for per-type diagnostics
#define  DARWIN_DIAG_PERTYPE


C debugging

C turn on debugging output
#undef DARWIN_DEBUG

C compute and print global element totals
#define DARWIN_ALLOW_CONS

C value for unused traits
#define DARWIN_UNUSED 0

C fill diagnostics for most tendency terms
#define  DARWIN_DIAG_TENDENCIES


C deprecated

C base particle scavenging on POP as in darwin2
#define  DARWIN_PART_SCAV_POP


C random trait generation

C assign traits based on random numbers as in darwin2
#undef  DARWIN_RANDOM_TRAITS

C set traits for darwin2 2-species setup (requires DARWIN_RANDOM_TRAITS)
#undef  DARWIN_TWO_SPECIES_SETUP

C set traits for darwin2 9-species setup (requires DARWIN_RANDOM_TRAITS)
#undef  DARWIN_NINE_SPECIES_SETUP

C enable diazotrophy when using (requires DARWIN_RANDOM_TRAITS)
#undef  DARWIN_ALLOW_DIAZ

#endif /* ALLOW_DARWIN */
#endif /* DARWIN_OPTIONS_H */
