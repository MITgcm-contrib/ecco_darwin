C $Header: /u/gcmpack/MITgcm/pkg/obcs/OBCS_OPTIONS.h,v 1.18 2014/09/27 00:13:40 dimitri Exp $
C $Name:  $

C CPP options file for OBCS package
C Use this file for selecting options within the OBCS package

#ifndef OBCS_OPTIONS_H
#define OBCS_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

#ifdef ALLOW_OBCS
C Package-specific Options & Macros go here

C Enable individual open boundaries
C
C 2026-10-02: WEST turned ON. It was #undef here -- a leftover from the
C square-domain era, when the model array's PERIMETER was taken to be the
C open boundary and the west/east edges of that square were found to hold
C no wet cells. pan_AO's boundary is not the perimeter: it is the DISC RIM,
C defined in grid/disc_ob_index.npz, and the rim runs through tiles in all
C four local orientations. All four directions are genuinely populated --
C grid/disc_ob_namelist.txt sets OB_Iwest, OB_Ieast, OB_Jnorth and
C OB_Jsouth, and input/data.obcs carries a full OBW* fileset (OBWt, OBWs,
C OBWu, OBWv, OBWa, OBWh, OBWsn, OBWuice, OBWvice and 33 OBWptr) exactly
C as it does for the other three.
C
C Found by dt test 25221133, which reached OBCS_CHECK and died reporting
C that the WEST flag was off while tile bi,bj= 1 1 had a Western OB --
C 31 errors in tile OB set-up.
C (That message is paraphrased deliberately. This header is #included into
C every obcs .f, so a literal CPP token quoted here would be matched by any
C grep for the real directive -- the NOTE TO EDITORS hazard, which bit this
C project on 2026-09-30 and which the first version of this comment hit:
C `grep -c` for the token in build/obcs_check.f returned 1, and it was
C this comment, not code.)
C
C Worth noting WHY this one was caught cheaply while the disc's other
C boundary traps were not: OBCS_CHECK compares the namelist's OB indices
C against the compiled CPP flags and fires loudly. The failure mode this
C project keeps hitting -- toolkit code that slices [:, :1] / [-1:, :] for
C a boundary -- produces right-shaped arrays and no error at all. A check
C that can see the disagreement is the difference.
#define ALLOW_OBCS_NORTH
#define ALLOW_OBCS_SOUTH
#define ALLOW_OBCS_EAST
#define ALLOW_OBCS_WEST

C This include hooks to the Orlanski Open Boundary Radiation code
#undef ALLOW_ORLANSKI

C Enable OB values to be prescribed via external fields that are read
C from a file
#define ALLOW_OBCS_PRESCRIBE

C Enable OB conditions following Stevens (1990)
#undef ALLOW_OBCS_STEVENS

C ===========================================================================
C ALL THREE ENABLED 2026-10-02 18:30. COMPILED IN, NOT NECESSARILY SWITCHED
C ON -- each is a runtime switch in data.obcs once the flag is defined. This
C is the same pattern already used for GMRedi/CTRL/SMOOTH/SALT_PLUME: pay the
C build cost once, keep the decision reversible without another build.
C
C The trap being closed: with these #undef, setting `useOBCSbalance=.TRUE.`
C in data.obcs is READ AND SILENTLY IGNORED. No error, no warning. The
C namelist would claim a volume constraint the binary cannot apply. That is
C the third member of a family this project has now hit repeatedly (runoff
C UNSET_RL, blank diagnostics filename) -- configuration that looks complete,
C reports nothing, and does nothing. job_build.sh now carries a pre-flight
C guard asserting that a .TRUE. in data.obcs implies the matching #define.
C
C WHY NOW -- the evidence, from 25224224 (first full model day):
C   dynstat_eta_mean -1.0022897 -> -1.0562699 over one day. Over 1.755e13 m^2
C   of wet area that is -0.054 m/day = ~-11 Sv mean, against a true expected
C   net of ~0.1 Sv (P-E-R). It does NOT settle: -14.4 Sv (step 30) decaying
C   to -7.8 (step 120) then RE-ACCELERATING to -15.2 (step 280). Sustained,
C   -19.7 m/yr -- the model drains.
C
C   This is the failure mode a circular open boundary makes unavoidable. The
C   rim is a closed loop at constant latitude (58.697-58.732 N, 1674 OB
C   cells), so a net flux imbalance has nowhere to radiate and goes straight
C   into the free surface. Until now the only volume control was
C   gen_obcs_transport_correction.py, which is NOT idempotent.
C
C BALANCE FACTOR SEMANTICS -- read obcs_balance_flow.F before changing these.
C   OBCS_balanceFac{N,S,E,W} default to 1. and SHOULD STAY THERE for a disc.
C     > 0 : boundary joins ONE global correction, inFlow/areaOB applied
C           uniformly across every participating boundary (lines ~376-380).
C     < 0 : boundary is balanced INDEPENDENTLY, its own net flow removed
C           from itself.
C     = 0 : boundary excluded from correction.
C   -1 is WRONG HERE. pan_AO's four "boundaries" are four index-space
C   segments of a single physical ring; forcing each to zero net flow
C   separately would suppress genuine sector-to-sector throughflow -- e.g.
C   Atlantic inflow in the Scotland sector balanced against Fram/CAA outflow
C   in another. The default of 1. everywhere is what makes the whole rim one
C   budget, which is what the geometry actually is.
C
C BALANCE IS A SAFETY NET, NOT A DIAGNOSIS. The correction is a UNIFORM
C velocity increment over the entire rim. If the ~15 Sv originates in one bad
C sector (a dv-extraction or transport-correction error), balance will smear
C the fix evenly around the ring and contaminate the sectors that were right.
C Enabling it does not remove the need to find where the imbalance comes
C from -- the per-boundary net transport is computable offline from the OBCS
C files and the grid metrics, and should be.
C
C SPONGE STAYS OFF AT RUNTIME for now. A sponge is a damper: switching it on
C would suppress rim noise whether or not the rim is the problem, and 25224389
C is running specifically to find out (daily-mean WVEL -- spin-up averages
C down, a stationary staircase artifact does not). The symptom that would
C have justified a sponge, advcfl_wvel 4.31, was shown to be cold-start
C ringing: it settles to 0.58-0.62 by step 180 of 288.
C
C If a sponge is later switched on, do NOT copy kerguelen's spongeThickness.
C obcs_sponge relaxes over N cells in the INDEX direction normal to each
C nominal boundary. On this staircase rim that is a clean N-cell band where
C the rim is flat in index space (`58*1` in OB_Iwest) but spans only ~0.71N
C cells normal to the rim in the 45-degree sectors -- a sponge ~30% thinner
C there, i.e. non-uniform around the ring, and a damper whose strength varies
C along a boundary can itself reflect.
C
C SEAICE sponge is included because the rim crosses the Labrador Sea and the
C Sea of Okhotsk, both seasonally ice-covered, so ice does advect across it.
C ===========================================================================

C Allow sponge layer treatment of open boundary conditions
#define ALLOW_OBCS_SPONGE

C Include hooks to sponge layer treatment of pkg/seaice variables
#define ALLOW_OBCS_SEAICE_SPONGE

C balance barotropic velocity
#define ALLOW_OBCS_BALANCE

C add tidal contributions to normal OB flow
#define ALLOW_OBCS_TIDES

C Use older implementation of obcs in seaice-dynamics
C note: most of the "experimental" options listed below have not yet
C       been implementated in new version.
#undef OBCS_UVICE_OLD

#ifdef OBCS_UVICE_OLD
C     The following five CPP options are experimental and aim to deal
C     with artifacts due to the low-frequency specification of sea-ice
C     boundary conditions compared to the model forcing frequency.
C     Ice convergence at edges can cause model to blow up.  The
C     following CPP option fixes this problem at the expense of less
C     accurate boundary conditions.
#undef OBCS_SEAICE_AVOID_CONVERGENCE

C     Smooth the component of sea-ice velocity perpendicular to the edge.
#undef OBCS_SEAICE_SMOOTH_UVICE_PERP

C     Smooth the component of sea ice velocity parallel to the edge.
#undef OBCS_SEAICE_SMOOTH_UVICE_PAR

C     Compute rather than specify seaice velocities at the edges.
#undef OBCS_SEAICE_COMPUTE_UVICE
#endif /* OBCS_UVICE_OLD */

C     Smooth the tracer sea-ice variables near the edges.
#undef OBCS_SEAICE_SMOOTH_EDGE

#endif /* ALLOW_OBCS */
#endif /* OBCS_OPTIONS_H */
