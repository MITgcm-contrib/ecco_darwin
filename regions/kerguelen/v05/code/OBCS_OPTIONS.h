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
#define ALLOW_OBCS_NORTH
#define ALLOW_OBCS_SOUTH
#define ALLOW_OBCS_EAST
C Kerguelen has a 4th (west) open boundary, unlike GoM_1km (this
C project's base template, land-bounded to the west) -- must be defined,
C confirmed via a real OBCS_CHECK fatal error (job 24947876) when left
C #undef with data.obcs's OB_Iwest set.
#define ALLOW_OBCS_WEST

C This include hooks to the Orlanski Open Boundary Radiation code
#undef ALLOW_ORLANSKI

C Enable OB values to be prescribed via external fields that are read
C from a file
#define ALLOW_OBCS_PRESCRIBE

C Enable OB conditions following Stevens (1990)
#undef ALLOW_OBCS_STEVENS

C Allow sponge layer treatment of open boundary conditions
C 2026-08-26: enabled to damp the near-boundary U/V(/T/S) flicker behind
C the FeT "reflection" artifact seen at the east (and north) open
C boundaries -- OBCS_u1_adv_Tr's per-timestep outflow/inflow test is
C driven by the sign of the local normal velocity, which can flicker at
C the edge in an eddying flow; relaxing U/V toward the prescribed OB
C profile over a short sponge zone should stabilize that sign test.
C See DECISIONS.md and the diagnostic behind this (check_boundary_reflection.py).
#define ALLOW_OBCS_SPONGE

C Include hooks to sponge layer treatment of pkg/seaice variables
#undef ALLOW_OBCS_SEAICE_SPONGE

C balance barotropic velocity
C 2026-08-20: enabled -- see DECISIONS.md. Even after fixing the
C boundary-length overcount in gen_obcs_transport_correction.py, this
C is cheap insurance against secular ETAN drift from any residual net
C imbalance across the 4 open boundaries over a ~28-day integration.
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
