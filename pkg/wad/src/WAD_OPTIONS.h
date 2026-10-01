#ifndef WAD_OPTIONS_H
#define WAD_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

CBOP
C !ROUTINE: WAD_OPTIONS.h
C !INTERFACE:
C #include "WAD_OPTIONS.h"

C !DESCRIPTION:
C *==================================================================*
C | CPP options file for pkg "wad" (wetting and drying):
C | Control which optional features to compile in this package code.
C *==================================================================*
CEOP

#ifdef ALLOW_WAD

C Print extra per-tile information on face flips and clipped cells
#undef WAD_DEBUG

C Smooth (tanh) face mask for adjoint differentiability (Phase 3,
C  not implemented yet: wad_check stops if wadSmoothWidth > 0)
#undef WAD_SMOOTH_MASK

#endif /* ALLOW_WAD */
#endif /* WAD_OPTIONS_H */
