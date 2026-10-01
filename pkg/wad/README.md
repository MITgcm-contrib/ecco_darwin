# pkg/wad — wetting and drying for MITgcm

`pkg/wad` lets MITgcm water columns dry to a thin film and re-wet: tidal
flats, river deltas, marshes and beaches can be modelled with the
non-linear free surface (r\* or the surface-level variant), with volume and
tracers conserved to round-off. It works with KPP, GM/Redi, pTracers,
OBCS, river sources (`addMass`) and the real fresh-water flux. This version
does not support sea ice (pkg/seaice); coupling to sea ice (grounded and
landfast ice on drying flats) is being developed and will be added later.

Author: Dustin Carroll (SJSU / JPL), 2026. Developed against MITgcm
`master` at commit `d861cd501` ("Fix pkg/cheapaml flux option LANL",
#1001).

## Contents

```
wad/
  README.md            this file
  src/                 the package itself (pkg/wad): 16 .F files, WAD.h,
                       WAD_OPTIONS.h
  mitgcm_changes/      MITgcm files changed to call the package
    wad_hooks.patch    the changes as one patch against d861cd501
    model/ pkg/        the same 16 files, complete, at their MITgcm paths
  verification/        5 test experiments, 16 runs with their variants
                       (code/, input*/, results/)
  doc/
    wad_formulation.md the method, step by step, with code excerpts
    figures/           figures used there
```

## Install

From the top of an MITgcm checkout (tested with `d861cd501`; the patch
touches 16 files and usually applies to nearby versions too):

```sh
cp -r <path>/wad/src  pkg/wad
git apply <path>/wad/mitgcm_changes/wad_hooks.patch      # or: patch -p1 < ...
cp -r <path>/wad/verification/wad_* verification/
```

If the patch does not apply, copy the files from `mitgcm_changes/model` and
`mitgcm_changes/pkg` instead (they are complete files at d861cd501) and
merge by hand. Then add `wad` to your experiment's `code/packages.conf`,
`useWAD=.TRUE.` to `data.pkg`, and a `data.wad` (namelist `WAD_PARM01`).

Run the tests with testreport, e.g.

```sh
cd verification
./testreport -t "wad_thacker_1d wad_balzano wad_flat_xz wad_estuary_3d wad_mudflat"
```

With `pkg/wad` compiled in and `useWAD=.FALSE.`, model results are
unchanged: advect_xz, internal_wave, global_ocean.cs32x15 and seaice_obcs
give the same testreport results as stock MITgcm (on the test machine,
internal_wave and global_ocean.cs32x15.in_p/.seaice fail in the stock build
too).

## Required model settings

`wad_check` stops the run unless:

- `#define NONLIN_FRSURF` (CPP_OPTIONS.h), `nonlinFreeSurf >= 3`,
  `select_rStar >= 0` (r\* or surface-level NLFS), no sigma coordinates,
  height (z) coordinates
- `exactConserv=.TRUE.`, `implicitDiffusion=.TRUE.`, `implicDiv2DFlow=1`
- `wadCritDepth > wadMinDepth > 0`
- not used with: pkg/seaice, pkg/thsice, pkg/salt_plume, shelfice
  remeshing; `applyExchUV_early=.FALSE.`

and warns about: flux-form momentum (use `vectorInvariantMomentum`), no
`upwindShear` with Nr > 1, explicit quadratic bottom drag with Nr > 1 (use
`selectImplicitDrag=2` with `#define ALLOW_SOLVE4_PS_AND_DRAG`), and
Adams-Bashforth tracer stepping (use a two-level scheme, e.g. 33 or 77).

A working set (verification/wad_mudflat):

```
 nonlinFreeSurf=4, select_rStar=2, exactConserv=.TRUE.,
 implicitFreeSurface=.TRUE., hFacMin=0.01, hFacMinDr=0.01, hFacInf=0.001,
 vectorInvariantMomentum=.TRUE., upwindShear=.TRUE.,
 bottomDragQuadratic=2.5E-3, selectImplicitDrag=2,
 implicitDiffusion=.TRUE., implicitViscosity=.TRUE.,
 saltAdvScheme=33, staggerTimeStep=.TRUE.,
```

## Parameters (`data.wad`, namelist `WAD_PARM01`)

| Parameter | Default | Meaning |
|---|---|---|
| `wadMinDepth` | 0.05 m | water film kept in a dry column |
| `wadCritDepth` | 0.10 m | a face is open while the water over its crest (higher surface minus higher bed) exceeds this |
| `wadAdvDepth` | 0.5 m | momentum advection ramped off below this face flow depth |
| `wadGMDepth` | 2 m | GM/Redi tensor tapered off below this column depth |
| `wadKPPCapDepth` | 20 m | KPP non-local flux fraction capped at 1 in shallower columns |
| `wadKPPDepth` | 0 (off) | optional ramp of the KPP non-local flux |
| `wadDryForcing` | T | no surface heat, salt, short-wave or tracer forcing on dry columns |
| `wadUpwindFace` | T | open-face thickness from the upwind surface |
| `wadCarryVel` | F | a face that opens takes the upstream depth-mean u\* (needs `momImplVertAdv` if Nr > 1) |
| `wadDragDepth` | 0 (off) | velocity relaxation in the last centimetres (Oey 2005) |
| `wadManningN` | 0 (off) | depth-dependent drag Cd = g n² / D^(1/3) (needs `ALLOW_BOTTOMDRAG_ROUGHNESS`) |
| `wadDragMax` | 0.05 | upper limit of that drag coefficient |
| `wadManningFile`, `wadManningT1/T2` | none | extra 2-D roughness, decaying linearly between two times |
| `wadMaxFroude` | 0 (off) | cap \|u\| ≤ Fr·√(gD) on open faces (thin sheets off banks) |
| `wadMaxSpeed` | 30 m/s | stop, and report where, above this speed |
| `wadConserveVol` | T | take the floor's added volume from neighbours |
| `wadMonFreq` | monitorFreq | `%WAD_MON` interval |
| `wadSmoothWidth` | 0 | reserved (smooth mask for the adjoint; not implemented) |

## Diagnostics and monitor

`WADdepth` (water depth), `WADdryC` (1 at the film), `WADmaskW/S` (face
masks), `WADEVLIM` (evaporation removed on dry columns). `%WAD_MON` prints
the dry-column count, face flips, limited columns, minimum depth, volume and
the volume budget error (`budgetErr/vol0`, ~1e-15 in all tests). Pickup
files `pickup_wad` and `pickup_wadface` give
bit-identical restarts.

## Verification experiments

Each experiment folder has its own `README.md`: what it tests, its set-up
and variants, how to build and run it (testreport or by hand), how to
regenerate its inputs, and what to check in the output.

| Experiment | Variants | What it tests |
|---|---|---|
| `wad_thacker_1d` | carry, zstar | Thacker (1981) oscillating basin with an exact moving shoreline |
| `wad_balzano` | pool, step, zstar | Balzano (1998) slope, bed step and trapped pool |
| `wad_flat_xz` | tide | stratified beach at rest (spurious currents) and under an M2 tide |
| `wad_estuary_3d` | carry, ptr, seiche | 3-D tidal estuary with drying flats, KPP, GM/Redi; passive tracers; a cross-estuary seiche over the flats |
| `wad_mudflat` | rain, rfwf0 | macrotidal mudflat with creeks and a river; real fresh-water flux |

All conserve volume to ~1e-15 of the total and keep uniform tracers
uniform to round-off.

## Cost and time step

With nothing drying, `useWAD=.TRUE.` adds ~9 % per time step (outflow
limiter 3 %, masks 2 %, the rest small) and does not lower the stable time
step. Drying itself can: in the 100 m estuary test the stable step fell
from 10 s to 6 s (limited at a 6 m bed step in thin bottom cells).

## Limitations

- No adjoint (the face mask is a step); `wadSmoothWidth` is reserved.
- The error against exact solutions is set by the film and the face
  threshold, not by grid size.
- Opening faces launch small waves at the shoreline (up to ~0.3 m of
  step-to-step ripple within 2–3 cells of the shoreline in the frictionless
  Thacker test; halved by bottom drag).
- The surface-level NLFS variant (`select_rStar=0`) only dries columns one
  level deep.

## References

- Adcroft, A., & Campin, J.-M. (2004). Rescaled height coordinates for accurate representation of free-surface flows in ocean circulation models. Ocean Modelling, 7, 269–284.
- Balzano, A. (1998). Evaluation of methods for numerical simulation of wetting and drying in shallow water flow models. Coastal Engineering, 34, 83–107.
- Casulli, V. (2009). A high-resolution wetting and drying algorithm for free-surface hydrodynamics. Int. J. Numer. Meth. Fluids, 60, 391–408.
- Casulli, V., & Cheng, R. T. (1992). Semi-implicit finite difference methods for three-dimensional shallow water flow. Int. J. Numer. Meth. Fluids, 15, 629–648.
- Medeiros, S. C., & Hagen, S. C. (2013). Review of wetting and drying algorithms for numerical tidal flow models. Int. J. Numer. Meth. Fluids, 71, 473–487.
- Oey, L.-Y. (2005). A wetting and drying scheme for POM. Ocean Modelling, 9, 133–150.
- Stelling, G. S., & Duinmeijer, S. P. A. (2003). A staggered conservative scheme for every Froude number in rapidly varied shallow water flows. Int. J. Numer. Meth. Fluids, 43, 1329–1354.
- Thacker, W. C. (1981). Some exact solutions to the nonlinear shallow-water wave equations. J. Fluid Mech., 107, 499–508.
- Warner, J. C., Defne, Z., Haas, K., & Arango, H. G. (2013). A wetting and drying scheme for ROMS. Computers & Geosciences, 58, 54–61.
