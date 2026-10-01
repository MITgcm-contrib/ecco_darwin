# wad_thacker_1d: Thacker (1981) oscillating basin

A frictionless parabolic basin in 1-D (x-z section) whose planar surface
sloshes with an exact analytic solution (Thacker 1981): the shoreline moves
over the dry beach on both sides every period. It is the accuracy test of
the package: volume and salt must be exact, and the surface error against
the exact solution measures the wetting-drying scheme alone.

## Set-up

- 250 columns of 40 m (5 tiles of 50), 1 level of 20 m, r* (`select_rStar=2`).
- Basin half-width 3 km, centre depth 10 m below Thacker's mean level; the
  model rest level is 10 m above it (`DATUM`), so the whole beach is ocean.
- Δt = T/272 ≈ 4.95 s (period T = 1346 s), 1360 steps = 5 periods; surface
  written every T/8.
- Packages: `gfd wad diagnostics`.

## Variants

| Variant | Change |
|---|---|
| `carry` | `wadCarryVel=.TRUE.`: a face that opens takes the upstream depth-mean velocity |
| `zstar` | surface-level non-linear free surface (`select_rStar=0`, `hFacInf=0.001`) |

## Inputs

All binary inputs are in the repository; to regenerate them: `input/analytic.py gen` writes `bathy_thacker.bin` and `eta_thacker.bin`
(run it in `input/`).

## Build and run

With testreport (runs the main case and every `input.*` variant and compares
with `results/`), from MITgcm's `verification/` directory:

```sh
./testreport -t wad_thacker_1d            # add -of <optfile>, -mpi as usual
```

By hand:

```sh
cd verification/wad_thacker_1d
mkdir -p build run
cd build
../../../tools/genmake2 -mods ../code -of <optfile>
make depend && make
cd ../run
ln -s ../input/* .
../build/mitgcmuv > output.txt
```

For a variant `input.X`, link its files first and then the rest of `input/`:

```sh
cd ../run_X        # a fresh directory
ln -s ../input.X/* .
for f in ../input/*; do [ -e $(basename $f) ] || ln -s $f .; done
../build/mitgcmuv > output.txt
```

## What to check

Every run prints `%WAD_MON` lines (every `wadMonFreq`):

```
%WAD_MON: iter= ... nDry= ... nFlips= ... nLimit= ... minDepth= ...
%WAD_MON: vol,dVolAdd,dVolTake,cumVolAdj= ...
%WAD_MON: cumInflow,budgetErr,budgetErr/vol0= ...
```

`minDepth` must stay at or above `wadMinDepth` (0.05 m) and
`budgetErr/vol0` (the volume budget error relative to the initial volume)
at round-off, ~1e-15. Compare the `%MON` lines with `results/output*.txt`
(testreport does this).

Specific to this experiment: `input/analytic.py cmp RUNDIR [--plot] [--exact]` compares `Eta.*.data`
with the time-discrete exact solution (or the continuous one with
`--exact`): expect an L2 error of ~7 % of the initial amplitude after 5
periods (r\*), volume and salt exact.
