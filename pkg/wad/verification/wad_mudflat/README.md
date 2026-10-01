# wad_mudflat: Macrotidal mudflat with tidal creeks and a river

An idealized macrotidal mudflat (2 m M2 tide) cut by a network of tidal
creeks, with a 20 m³/s river at the head of the main creek and passive
tracers marking river and flat water. Most of the flat dries every tide.
The variants add the real fresh-water flux.

## Set-up

- 100 × 60 columns of 100 m (4 × 2 tiles of 25 × 30), 5 levels of 2.5 m, r*.
- Offshore basin 10 m deep (x < 2 km), slope to −3 m, flat rising to +2 m
  at 10 km (~1:1400); a main creek along y = 3 km with 6 side branches.
- Tide at the western boundary (Flather, `code/obcs_calc.F`); river as
  `addMass` (`selectAddFluid=1`, fresh); river tracer relaxed by pkg/rbcs.
- Δt = 6 s, 7452 steps = 1 M2 cycle; implicit bottom drag
  (`selectImplicitDrag=2`).
- Packages: `gfd obcs kpp ptracers rbcs diagnostics wad`.

## Variants

| Variant | Change |
|---|---|
| `rain` | `useRealFreshWaterFlux` with extreme surface fresh water (`empmr.bin`: 20 mm/h evaporation, a 50 mm/h storm in hours 5–7): tests the evaporation limiter on the drying flats |
| `rfwf0` | `useRealFreshWaterFlux` with zero fresh-water flux: must reproduce the main case |

## Inputs

All binary inputs are in the repository; to regenerate them: `input/gendata.py` writes the grid, initial fields, river (`addMass.bin`),
tracers and rbcs files in `input/`.

## Build and run

With testreport (runs the main case and every `input.*` variant and compares
with `results/`), from MITgcm's `verification/` directory:

```sh
./testreport -t wad_mudflat            # add -of <optfile>, -mpi as usual
```

By hand:

```sh
cd verification/wad_mudflat
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

Specific to this experiment: Tracers stay within their initial ranges; `rain`: salinity
bounded, budget at round-off with the fresh-water flux counted;
`rfwf0`: same `%MON` as the main case.
