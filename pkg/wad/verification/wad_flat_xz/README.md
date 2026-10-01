# wad_flat_xz: Stratified beach at rest, and under a tide

An x-z section of a stratified ocean over a shelf and a beach whose top is
dry at rest. With no forcing every velocity is spurious: the test measures
the pressure-gradient error of thin r* columns next to dry cells (expect
max |u| below 1 mm/s, decaying). The `tide` variant drives the same section
with an M2 tide and KPP.

## Set-up

- 100 columns of 100 m (5 tiles of 20), 10 levels of 2.5 m, r*.
- 20 m shelf for x < 5 km, beach rising to +2 m at 10 km; model rest level
  3 m above sea level.
- Temperature stratified with depth only (20 °C at the surface to 10 °C at
  20 m).
- Main case: Δt = 30 s, 2880 steps = 1 day at rest.
- Packages: `gfd obcs kpp wad`.

## Variants

| Variant | Change |
|---|---|
| `tide` | 1.5 m M2 tide at a western Flather boundary, KPP on; Δt = 6 s, one tidal cycle |

## Inputs

All binary inputs are in the repository; to regenerate them: `input/gendata.py` writes `bathy.bin`, `eta.bin` and `T.bin`.

## Build and run

With testreport (runs the main case and every `input.*` variant and compares
with `results/`), from MITgcm's `verification/` directory:

```sh
./testreport -t wad_flat_xz            # add -of <optfile>, -mpi as usual
```

By hand:

```sh
cd verification/wad_flat_xz
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

Specific to this experiment: Main case: `dynstat_uvel_max/min` in `%MON` must stay below ~1 mm/s
(0.53 → 0.21 mm/s over the day, the same as without dry cells). `tide`: T
stays within its initial range.
