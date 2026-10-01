# wad_balzano: Balzano (1998) tidal-flat tests: slope, step, pool

A sloping beach driven by an M2-like tide at an open boundary, after
Balzano's (1998) benchmark tests: the shoreline must follow the tide over a
plain slope, flow over a 1 m bed step without negative depths, and leave a
pool trapped behind a sill at the sill crest (plus at most `wadCritDepth`)
at low water.

## Set-up

- 140 columns of 100 m (5 tiles of 28), 1 level of 8 m, r*.
- Bed from −5 m at the boundary to 0 m at 13.8 km; land at the end.
- A 2 m, 12 h tide imposed at the western boundary (`code/obcs_calc.F`,
  `code/balzano_get_etan.F`); model rest level 2.5 m above mean sea level.
- Δt = 10 s, 12960 steps = 3 tidal cycles of 12 h; starts at high water.
- Packages: `gfd obcs wad diagnostics`.

## Variants

| Variant | Change |
|---|---|
| `step` | a 1 m bed step half way up the beach (`input.step/bathy.bin`) |
| `pool` | a 0.5 m sill and a 1 m deep pool on the upper beach (`input.pool/bathy.bin`) |
| `zstar` | surface-level non-linear free surface (`select_rStar=0`, `hFacInf=0.005`) |

## Inputs

All binary inputs are in the repository; to regenerate them: `input/gendata.py` writes `input/bathy.bin`, `input.step/bathy.bin`,
`input.pool/bathy.bin` and `input/eta_init.bin`.

## Build and run

With testreport (runs the main case and every `input.*` variant and compares
with `results/`), from MITgcm's `verification/` directory:

```sh
./testreport -t wad_balzano            # add -of <optfile>, -mpi as usual
```

By hand:

```sh
cd verification/wad_balzano
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

Specific to this experiment: `input/balzano.py RUNDIR [--plot]` reports the minimum depth, the volume
budget, the shoreline at each low and high water against a bathtub estimate,
and for `pool` the lowest pool level against the sill crest (expect crest +
~0.100 m).
