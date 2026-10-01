# wad_estuary_3d: 3-D tidal estuary with drying flats

A 3-D estuary: a 10 m deep channel shoaling to 3 m at its head, with tidal
flats on both sides that dry at low water, stratified salinity, Coriolis,
KPP and GM/Redi, forced by a 1.5 m M2 tide through a Flather boundary. The
variants add passive tracers, the carry option and a cross-estuary
seiche.

## Set-up

- 40 × 30 columns of 200 m (2 × 2 tiles of 20 × 15), 5 levels of 2.6 m, r*.
- Channel along x (|y − 3 km| < 600 m); flats rising to +1.3 m at the side
  walls; no flats in the first 1 km so the open boundary is in deep water.
- Salinity 30 at the surface to 34 at 10 m; f = 1e-4 s⁻¹; model rest level
  3 m above sea level; starts at high water, at rest.
- Δt = 15 s (the Flather boundary needs c Δt/Δx < 1), 2981 steps = 1 M2
  cycle.
- Packages: `gfd obcs kpp gmredi ptracers diagnostics wad`
  (`code/obcs_calc.F`: Flather tide; `code/obcs_wad_offshore.F`: offshore
  water entering on the flood).

## Variants

| Variant | Change |
|---|---|
| `carry` | `wadCarryVel=.TRUE.` with `momImplVertAdv=.TRUE.` |
| `ptr` | two passive tracers: a uniform one (must stay uniform) and a dye patch |
| `seiche` | the high-water surface tilted by 1e-4 across the estuary (±0.29 m at the walls) and released: a seiche over the drying flats; 480 steps (2 h), monitor every 10 min |

## Inputs

All binary inputs are in the repository; to regenerate them: `input/gendata.py` writes `input/` (`bathy.bin`, `eta.bin`, `S.bin`),
`input.ptr/` (`ptr1.bin`, `ptr2.bin`) and `input.seiche/eta_seiche.bin`. `WAD_NX`, `WAD_NY`,
`WAD_NR` refine the grid for runs outside testreport.

## Build and run

With testreport (runs the main case and every `input.*` variant and compares
with `results/`), from MITgcm's `verification/` directory:

```sh
./testreport -t wad_estuary_3d            # add -of <optfile>, -mpi as usual
```

By hand:

```sh
cd verification/wad_estuary_3d
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

Specific to this experiment: Salinity stays within its initial range (GM/Redi aside: ~−0.05 psu);
`ptr`: tracer 1 stays 1 to ~2e-13;
`seiche`: the cross-channel tilt sloshes and decays over the 2 hours.
