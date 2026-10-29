# sediment_seiche: a seiche over a mud and sand bed (no pkg/wetdry)

The standalone test of `pkg/sediment`. It needs only pkg/ptracers: no
pkg/wetdry, no open boundaries, no river. A 0.5 m surface tilt in a closed
basin starts a seiche. Its currents (up to 0.52 m/s) erode a uniform bed of
mud and fine sand, which then settles again. Nothing in the column dries.

## Set-up

- 40 × 1 columns of 500 m (2 tiles of 20 × 1), land at both ends: a
  closed basin 19 km long and 5 m deep, 10 levels of 0.5 m.
- Initial surface: a linear tilt from +0.5 m (west) to −0.5 m (east).
  Seiche period ~1.5 h (2L/√(gH)).
- Uniform T and S, linear EOS, no T/S stepping: the flow is barotropic.
- Quadratic bottom drag 2.5e-3, which also sets the bed stress.
- Δt = 20 s, 360 steps = 2 h. Linear free surface.
- Two classes (`data.sediment`): mud (pTracer 1, w_s = 0.5 mm/s,
  τ_ce = 0.05 Pa) and fine sand (pTracer 2, w_s = 10 mm/s,
  τ_ce = 0.15 Pa), 5 kg/m^2 of each on the bed at t = 0, clear water.
- Packages: `gfd ptracers sediment`.

## Variants

| Variant | Change |
|---|---|
| `rstar` | non-linear free surface with r\* (`nonlinFreeSurf=4`, `select_rStar=2`) |
| `layers` | layered bed (`SEDIMENT_PARM02`): 6 layers, porosity 0.9, mud cohesive, consolidation and swelling of τ_cb (1 h and 10 h, short for this test) |

## Inputs

All binary inputs are in the repository. To make them again, run
`python3 gendata.py` in `input/`. It writes `bathy.bin`, `eta.bin`,
`bed_mud.bin` and `bed_sand.bin`.

## Build and run

With testreport, from MITgcm's `verification/` directory:

```sh
./testreport -t sediment_seiche
```

Add `-of <optfile>` and `-mpi` as usual. By hand:

```sh
cd verification/sediment_seiche
mkdir -p build run
cd build
../../../tools/genmake2 -mods ../code -of <optfile>
make depend && make
cd ../run
ln -s ../input/* .
../build/mitgcmuv > output.txt
```

For a variant `input.X`, link its files first, then the rest of `input/`:

```sh
mkdir ../run_X && cd ../run_X
ln -s ../input.X/* .
for f in ../input/*; do [ -e $(basename $f) ] || ln -s $f .; done
../build/mitgcmuv > output.txt
```

## What to check

`%SEDIMENT_MON` lines (every `monitorFreq`, 30 min):

```
%SEDIMENT_MON: iter ... class n susp,bed [kg]= ...
%SEDIMENT_MON: iter ... class n cumRiver,cumDep,cumEro [kg]= ...
%SEDIMENT_MON: iter ... max step budget error [kg]= ...
```

- For each class, suspended + bed mass stays at 4.75e7 kg (5 kg/m^2 on
  38 wet columns of 500 m × 500 m).
- The step budget error stays at round-off (~1e-9 kg, ~1e-17 of the mass).
- Both classes erode. After 2 h, main run: mud cumulative erosion
  3.4e6 kg, sand 3.4e6 kg. `rstar` gives nearly the same values.
- `layers`: the bed armors. Mud erosion falls to 0.23e6 kg, because the
  50 % mud bed is cohesive and its bulk critical stress controls.

Compare the `%MON` lines with `results/output*.txt` (testreport does
this).

Restart check (done for `layers`, 180 + 180 steps against 360 steps): the
`%SEDIMENT_MON` lines are identical. The dynamics statistics differ at
~1e-9 relative after a restart in this set-up, and the pTracer statistics
follow them (~1e-10). The dynamics difference is the same with
`useSEDIMENT=.FALSE.`, so it does not come from pkg/sediment.

Reference outputs were made with testreport (default `-ieee`) on arm64
macOS with conda-forge gfortran (optfile `darwin_arm64_gfortran_conda`).
