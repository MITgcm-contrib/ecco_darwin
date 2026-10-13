Package `BBL` is a simple bottom boundary layer scheme.  The initial
motivation is to allow dense water that forms on the continental shelf around
Antarctica (High Salinity Shelf Water) in the CS510 configuration to sink to
the bottom of the model domain and to become a source of Antarctic Bottom
Water.  The bbl package aims to address the following two limitations of
package down_slope:

(i) In `pkg/down_slope`, dense water cannot flow down-slope unless there is a
step, i.e., a change of vertical level in the bathymetry.  In pkg/bbl, dense
water can flow even on a slight incline or flat bottom.

(ii) In `pkg/down_slope`, dense water is diluted as it flows into grid cells
whose thickness depends on model configuration, typically much thicker than a
bottom boundary layer.  In pkg/bbl, dense water is contained in a thin
sub-layer and hence able to preserve its tracer properties.

Specifically, the bottommost wet grid cell of thickness

    thk = hFacC(kBot) * drF(kBot),

with properties `tracer`, and density `rho` is divided in two sub-levels:

1. A bottom boundary layer with T/S tracer properties `bbl_tracer`,
density `bbl_rho`, and thickness `bbl_eta`.

2. A residual thickness `resThk = thk - bbl_eta` with tracer properties

    resTracer = ( tracer * thk - bbl_tracer * bbl_eta ) / resThk

such that the volume integral of `bbl_tracer` and `resTracer` is consistent with
the Tracer properties of bottommost wet grid cell.

The bottom boundary layer properties `bbl_tracer` evolve as follows:

I. Vertical detrainment between BBL and residual volume of bottommost wet cell:

(i) If `bbl_rho > rho` , the T/S properties of the BBL are detrained into
the residual volume with a vertical velocity `bbl_wvel` such that:

    bbl_eta(T+deltaT) = max( 0, bbl_eta(T) - bbl_wvel * deltaT )

(ii) If `rho >= bbl_rho` or `bbl_eta(T) = 0` then set `bbl_tracer = Tracer`

Above operations do not change tracer properties of bottommost wet cell.

II. Horizontal exchange between adjacent bottom boundary layer cells when
heavy BBL water is above or at the same level as lighter BBL water.  The
strength of horizontal exchange is controlled by velocity parameter `bbl_hvel`.

(i) Replenish thickness of donor cell, `bbl_eta(d)`, if needed:

    IF bbl_tracer(d) .EQ. Tracer(d) .AND. bbl_eta(d) .LT. bbl_initEta THEN
       bbl_eta(d) = min ( bbl_initEta, Thk(d) )
    ENDIF

(ii) Heavy water flowing out of a cell lowers `bbl_eta(d)` in donor cell and
increases `bbl_eta(r)` in receiving cell.  The volume transport `dVol` out
of the donor cell is:

    dVol = min( bbl_eta(d) * rA(d) / 2,
                (thk(r) - bbl_eta(r)) * rA(r) / 2,
                d(x/y)G * bbl_eta(d) * bbl_hvel * deltaT )

(iii) Compute and accumulate tracer tendencies for donor and receiving cell.
These will be applied to the model's tracer quantities by `bbl_tendency_apply`.

    bbl_tend = dVol * (bbl_tracer(d) - resTracer(r)) / deltaT

    bbl_tendTracer(d) = bbl_tendTracer(d) - bbl_tend / rA(d) / thk(d)

    bbl_tendTracer(r) = bbl_tendTracer(r) + bbl_tend / rA(r) / thk(r)

(iv) Adjust bbl thickness and tracer properties.

    bbl_eta(d,T+deltaT) = bbl_eta(d,T) - dVol / rA(d)
    IF bbl_eta(d) .LT. 0.0001 THEN
       bbl_eta(d) = 0.0
       bbl_tracer(d) = tracer(d)
    ENDIF

    bbl_eta(r,T+deltaT) = bbl_eta(r,T) + dVol / rA(r)

    bbl_tracer(r,T+deltaT) = ( dVol * bbl_tracer(d,T) +
                               bbl_eta(r,T) * rA(r) * bbl_tracer(r,T) ) /
                             ( bbl_eta(r,T+deltaT) * rA(r) )

The horizontal exchange is computed in flux form, so that the result does
not depend on the order in which faces are visited nor on the tiling:

1. Steps I and II(i) are done first: the direction of every active face is
   decided from the state after vertical detrainment, and donor cells are
   replenished.
2. `dVol` is computed for every face from this frozen state, then limited so
   that, summed over its faces, a cell gives away at most half of its bbl
   volume and receives at most half of its residual volume.
3. Tendencies, `bbl_eta` and `bbl_tracer` of all cells are updated at once,
   with the receiving cell's `bbl_tracer` mixing all inflows.

Faces at tile edges are computed by both neighbouring tiles from identical
halo data, which requires `OLx`, `OLy` >= 2.

Lateral BBL speed (optional, `bbl_useNofSpeed = .TRUE.`): instead of the
constant `bbl_hvel`, the speed through each active face is

    u = min( bbl_hvel, bbl_nofCoeff * g' * slope / max(|f|, bbl_fMin) )

with `g' = gravity * |deltaRho| / rhoConst` the density contrast between the
heavy BBL water and the lighter water it is compared with, `slope` the
difference in BBL-top depth between the two cells divided by `dxC` (`dyC`),
and `f` averaged onto the face. `g'*slope/f` is the Nof (1983) speed of a
dense plume on a slope, which is mostly along-slope; `bbl_nofCoeff` (default
1; the Campin & Goosse 1999 scheme in MOM uses 1/3) is the fraction that
drains downslope. `bbl_hvel` becomes an upper bound and `bbl_fMin` (default 1e-5 1/s)
keeps the speed finite near the equator. Default `.FALSE.` keeps the constant
speed.

Entrainment (optional, `bbl_selectEntrain = 1`): after the lateral exchange,
a BBL that is denser than its cell and moving (speed `U` = fastest active
face) entrains residual water of the bottommost cell at

    w_E = E(Ri) * U,  E = (0.08 - 0.1 Ri) / (1 + 5 Ri) for Ri < 0.8, else 0

(Turner 1986), with bulk Richardson number `Ri = g' * bbl_eta / U^2` and `g'`
from the BBL density minus the cell density. `bbl_eta` grows by
`w_E * dTtracerLev(kBot)` (at most up to the cell thickness) and the BBL
tracers mix with the residual tracers after this step's BBL tendency. This
moves water only within the bottommost cell, so it conserves the cell's
tracer content. Detrainment at `bbl_wvel` still applies; the net exchange is
entrainment minus detrainment. Default 0: no entrainment. Use it together
with `bbl_useNofSpeed`: with a constant `bbl_hvel` of O(1 m/s), `w_E` is
O(0.1 m/s) and the BBL thickens without bound (in `global_ocean.cs32x15`
tests, `bbl_eta` ~150 m and ~500 Sv of summed BBL transport south of 60S);
`BBL_CHECK` prints a warning for this combination.

Fading out near neutral density (optional, `bbl_gpTaper`, `bbl_tauMix`): by
default a BBL that is not denser than its cell (`g' <= 0`) is mixed into the
cell at once, while Turner entrainment is largest (`E = 0.08`) as `g' -> 0+`,
so the BBL thickness jumps where it crosses neutral density. This is
harmless in forward runs but makes the model discontinuous for the adjoint.
With `bbl_gpTaper > 0` (m/s^2) entrainment is scaled by
`min(1, g'/bbl_gpTaper)` and, each step, a fraction
`(1 - min(1, max(0, g'/bbl_gpTaper))) * min(1, dt/bbl_tauMix)` of the BBL
dissolves into its cell (BBL T/S unchanged, so the cell content is
conserved). Use both together, with `bbl_tauMix` well above the time step:
in the cone test with `bbl_tauMix = 0` (dissolve within one step) round-off
asymmetries grew to O(1) in 150 steps, and 3 h (3 steps) still grew; 1 day
was stable. Defaults 0 keep the original switch (bit-identical results).
In the cone test (Nof speed + entrainment, `bbl_gpTaper = 1e-5`,
`bbl_tauMix = 1 day`) a longer-lived light BBL reduced the dye beyond the
summit from 26% to 21%; without entrainment the change was negligible.

Plume-like detrainment (optional, `bbl_selectDetrain = 1`): by default a BBL
denser than its cell thins at the constant `bbl_wvel` everywhere, also while
it descends. With `bbl_selectDetrain = 1` it detrains at
`bbl_wvel * max(0, 1 - uOut/bbl_uDet)`, where `uOut` is the fastest speed at
which the cell's BBL flows out (`bbl_uDet`, default 0.01 m/s): a stagnant BBL
detrains at `bbl_wvel`, a descending plume keeps its water and loses it only
downslope and, at neutral density, by mixing into the cell (cf. Price and
Baringer 1994). Continuous in the outflow speed when `bbl_useNofSpeed` is on.
With entrainment, combine it with `bbl_gpTaper`/`bbl_tauMix`: without them the
plume thickens until it reaches the neutral-density switch, and round-off
asymmetries in the cone test grew to 3e-2; with `bbl_gpTaper = 1e-5`,
`bbl_tauMix = 1 day` it stays mirror-symmetric (2e-13), conserves, and is
bit-identical across tilings. Effect in the cone test: more, thinner BBL cells
(756 active cells, mean 17.6 m vs 264 and 43.6 m); dye beyond the summit
20.5% vs 21.0%. Default 0 keeps the original scheme (bit-identical).

Passive tracers (`pkg/ptracers`, e.g. DIC, O2 and nutrients of a
biogeochemical model): with `bbl_doPtracers = .TRUE.` (default when
`usePTRACERS = .TRUE.`) the BBL also carries every passive tracer, with the
same exchange, mixing, entrainment and reset as T and S (hidden BBL
concentrations `bbl_ptr`, applied through `PTRACERS_APPLY_FORCING` over
`PTRACERS_dTLev`, and saved in `pickup_bbl` as `bblPt001`, ...). Without it,
dense water moved by the BBL would leave its tracers behind. A passive tracer
initialised equal to theta (or salinity) stays bit-identical to theta in a
test where pkg/bbl is the only process.

Restarting from a pickup without `pickup_bbl` (e.g. switching the BBL on in
a spun-up run): set `bbl_coldStart = .TRUE.` to initialise the BBL as at
iteration 0 (BBL T/S from the bottom cell, zero thickness, or the
`bbl_*File` fields if set) instead of reading `pickup_bbl`.

Diagnostics: `BBLtheta`, `BBLsalt`, `BBLeta` (BBL state, filled at the start
of the time step), `BBLtendT`/`BBLtendS` (BBL tendency of the bottom cell),
`BBLwEnt` (entrainment velocity, m/s), `BBLuTrp`/`BBLvTrp` (exchanged volume
per time step, m^3/s) and `BBLuVel`/`BBLvVel` (speed before limiter, m/s),
positive towards increasing i/j.

Adjoint (TAF): `bbl_calc_rhs`, `bbl_calc_rho`, `bbl_init_varia` and
`bbl_tendency_apply` are in `bbl_ad_diff.list`; the BBL state is stored at the
start of `BBL_CALC_RHS` and in the checkpoint-level directives, and
initialised after pkg/ctrl so that initial-condition controls reach the BBL.
TAF differentiates the branch taken at each switch, so use the fade-out
(`bbl_gpTaper`, `bbl_tauMix`) for smoother gradients. Gradient checks in
`global_oce_latlon` (code_ad + bbl, 1 - finite difference/adjoint): 1e-10 to
2e-4 with the constant speed, Nof speed, entrainment, fade-out and plume-like
detrainment (4 MPI ranks identical to serial; z* up to 3e-3); 8e-4 for passive
tracers (as the stock ptracers adjoint); the tangent-linear model agrees too.
In ECCO v4r6 (LLC90, checkpoint68g) a 24-step check at Antarctic shelf cells
agrees to 0.3% with GM/Redi off; with GM/Redi on it is 3-11% off because
ECCO's `GMREDI_WITH_STABLE_ADJOINT` zeroes the adjoint of the isopycnal slopes
each step, and the BBL fade-out feels the bottom-cell density that GM bolus
advection changes. The BBL state adds tape records: `nWh` in `MDSIO_BUFF_WH.h`
may need raising (727 -> 1000 in `tutorial_tracer_adjsens`).

Recommended settings (from the tests summarised below; the defaults keep
the original constant-speed scheme so that existing set-ups are unchanged):

     bbl_useNofSpeed   = .TRUE.,
     bbl_nofCoeff      = 1.,
     bbl_selectEntrain = 1,
     bbl_hvel          = 1.,     (only an upper bound on the speed)

- `global_ocean.cs32x15` with sea ice (10 years): every BBL variant made the
  deep Southern Ocean (south of 60S, below 2 km) denser; with the settings
  above, by about half as much as the default constant 1 m/s, with 4 Sv
  instead of 221 Sv of summed BBL transport. `bbl_nofCoeff = 1` gave twice
  the effect of 1/3; entrainment with a constant speed ran away (see above).
- Idealised polynya on a Southern Ocean slope (5 km, 120 days): with these
  settings ~5 times more dense-water mass reached below 1500 m than without
  the BBL (35 vs 7 Mt), and water denser than any initial water reached
  2000-2500 m (none below 1300 m without the BBL); a constant 0.1 m/s BBL
  gave only ~1.5 times.
- All variants conserve heat and salt to round-off and give bit-identical
  results across tilings, MPI/OpenMP and restarts.

Some key limitations of package as it stands:

- Detrainment velocity (`bbl_wvel`) is prescribed rather than computed, as is
  the horizontal velocity (`bbl_hvel`) unless `bbl_useNofSpeed` is set.
- bbl depth (`bbl_eta`) cannot be thicker than thickness of bottommost wet grid cell.
- The BBL density of the upslope cell is compared with the resident density of
  the downslope cell at the downslope cell's own level. When the upslope cell
  has the deeper bottom level (a thick BBL over a slightly shallower level
  boundary) the two are at different pressures; this did not occur in
  `global_oce_latlon` (0 of 3684 candidate faces).
- The entrainment speed `U` is the BBL speed relative to an ambient fluid at
  rest; the resolved bottom velocity is not used.
- The adjoint differentiates through hard switches (heavy-over-light test,
  exact-equality replenish test, `bbl_eta < 1e-4` reset, limiter and speed
  cap); finite differences that cross one see a jump the adjoint does not.
