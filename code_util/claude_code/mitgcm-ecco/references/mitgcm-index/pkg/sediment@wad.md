# pkg/sediment  (wad: ~/Documents/research/ECCO/wetting_drying/MITgcm)

---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+---- BOP

**vs its upstream base (merge-base, see README):** new (not in its upstream base)

**runtime switch:** `useSEDIMENT`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.sediment`

## Namelist parameters
### SEDIMENT_PARM01
- `sedNcl` — number of sediment classes in use (<= sedMaxCl)
- `sedPtrIdx` — pTracer index holding each class (kg/m^3)
- `sedWs` — settling velocity [m/s] (> 0 downward)
- `sedTauCE` — critical bed shear stress for erosion [Pa]
- `sedEroM` — erosion rate constant [kg/m^2/s]: E = sedEroM*(tau_b/sedTauCE - 1) when tau_b > sedTauCE
- `sedTauCD` — critical stress for deposition [Pa] (Krone): deposited fraction of the settling flux out of the bottom layer is 1 - tau_b/sedTauCD ; <= 0: always 1
- `sedRiverC` — concentration of the added fluid (addMass) [kg/m^3]
- `sedRiverCFile` — optional 2-D map of the concentration of the added fluid [kg/m^3] (e.g. one value per river at its source cell); replaces sedRiverC where given
- `sedBulkRho` — dry bulk density of the deposit [kg/m^3] (only for the equivalent-thickness diagnostic and monitor)
- `sedBedFile` — initial bed mass file [kg/m^2] (' ': clean bed)
- `sedDryDepth` — no erosion from columns shallower than this [m]; deposition ramps from 1 at 2*sedDryDepth to 0 here
- `sedDryDeposit` — .TRUE.: deposit from drying columns (films) too
- `sedRefDepth` — the bed stress uses the mean velocity of the lowest sedRefDepth [m] of the column (not just the bottom cell, which can be a thin partial cell damped by drag); <= 0: bottom wet level only
- `sedMonFreq` — frequency of the %SEDIMENT_MON monitor [s]
- `sedConsTime` — fluff layer: deposits enter an easily eroded fluff layer that consolidates into the bed with this e-folding time [s]; <= 0 (default): no fluff layer
- `sedTauCEf` — critical stress for erosion of the fluff [Pa]
- `sedEroMf` — erosion rate constant of the fluff [kg/m^2/s]; the consolidated bed only erodes once the fluff of the column is gone
### SEDIMENT_PARM02
- `sedNlay` — number of bed layers; 0 (default): one bed store per class (with the optional fluff layer above). Layer 1 is the active layer, the only one that erodes
- `sedCohesive` — .TRUE. for cohesive (mud) classes
- `sedPorosity` — bed porosity phi: E_i = sedEroM*(1-phi)*f_i *(tau_b/tau_ce,i - 1), f_i the mass fraction of class i in the active layer (Moriarty et al. Eq. 2)
- `sedFcLo` — cohesive mass fraction of the active layer at which the bed starts / ends turning cohesive: the cohesive-behaviour factor Pc goes linearly 0 -> 1, tau_ce,i = max( Pc*tau_cb + (1-Pc)*sedTauCE(i), sedTauCE(i) ) (Sherwood et al. Eq. 6); sedTauCE is
- `sedFcHi` — cohesive mass fraction of the active layer at which the bed starts / ends turning cohesive: the cohesive-behaviour factor Pc goes linearly 0 -> 1, tau_ce,i = max( Pc*tau_cb + (1-Pc)*sedTauCE(i), sedTauCE(i) ) (Sherwood et al. Eq. 6); sedTauCE is
- `sedTcbNew` — bulk critical stress of fresh deposit [Pa]
- `sedTcbOffset` — equilibrium bulk critical stress tau_cb,eq = exp( (ln z_m - offset)/slope ) [Pa] at mass depth z_m [kg/m^2] (Sherwood et al. Eq. 4), limited to sedTcbMin .. sedTcbMax
- `sedTcbSlope` — equilibrium bulk critical stress tau_cb,eq = exp( (ln z_m - offset)/slope ) [Pa] at mass depth z_m [kg/m^2] (Sherwood et al. Eq. 4), limited to sedTcbMin .. sedTcbMax
- `sedTcbMin`
- `sedTcbMax`
- `sedTcons` — time scales [s] on which tau_cb relaxes toward tau_cb,eq from below (consolidation) and from above (swelling)
- `sedTswell` — time scales [s] on which tau_cb relaxes toward tau_cb,eq from below (consolidation) and from above (swelling)
- `sedZaCoef` — active layer thickness za = max( sedZaMin, sedZaCoef*(tau_b - tau_ce) ) [m, m/Pa] (Harris and Wiberg 1997)
- `sedZaMin` — active layer thickness za = max( sedZaMin, sedZaCoef*(tau_b - tau_ce) ) [m, m/Pa] (Harris and Wiberg 1997)
- `sedLayMass` — mass per unit area [kg/m^2] above which material moved below the active layer starts a new layer
### SEDIMENT_PARM03
- `sedWaveMode` — 0 (default): no waves; 1: bottom orbital velocity [m/s], period [s] and (optional) direction [deg, toward, counter-clockwise from +x] from files, e.g. WAVEWATCH III UBR and a mean period; 2: parametric wind sea (Young and Verhagen 1996) from sedWindSpeed,
- `sedWaveUbFile` — files of mode 1,
- `sedWavePerFile` — files of mode 1,
- `sedWaveDirFile` — files of mode 1,
- `sedWaveHsFile` — significant wave height file of mode 3 (with sedWavePerFile and sedWaveDirFile), sedWaveNrec records every sedWaveDt [s], cyclic, linear in time between records
- `sedWaveNrec`
- `sedWaveDt`
- `sedWindSpeed` — 10 m wind speed [m/s] (mode 2)
- `sedWindDir` — direction the wind blows toward [deg ccw from +x]
- `sedFetch` — fetch [m] (mode 2)
- `sedWaveGamma` — depth-limited breaking, Hs <= sedWaveGamma*depth
- `sedZ0` — bed roughness length for the wave friction factor f_w = 1.39*(A/sedZ0)**-0.52 [m] (Soulsby 1997)
- `sedWaveIce` — .TRUE.: orbital velocity times (1 - sea-ice area)

## Headers
- `SEDIMENT.h` — BOP
- `SEDIMENT_OPTIONS.h` — BOP

## Routines (13)
`sediment_add_tendency.F`, `sediment_bed.F`, `sediment_check.F`, `sediment_diagnostics_init.F`, `sediment_init_fixed.F`, `sediment_init_varia.F`, `sediment_read_pickup.F`, `sediment_readparms.F`, `sediment_step.F`, `sediment_waves.F`, `sediment_write_pickup.F`

## Called from outside the package
- `SEDIMENT_STEP` ← `model/src/forward_step.F:1133`
- `SEDIMENT_CHECK` ← `model/src/packages_check.F:514`
- `SEDIMENT_INIT_FIXED` ← `model/src/packages_init_fixed.F:655`
- `SEDIMENT_INIT_VARIA` ← `model/src/packages_init_variables.F:556`
- `SEDIMENT_READPARMS` ← `model/src/packages_readparms.F:428`
- `SEDIMENT_WRITE_PICKUP` ← `model/src/packages_write_pickup.F:259`
- `SEDIMENT_ADD_TENDENCY` ← `pkg/ptracers/ptracers_apply_forcing.F:83`

## Verification experiments compiling it (2)
`wad_mangrove@wad` `wad_mudflat@wad`
