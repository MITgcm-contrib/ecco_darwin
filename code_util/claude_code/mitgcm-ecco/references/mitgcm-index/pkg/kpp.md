# pkg/kpp

K-Profile Parameterization (Large et al. 1994) vertical mixing with nonlocal transport.

**in groups:** oceanic
**runtime switch:** `useKPP`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.kpp`
**manual:** `doc/ocean_state_est/ocean_state_est.rst`, `doc/outp_pkgs/outp_pkgs.rst`, `doc/phys_pkgs/shelfice.rst`
**adjoint support files:** kpp_ad_diff.list

## Namelist parameters
### KPP_PARM01
- `kpp_freq` — Re-computation frequency for KPP parameters     (s)
- `kpp_dumpFreq` — KPP dump frequency.                             (s)
- `kpp_taveFreq`
- `KPPwriteState` — if true, write KPP state to file
- `KPP_ghatUseTotalDiffus` —  if T : Compute the non-local term using the total vertical diffusivity ; if F (=default): use KPP vertical diffusivity
- `KPPuseDoubleDiff` — if TRUE, include double diffusive contributions
- `LimitHblStable` — if TRUE (the default), limits the depth of the hbl under stable conditions.
- `KPPuseSWfrac3D`
- `minKPPhbl` — KPPhbl minimum value               (m)
- `epsln`
- `phepsi`
- `epsilon`
- `vonk`
- `dB_dz`
- `conc1`
- `conam`
- `concm`
- `conc2`
- `zetam`
- `conas`
- `concs`
- `conc3`
- `zetas`
- `Ricr`
- `cekman`
- `cmonob`
- `concv`
- `hbf`
- `zmin`
- `zmax`
- `umin`
- `umax`
- `num_v_smooth_Ri`
- `Riinfty`
- `BVSQcon`
- `difm0`
- `difs0`
- `dift0`
- `difmcon`
- `difscon`
- `diftcon`
- `Rrho0`
- `dsfmax`
- `cstar`
- `KPPmixingMaps`
- `num_v_smooth_BV`
- `num_z_smooth_sh`
- `num_m_smooth_sh`

## CPP options (defaults as shipped)
- `KPP_SMOOTH_SHSQ` (define, KPP_OPTIONS.h) — o When set, smooth shear horizontally with 121 filters
- `KPP_SMOOTH_DVSQ` (undef, KPP_OPTIONS.h)
- `KPP_SMOOTH_DBLOC` (define, KPP_OPTIONS.h) — o When set, smooth dbloc KPP variable horizontally
- `KPP_SMOOTH_DENS` (undef, KPP_OPTIONS.h) — o When set, smooth all KPP density variables horizontally
- `KPP_SMOOTH_DBLOC` (define, KPP_OPTIONS.h)
- `KPP_SMOOTH_VISC` (undef, KPP_OPTIONS.h) — o When set, smooth vertical viscosity horizontally
- `KPP_SMOOTH_DIFF` (undef, KPP_OPTIONS.h) — o When set, smooth vertical diffusivity horizontally
- `KPP_ESTIMATE_UREF` (undef, KPP_OPTIONS.h) — o Get rid of vertical resolution dependence of dVsq term by estimating a surface velocity that is independent of first level thickness in the model.
- `KPP_DO_NOT_MATCH_DIFFUSIVITIES` (undef, KPP_OPTIONS.h) — at the bottom of the mixing layer when interior mixing is noisy. This is documented somehow in van Roekel et al. (2018), 10.1029/2018MS001336 For better backward compatibility, the flags are defined as negative
- `KPP_DO_NOT_MATCH_DERIVATIVES` (undef, KPP_OPTIONS.h) — only makes sense if the diffusitivies are matched
- `KPP_SMOOTH_REGULARISATION` (undef, KPP_OPTIONS.h) — o Include/exclude smooth regularization at the cost of changed results. With this flag defined, some MAX(var,phepsi) are replaced by var+phepsi
- `KPP_SCALE_SHEARMIXING` (undef, KPP_OPTIONS.h) — o reduce shear mxing by shsq**2/(shsq**2+1e-16) according to Polzin (1996), JPO, 1409-1425), so that there will be no shear mixing with very small shear
- `KPP_GHAT` (define, KPP_OPTIONS.h) — o Include/exclude KPP non/local transport terms
- `EXCLUDE_KPP_SHEAR_MIX` (undef, KPP_OPTIONS.h) — o Exclude Interior shear instability mixing
- `EXCLUDE_KPP_DOUBLEDIFF` (undef, KPP_OPTIONS.h) — o Exclude double diffusive mixing in the interior
- `KPP_AUTODIFF_EXCESSIVE_STORE` (undef, KPP_OPTIONS.h) — o Avoid as many as possible AD recomputations usually not necessary, but useful for testing
- `ALLOW_KPP_VERTICALLY_SMOOTH` (undef, KPP_OPTIONS.h) — o Vertically smooth Ri (for interior shear mixing)

## Headers
- `KPP.h` — BOP
- `KPP_OPTIONS.h` — CPP options file for KPP package. Use this file for selecting options within the KPP package.
- `KPP_PARAMS.h` — Basic parameter header for KPP vertical mixing parameterization.  These parameters are initialized by and/or read in from data.kpp file.

## Routines (27)
`kpp_calc.F`, `kpp_calc_diff_ptr.F`, `kpp_calc_diff_s.F`, `kpp_calc_diff_t.F`, `kpp_calc_visc.F`, `kpp_check.F`, `kpp_diagnostics_init.F`, `kpp_do_exch.F`, `kpp_forcing_surf.F`, `kpp_init_fixed.F`, `kpp_init_varia.F`, `kpp_output.F`, `kpp_readparms.F`, `kpp_routines.F`, `kpp_transport_ptr.F`, `kpp_transport_s.F`, `kpp_transport_t.F`

## Called from outside the package
- `KPP_CALC_DIFF_PTR` ← `model/src/calc_3d_diffusivity.F:188`
- `KPP_CALC_DIFF_S` ← `model/src/calc_3d_diffusivity.F:181`
- `KPP_CALC_DIFF_T` ← `model/src/calc_3d_diffusivity.F:176`
- `KPP_CALC_VISC` ← `model/src/calc_viscosity.F:78`
- `KPP_CALC` ← `model/src/do_oceanic_phys.F:956`
- `KPP_CALC_DUMMY` ← `model/src/do_oceanic_phys.F:961`
- `KPP_DO_EXCH` ← `model/src/do_oceanic_phys.F:1103`
- `KPP_OUTPUT` ← `model/src/do_the_model_io.F:145`
- `KPP_CHECK` ← `model/src/packages_check.F:246`
- `KPP_INIT_FIXED` ← `model/src/packages_init_fixed.F:321`
- `KPP_INIT_VARIA` ← `model/src/packages_init_variables.F:255`
- `KPP_READPARMS` ← `model/src/packages_readparms.F:206`
- `KPP_TRANSPORT_PTR` ← `pkg/generic_advdiff/gad_calc_rhs.F:676`
- `KPP_TRANSPORT_S` ← `pkg/generic_advdiff/gad_calc_rhs.F:670`
- `KPP_TRANSPORT_T` ← `pkg/generic_advdiff/gad_calc_rhs.F:665`

## Verification experiments compiling it (5)
`1D_ocean_ice_column` `lab_sea` `seaice_obcs` `tutorial_tracer_adjsens` `vermix`
