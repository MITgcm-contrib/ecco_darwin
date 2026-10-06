# pkg/gmredi  (wad: ~/Documents/research/ECCO/wetting_drying/MITgcm)

Calculates the 3D diffusivity as per Bates et al. (2014) \ev

**vs its upstream base (merge-base, see README):** changed: gmredi_calc_tensor.F

**in groups:** oceanic
**runtime switch:** `useGMREDI`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.gmredi`
**manual:** `doc/phys_pkgs/gmredi.rst`, `doc/algorithm/algorithm.rst`, `doc/contributing/contributing.rst`, `doc/examples/cfc_offline/cfc_offline.rst`, `doc/examples/examples.rst`
**adjoint support files:** gmredi_ad_diff.list

## README
```
heimbach@mit.edu, 16-Aug-2001

The following files (and probably only those) could
potentially cause problems in the context of AD,
i.e. may leed to large sensitivities:

- gmredi_calc_tensor:
  (all for the case GM_VISBECK_VARIABLE_K)
         N2=(-Gravity*recip_Rhonil)/dRdSigmaLtd(i,j)
         SN=sqrt(Ssq*N2)
        VisbeckK(i,j,bi,bj)=
     &     min(VisbeckK(i,j,bi,bj),GM_Visbeck_maxval_K)

- gmredi_slope_limit:
  gradSmod, Small_Number, dSigmaDrLtd, Smod
  Lrho=Cspd/abs(Fcori(i,j,bi,bj))
```

## Namelist parameters
### GM_PARM01
- `GM_AdvForm` — use Advective Form (instead of Skew-Flux form)
- `GM_AdvSeparate` — do separately advection by Eulerian and Bolus velocity
- `GM_InMomAsStress` — apply GM as a stress in momentum Eq.
- `GM_isopycK` — Isopycnal diffusivity [m^2/s] (Redi-tensor)
- `GM_background_K` — Thickness diffusivity [m^2/s] (GM bolus transport)
- `GM_iso2dFile` — input file for 2.D horiz scaling of Isopycnal diffusivity
- `GM_iso1dFile` — input file for 1.D vert. scaling of Isopycnal diffusivity
- `GM_bol2dFile` — input file for 2.D horiz scaling of Thickness diffusivity
- `GM_bol1dFile` — input file for 1.D vert. scaling of Thickness diffusivity
- `GM_K3dRediFile` — input file for background 3.D Isopycal(Redi) diffusivity
- `GM_K3dGMFile` — input file for background 3.D Thickness (GM) diffusivity
- `GM_background_K3dFile`
- `GM_isopycK3dFile`
- `GM_taper_scheme` — select which tapering/clipping scheme to use
- `GM_maxSlope` — maximum slope (tapering/clipping) [-]
- `GM_Kmin_horiz` — minimum horizontal diffusivity [m^2/s]
- `GM_Small_Number` — epsilon used in computing the slope
- `GM_slopeSqCutoff` — slope^2 cut-off value
- `GM_Scrit` — parameter for 'dm95' & 'ldd97' tapering fct
- `GM_Sd` — parameter for 'dm95' & 'ldd97' tapering fct
- `GM_isoFac_calcK` — add fraction of  dynamically computed variable K, e.g. Visbeck, also to Redi tensor (default = 1.)
- `GM_facTrL2dz` — minimum Trans. Layer Thick. as a factor of local dz
- `GM_facTrL2ML` — maximum Trans. Layer Thick. as a factor of Mix-Layer Depth
- `GM_maxTransLay` — maximum Trans. Layer Thick. [m]
- `GM_UseBVP` — use Boundary-Value-Problem method for Bolus transport
- `GM_BVP_cMin` — minimum value for wave speed parameter "c" in BVP [m/s]
- `GM_BVP_ModeNumber` — vertical mode number used for speed "c" in BVP transport
- `GM_useSubMeso` — use parameterization of mixed layer (Sub-Mesoscale) eddies
- `subMeso_Ceff` — efficiency coefficient of Mixed-Layer Eddies [-]
- `subMeso_invTau` — inverse of mixing time-scale in sub-meso parameteriz. [s^-1]
- `subMeso_LfMin` — minimum value for length-scale "Lf" [m]
- `subMeso_Lmax` — maximum horizontal grid-scale length [m]
- `GM_Visbeck_alpha`
- `GM_Visbeck_length`
- `GM_Visbeck_depth`
- `GM_Visbeck_minDepth`
- `GM_Visbeck_maxSlope`
- `GM_Visbeck_minVal_K`
- `GM_Visbeck_maxVal_K`
- `GM_useBatesK3d` — use Bates etal (2014) calculation for 3-d K
- `GM_Bates_smooth` — Expand PV closure in terms of baroclinic modes (=.FALSE. for debugging only!)
- `GM_Bates_use_constK` — Imposes a constant K for the eddy transport
- `GM_Bates_beta_eq_0` — Ignores the beta term when calculating grad(q)
- `GM_Bates_ThickSheet` — Use a thick PV sheet
- `GM_Bates_surfK` — Imposes a constant K in the surface layer
- `GM_Bates_constRedi` — Imposes a constant K for the Redi diffusivity
- `GM_Bates_gamma` — mixing efficiency for 3D eddy diffusivity [-]
- `GM_Bates_b1` — an empirically determined constant of O(1)
- `GM_Bates_EadyMinDepth` — upper depth for Eady calculation
- `GM_Bates_EadyMaxDepth` — lower depth for Eady calculation
- `GM_Bates_Lambda`
- `GM_Bates_smallK`
- `GM_Bates_maxK` — Upper bound on the diffusivity
- `GM_Bates_constK` — Constant diffusivity to use when GM_useBatesK3d=T and GM_Bates_use_constK=T and/or GM_Bates_constRedi=T
- `GM_Bates_maxC`
- `GM_Bates_Rmax` — Length scale upper bound used for calculating urms
- `GM_Bates_Rmin` — Length scale lower bound for calc. the eddy radius
- `GM_Bates_minCori` — minimum value for f (prevents Pb near the equator)
- `GM_Bates_minN2` — minimum value for the square of the buoyancy frequency
- `GM_Bates_surfMinDepth` — minimum value for the depth of the surface layer
- `GM_Bates_vecFreq` — Frequency at which to update the baroclinic modes
- `GM_Bates_minRenorm` — minimum value for the renormalisation factor
- `GM_Bates_maxRenorm` — maximum value for the renormalisation factor
- `GM_useLeithQG` — add Leith QG viscosity to GMRedi tensor
- `GM_useGEOM` — use the GEOME formulation to calculate kgm
- `GEOM_lmbda` — lin eddy energy dissipation rate
- `GEOM_alpha` — non-dim eddy efficiency param (=<1 in QG)
- `GEOM_ini_EKE` — initial value for depth-int eddy kinetic energy
- `GEOM_diffKh_EKE` — depth-int param eddy energy diffusion coeff
- `GEOM_vert_struc` — allow for N2 structure function
- `GEOM_vert_struc_min` — lower bound on N2/Nref vertical structure func
- `GEOM_vert_struc_max` — upper bound on N2/Nref vertical structure func
- `GEOM_minval_K` — lower bound on diffusivity
- `GEOM_maxval_K` — upper bound on diffusivity
- `GM_MNC`

## CPP options (defaults as shipped)
- `GM_READ_K3D_REDI` (undef, GMREDI_OPTIONS.h) — Allows to read-in background 3-D Redi and GM diffusivity coefficients Note: need these to be defined for use as control (pkg/ctrl) parameters
- `GM_READ_K3D_GM` (undef, GMREDI_OPTIONS.h)
- `GM_VISBECK_VARIABLE_K` (undef, GMREDI_OPTIONS.h) — This allows to use Visbeck et al formulation to compute K_GM+Redi
- `GM_GEOM_VARIABLE_K` (undef, GMREDI_OPTIONS.h) — This allows to use the GEOMETRIC formulation to compute K_GM
- `GM_BATES_K3D` (undef, GMREDI_OPTIONS.h) — This allows the Bates et al formulation to calculate the bolus transport and K for Redi
- `GM_BATES_PASSIVE` (undef, GMREDI_OPTIONS.h)
- `GM_NON_UNITY_DIAGONAL` (define, GMREDI_OPTIONS.h) — This allows the leading diagonal (top two rows) to be non-unity (a feature required when tapering adiabatically).
- `GM_EXTRA_DIAGONAL` (define, GMREDI_OPTIONS.h) — Allows to use different values of K_GM and K_Redi ; also to be used with the advective form (Bolus velocity) of GM
- `GM_BOLUS_ADVEC` (define, GMREDI_OPTIONS.h) — Allows to use the advective form (Bolus velocity) of GM instead of the Skew-Flux form (=default)
- `GM_BOLUS_BVP` (define, GMREDI_OPTIONS.h) — Allows to use the Boundary-Value-Problem method to evaluate GM Bolus transport
- `ALLOW_GM_LEITH_QG` (undef, GMREDI_OPTIONS.h) — Allow QG Leith variable viscosity to be added to GMRedi coefficient
- `GM_AUTODIFF_EXCESSIVE_STORE` (undef, GMREDI_OPTIONS.h) — Related to Adjoint-code:
- `GMREDI_MASK_SLOPES` (undef, GMREDI_OPTIONS.h)

## Headers
- `GMREDI.h` — BOP
- `GMREDI_OPTIONS.h` — BOP

## Routines (29)
`gmredi_calc_bates_k.F`, `gmredi_calc_diff.F`, `gmredi_calc_eigs.F`, `gmredi_calc_geom.F`, `gmredi_calc_psi_bolus.F`, `gmredi_calc_psi_bvp.F`, `gmredi_calc_qgleith.F`, `gmredi_calc_tensor.F`, `gmredi_calc_urms.F`, `gmredi_check.F`, `gmredi_diagnostics_fill.F`, `gmredi_diagnostics_impl.F`, `gmredi_diagnostics_init.F`, `gmredi_do_exch.F`, `gmredi_init_fixed.F`, `gmredi_init_varia.F`, `gmredi_mnc_init.F`, `gmredi_output.F`, `gmredi_read_pickup.F`, `gmredi_readparms.F`, `gmredi_residual_flow.F`, `gmredi_rtransport.F`, `gmredi_slope_limit.F`, `gmredi_slope_psi.F`, `gmredi_write_pickup.F`, `gmredi_xtransport.F`, `gmredi_ytransport.F`, `submeso_calc_psi.F`

## Called from outside the package
- `GMREDI_CALC_DIFF` ← `model/src/calc_3d_diffusivity.F:204`
- `GMREDI_CALC_TENSOR` ← `model/src/do_oceanic_phys.F:1034`
- `GMREDI_CALC_TENSOR_DUMMY` ← `model/src/do_oceanic_phys.F:1040`
- `GMREDI_DO_EXCH` ← `model/src/do_oceanic_phys.F:1097`
- `GMREDI_DIAGNOSTICS_IMPL` ← `model/src/do_statevars_diags.F:84`
- `GMREDI_OUTPUT` ← `model/src/do_the_model_io.F:139`
- `GMREDI_CHECK` ← `model/src/packages_check.F:252`
- `GMREDI_INIT_FIXED` ← `model/src/packages_init_fixed.F:331`
- `GMREDI_INIT_VARIA` ← `model/src/packages_init_variables.F:272`
- `GMREDI_READPARMS` ← `model/src/packages_readparms.F:211`
- `GMREDI_WRITE_PICKUP` ← `model/src/packages_write_pickup.F:125`
- `GMREDI_RESIDUAL_FLOW` ← `model/src/thermodynamics.F:272`
- `GMREDI_RTRANSPORT` ← `pkg/generic_advdiff/gad_calc_rhs.F:626`
- `GMREDI_XTRANSPORT` ← `pkg/generic_advdiff/gad_calc_rhs.F:346`
- `GMREDI_YTRANSPORT` ← `pkg/generic_advdiff/gad_calc_rhs.F:475`

## Verification experiments compiling it (1)
`wad_estuary_3d@wad`
