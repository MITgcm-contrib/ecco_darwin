C malibu_v01 [K2]: kelp blade optics for radtrans (data.kelpoptics)
C  kelp_optics_mode :: 0 = non-spectral absorption (mp_spp_katten B Q),
C                      1 = blades as absorbing/scattering surfaces
C  kelp_proj        :: mean projection of blade area on a plane (0.5 for
C                      randomly oriented blades)
C  kelp_sfa         :: blade area per dry weight [m2/gDW]
C  kelp_absBlade, kelp_refBlade, kelp_trnBlade :: blade absorptance,
C                      reflectance, transmittance per waveband (A+R+T = 1)
#ifdef ALLOW_RADTRANS
      INTEGER kelp_optics_mode
      _RL kelp_proj, kelp_sfa
      _RL kelp_absBlade(nlam)
      _RL kelp_refBlade(nlam)
      _RL kelp_trnBlade(nlam)
      COMMON /DARWIN_KELP_OPTICS_I/ kelp_optics_mode
      COMMON /DARWIN_KELP_OPTICS_R/ kelp_proj, kelp_sfa,
     &   kelp_absBlade, kelp_refBlade, kelp_trnBlade
C malibu_v01 [K4]: TSS (mineral particle) optics, mass-specific [m2/g]
      INTEGER tss_optics_mode
      _RL tss_aStar(nlam)
      _RL tss_bStar(nlam)
      _RL tss_bbRatio
      COMMON /DARWIN_TSS_OPTICS_I/ tss_optics_mode
      COMMON /DARWIN_TSS_OPTICS_R/ tss_aStar, tss_bStar, tss_bbRatio
#endif
C malibu_v01 [K6]: particle backscatter slope, bbp ~ (l/bbp_lambda0)^-bbp_eta
C malibu_v01 [K7]: river CDOM as a fraction of DOP runoff [mmol P/mmol P]
      _RL bbp_eta, bbp_lambda0, cdom_riv_frac
      COMMON /DARWIN_BBP_CDOM_R/ bbp_eta, bbp_lambda0, cdom_riv_frac
