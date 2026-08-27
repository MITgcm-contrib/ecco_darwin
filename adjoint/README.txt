Adjoint-based experiments and analyses built on MITgcm/darwin3.

tutorial_global_oce_biogeo/  Basic adjoint-of-DIC tutorial (TAF), the
                              starting template for the experiments below.
                              See tutorial_global_oce_biogeo/readme.txt.

global_oce_biogeo_bling_SOCAT/  Adjoint (TAF) optimization of BLING
                              biogeochemical initial conditions against a
                              SOCAT surface-pCO2 climatology, via an
                              offline M1QN3 driver. See its README.md.

global_oce_biogeo_darwin/    Forward-only Darwin ecosystem model, run on
                              the same grid/physics/forcing as
                              global_oce_biogeo_bling_SOCAT. See its
                              README.md.

darwin_bling_comparison/     Pilot study comparing Darwin's forward
                              finite-difference sensitivity against
                              BLING's adjoint sensitivity of the SOCAT
                              pCO2 cost function, at matched grid points.
                              See its README.md.

hybrid_darwin_bling/         Hybrid Darwin(forward)/BLING(adjoint)
                              surrogate-gradient optimization of DIC/ALK
                              initial conditions against SOCAT, using
                              BLING's adjoint as an approximate gradient
                              for Darwin's own (independently computed)
                              cost function. See its README.md.

hybrid_darwin_dic/           The same hybrid experiment driven by
                              pkg/dic's adjoint in place of BLING's, on
                              the 1-month SOCAT cost. Records both the
                              M1QN3-driven loop, which stalled, and the
                              manual steepest-descent loop that replaced
                              it, along with the genarr3d slot-order
                              constraint on data.ctrl. See its README.md.
