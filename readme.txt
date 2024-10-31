ECCO-Darwin Github Repository
git clone git@github.com:MITgcm-contrib/ecco_darwin.git

--code_util--
   /useful scripts
   
--dic--
   /ECCO-DIC set-up
   
--doc-- contains some documentation files

Simulations:

--regions--
  --bering_strait--
    /v05: Bering Strait llc270 cut-out
  --CCS--
    /llc270: California Current System llc270 cut-out for kelp model development
  --downscaling--
    /util for downscaling  
  --GoA--
    /v05: Gulf of Alaska llc270 cut-out for Takamitsu Ito
  --GoM--
    /: Gulf of Mexico lon-lat
  --GulfGuinea--
    /v05: Gulf of Guinea llc270 cut-out for COESSING  
  --LR17--
    /llc270: La Rochelle coastal Ssstem cut-out for sediment model development
  --mac_delta--
    /LatLon: Mackenzie Delta lon-lat  
    /llc270: Mackenzie Delta llc270 cut-out
    /llc4320: Mackenzie Delta llc4320 cut-out
  --Med--
    /v05: Mediteranean and Black Sea llc270 cut-out for HCMR
  --RedSea--
    /v05: Red Sea llc270 cut-out for KAUST
    /v05_coral: Red Sea llc270 cut-out including model for coral reefs
  --totten--
    /llc1080: east Antarctica totten lon-lat

--v02--
  /cs510_brix: optimized solution described in Brix et al. (2015)
  /manizza_2019: Manizza et al. (2019) simulation

--v03--
  /cs510_brix: unpublished ECCO-Darwin v3 with linear piston velocity and revised ICs
  /cs510_latest: Same as above with newer (checkpoint66o 2018/01/30) MITgcm code

--v04--
  /3deg: 3 degree verification experiment with budget scripts
  /llc270_devel: test solution with nonlinear dissolution and benthic DIC/ALK fluxes
  /llc270_JAMES_paper: Carroll et al. 2020 (JAMES) solution
  /llc270_JAMES_budget: same as above with diagnostics for closing budgets 

--v05--
  /1deg: 1 degree solution based on ECCO V4r4/V4r5
  /1deg_CDR: 1 degree solution based on ECCO V4r4/V4r5 for mCDR simulations
  /3deg: 3 degree verification experiment w/ Darwin 3 and budget scripts
  /3deg_CDR_MacroA: Kay's macroalgae (kelp) code development set-up
  /llc270: Carroll et al. 2022 (GBC) solution w/ Darwin 3
  /llc270_jra55do: same as above but with JRA55-do daily, point-source runoff
  /llc270_jra55do_mangroves: same as above but with mangrove carbon export
  /llc270_jra55do_nutrients: same as above but with DOC, DON, and DOP surface forcing
  /llc270_N2O: same as llc270 but with total and thermal N2O, argon, and thermal oxygen
  /llc270_sediment: same as llc270 above but with RADI sediment model
 
--v06--
  /1deg: 1 degree solution based on ECCO V4r5
  /llc270: v05 llc270 w/ Darwin upgrades and new ecosystem, BGC runoff, RADI sediment model, 
   pH fix and updated solver, hydrothermal vent iron forcing, and radiative transfer package

**Notes**
v02 to v04 use Darwin 1, v05 onward uses Darwin 3

**Team Github Guidelines**

When updating ECCO-Darwin code/setup, or anything that may impact ECCO-Darwin: 

1.) Check in readme.txt file(s) with updated instructions, along with any other input files that are needed.

2.) Test updated ECCO-Darwin from scratch by checking out fresh MITgcm code, 
and following the instructions line-by-line in the readme.txt file(s). 
