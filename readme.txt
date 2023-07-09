ECCO-Darwin Github Repository
git clone git@github.com:MITgcm-contrib/ecco_darwin.git

Simulations:

--doc-- contains some documentation files, including tag-index

--regions--
  --CCS--
    /llc270: California Current System llc270 cut-out for kelp model development
  --GoM--
    /: Gulf of Mexico lon-lat
  --mac_delta--
    /LatLon: Mackenzie Delta lon-lat  
    /llc270: Mackenzie Delta llc270 cut-out
    /llc4320: Mackenzie Delta llc4320 cut-out
  --totten--
    /llc1080: totten lon-lat

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
  /llc270_jra55do_nutrients: same as above but with DOC, DON, and DOP surface forcing
  /llc270_sediment: same as llc270 above but with RADI sediment model
  
--v06--
  /1deg: 1 degree solution based on ECCO V4r5
  /llc270: v05 llc270 w/ Darwin upgrades and new ecosystem, BGC runoff, RADI sediment model, 
   pH fix and updated solver, hydrothermal vent iron forcing, and radiative transfer package

--code_util--
  /LOAC/GlobalNews: biogeochimical exports from rivers (DIN, DON, DIP, DOP, DOC, PN, PP, POC, DSi) based on Mayorga et al., 2010

**Notes**
v02 to v04 use Darwin 1, v05 onward uses Darwin 3

**Team Github Guidelines**

When updating ECCO-Darwin code/setup, or anything that may impact ECCO-Darwin: 

1.) Check in readme.txt file(s) with updated instructions, along with any other input files that are needed.

2.) Test updated ECCO-Darwin from scratch by checking out fresh MITgcm code, 
and following the instructions line-by-line in the readme.txt file(s). 

3.) Once test is successful, if needed, add a full description of changes to MITgcm_contrib/ecco_darwin/doc/tag-index. 
Include the date and list all specific changes.

Example: 

05/03/19
- To make initialization and pickups more robust:
  - remove pickupSuff from data
  - remove PTRACERS_Iter0 and PTRACERS_initialFile from data.ptracers
  - delete pickup_dic.0000000001.data so that it does not get used for
    initialization because of above changes
