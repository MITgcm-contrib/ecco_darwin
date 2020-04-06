--regions--
  --mac_delta--
    /mac_delta_llc4320: mackenzie delta llc4320 cut-out
    /mac_delta_llc270: mackenzie delta llc270 cut-out
  --totten--

--v02--
  /v2_cs510_brix: optimized solution described in Brix et al. (2015)

--v03--
  /v3_cs510_Brix: unpublished ECCO-Darwin v3 with linear piston velocity and revised ICs
  /v3_cs510_latest: Same as above with newer (checkpoint66o 2018/01/30) MITgcm code

--v04--
  /v4_3deg: 3 degree verification experiment
  /v4_llc270_devel: test solution with nonlinear dissolution and benthic DIC/ALK fluxes
  /v4_llc270_JAMES_paper: Carroll et al. 2020 solution
  /v4_llc270_JAMES_budget: same as above with diagnostics for closing budgets 

--v05--
  /v5_llc270_jra55do: Carroll et al. 2020 solution with daily jra55 DO runoff


**Notes**
v02 to v04 use Darwin 1, v05 uses Darwin 3

**Team Github Guidelines**

When updating ECCO-Darwin code/setup, or anything that may impact ECCO-Darwin: 

1.) Check in readme.txt file(s) with updated instructions, along with any other input files that are needed.

2.) Test updated ECCO-Darwin from scratch by checking out fresh MITgcm code, 
and following the instructions line-by-line in the readme.txt file(s). 

3.) Once test is successful, add a full description of changes to MITgcm_contrib/ecco_darwin/tag-index. 
Include the date and list all specific changes.

Example: 

05/03/19
- To make initialization and pickups more robust:
  - remove pickupSuff from data
  - remove PTRACERS_Iter0 and PTRACERS_initialFile from data.ptracers
  - delete pickup_dic.0000000001.data so that it does not get used for
    initialization because of above changes

