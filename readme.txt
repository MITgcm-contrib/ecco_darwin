v2_cs510_Brix
Optimized solution described in Brix et al. (2015)

v3_cs510_Brix
Unpublished ECCO-Darwin v3 with linear piston velocity

v3_cs510_latest
Same as above with newer (checkpoint66o 2018/01/30) MITgcm code

v4_llc270
Optimized solution described in Carroll et al. (submitted 2019)

v5_llc270
Inclusion of LOAC (Land-Ocean-Aquatic-Coastal) forcing and
other improvements

**Team Development Process**

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

test from hzh@pfe
