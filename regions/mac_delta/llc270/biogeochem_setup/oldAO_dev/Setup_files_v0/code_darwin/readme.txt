From Steph (4/20/2021):

Changes made to

·        DARWIN_OPTIONS.h   - made several changes, including adding in ptracer CDOM (will need initial condition for this).

your code looks bit old – has “GUD” references in it still, so you might want to check

 

I added two options that were not in your code:

--DARWIN_ION_SED_SOURCE_POP (undef)

-- DARWIN_DIAG_PERTYPE (define)

And removed one which was in your code, but not mine: DARWIN_DIAG_TENDENCIES

But left your DARWIN_NUTRIENT_RUNOFF

NOTE: did not include DARWIN_ALLOW_DENIT – seems not important in AO? But also no nitrogen fixers, so might as well make nitrogen conversed. But will need to change if shifting to global.

 

·        DARWIN_SIZE.h : changes for new number plankton, groups and optical types

Note that nPPplank and nGRplank ae obsolete in newest code

 

·        PTRACERS_SIZE.h: add two new ptracers (CDOM and bacteria)

 

·        RADTRANS_SIZE.h: new file for this mod

 

·        Packages_conf: added radtrans and sun ​