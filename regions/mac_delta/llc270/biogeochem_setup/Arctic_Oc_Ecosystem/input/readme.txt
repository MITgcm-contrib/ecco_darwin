From Steph (4/20/2021)

Changes to:

·        data.pkg  - to include radtrans

·        data.ptracers – additional tracers CDOM (ptracer20) and bacteria (ptracer 28)

NOTE: these will need initial files, and that ptracer20 onward from your last set up to match new

·        data.darwin – I used my original file and added in you DARWIN_FORCING_PARAMS and DARWIN_INTERP_PARAMS.

Search for CSD for couple other comments

See the DARWIN_RADTRANS_PARAMS – you’ll need new pathways (DIRNAME) for the three optics files (also attached)
·        data.traits: completely new contents

·        data.radtrans: new file – need to add in pathways and filenames for OASIM input, and your own periodicity (currently just month climatology – would be better to use actual monthly if that is what Oliver has)​