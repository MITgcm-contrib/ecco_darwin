A set of 1-day ("endtime=86400," in data) verification experiments were
carried out to check the swfrac2d code with pkg/kpp:

STDOUT.0000_kpp
is from running ecco_darwin/regions/RedSea/v05_kpp
with "git checkout 7d24893", that is after replacement of
pkg/ggl90 with pkg/kpp but before inclusion of swfrac2d code.

STDOUT.0000_undef_ALLOW_WATERTYP
is from running ecco_darwin/regions/RedSea/v05_kpp
but with "#undef ALLOW_WATERTYP" in CPP_OPTIONS.h

STDOUT.0000_Jerlov_2
is from running ecco_darwin/regions/RedSea/v05_kpp
with "waterTypFile = 'Jerlov_2'," in data.exf
where Jerlov_2 is a 12-month climatology filled
with constant value 2.

STDOUT.0000_Jerlov_Kdclim12
is from running ecco_darwin/regions/RedSea/v05_kpp
with "waterTypFile = 'Jerlov_Kdclim12'," in data.exf
where Jerlov_2 is a 12-month climatology of Jerlov
water type for the RedSea region.

Since the first three verification experiments default to
Jerlov water type IA (jwtype=2), they are bit identical.
Only the last experiment gives different results.

JerlovKdclim12.jpg
shows January-16 Jerlov water type from Jerlov_Kdclim12

SSTdiff.jpg
displays the SST difference in mdeg C for the Jerlov_Kdclim12
simulation minus a simulation with Jerlov=2 after a 1-day
model integration starting on January 16, 1992.
