# This code replaces the single-column swfrac with swfrac2d to add
# capability to read space-time varying Jerlov water type.
# At present, swfrac2d is not consistent with pkg/layers and pkg/seaice
# See ../v05_kpp for a version that is compatible with pkg/kpp

# Running code with "#undef ALLOW_WATERTYP" in CPP_OPTIONS.h
# is identical to: (1) running v05,
# to (2) running with "waterTypFile = ' '," in data.exf,
# and to (3) running with "waterTypFile = 'Jerlov_2',
# that is, they all default to Jerlov water type IA (jwtype=2).

# Also note https://github.com/MITgcm/MITgcm/pull/750
# which aims to implement swfrac3d, a 3d version of swfrac.
