# This code replaces the single-column swfrac with swfrac2d to add
# capability to read space-time varying Jerlov water type.
# At present, swfrac2d is not consistent with pkg/layers and pkg/seaice

# To use space-time-varying Jerlov water type:
# "#define ALLOW_WATERTYP" in CPP_OPTIONS.h and specify waterTypFile
# and other waterTyp* file description parameters in data.exf

# Running code with "#undef ALLOW_WATERTYP" in CPP_OPTIONS.h
# is identical to: (1) running v05,
# to (2) running with "waterTypFile = ' '," in data.exf, 
# and to (3) running with "waterTypFile = 'Jerlov_2',
# that is, they all default to Jerlov water type IA (jwtype=2).

# Also note https://github.com/MITgcm/MITgcm/pull/750
# which aims to implement swfrac3d, a 3d version of swfrac.

######################

Manfredi Manizza did PhD with Corrine Le Quere on this topic:

Modelling the impact of phytoplankton on upper ocean physics on global scale
M Manizza, C Le Quéré, AJ Watson, E Buitenhuis - 2003

Early steps towards modelling the phytoplankton-heat feedback
M Manizza, C Le Quéré, AJ Watson, E Buitenhuis… - EGS-AGU-EUG Joint Assembly,
2003

https://doi.org/10.1029/2004GL020778

Modelling phytoplankton-light feedback and its ocean biogeochemical implications
M Manizza - 2006 - University of East Anglia
