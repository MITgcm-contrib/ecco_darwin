# Build executable for ECCO-Darwin version 2 (ag4)
# corresponding to optimized solution in Brix et al. (2015)
# This reproduces executable with circa-2011 MITgcm code from
# /nobackup/hbrix/MITgcm_110502/build_darwin_ag4/mitgcmuv

cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co -D "2011-05-12 10:30:47" MITgcm_code
cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co -D "2011-05-12 10:30:47" MITgcm_contrib/darwin/pkg/darwin
cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co MITgcm_contrib/ecco_darwin/v2_cs510_Brix
cd MITgcm/pkg
ln -sf ../../MITgcm_contrib/darwin/pkg/darwin .
cd ..
mkdir build
cd build
module purge
module load comp-intel/2016.2.181 mpi-sgi/mpt.2.14r19 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
../tools/genmake2 -of \
    ../../MITgcm_contrib/ecco_darwin/v2_cs510_Brix/code/linux_amd64_ifort+mpi_ice_nas \
    -mo ../../MITgcm_contrib/ecco_darwin/v2_cs510_Brix/code
make depend
make -j 16

==============
# Instructions for running Version 2 (ag4) to 2009-2010.
# lou:~hbrix/ECCO2/darwin_ag4
#
# darwin_ag4:
# Starting beginning of 2009 from a weighted combination of the ptracers
# and dic pickup files from darwin runs 'a4v' 'a5v' 'a6v' 'a7v' 'a9v'
# with fact1 = 0.16033; fact2 = 0.56404; fact3 = -0.077899; fact4 = 0.10418;
# fact5 = 0.24935; Base run is a5v. Factors are from Green's
# function analysis of these four runs for 2009/10 using LDEO surface
# pCO2, Takahashi pCO2 climatology (variability, not absolute values and
# the global mean CO2 flux for 2010. pickup and pickup_seaice file are
# from 2009 adjoint integration.
#
# runs={'a4v' 'a5v' 'a6v' 'a7v' 'a9v'}; gfrun='ag4';
# fact1 = 0.16033; 
# fact2 = 0.56404; 
# fact3 = -0.077899; 
# fact4 = 0.10418;
# fact5 = 0.24935;
# piston velocity factor in dic_surfforcing.F from 0.337 to 0.34822
# R_PICPOC in darwin_generate_phyto.F from 0.04 to 0.133.
# optimization based on:
#     'a4v'    'a6v'    'a7v'    'a9v'    'a10v'    'a13v'
# Globflux_err = 0.01
# a5v -0.57061
# a4v 0.16033 +/- 0.010098
# a6v -0.077899 +/- 0.035022
# a7v 0.10418 +/- 0.0089553
# a9v 0.24935 +/- 0.0076539
# a10v -0.41537 +/- 0.047835
#  -> piston factor = 0.34822
# a13v 1.55 +/- 0.064771
#  -> R_PICPOC = 0.133
# global-mean flux: -2.5471

cd ..
mkdir run
cd run
ln -sf ../build/mitgcmuv .
cp ../../MITgcm_contrib/ecco_darwin/v2_cs510_Brix/input/* .
ln -sf ~dmenemen/CMS/run_template_cg1/darwin* .
ln -sf ~dmenemen/CMS/run_template_cg1/DIFFKR_2_20_1_lat6070_cube81 .
ln -sf /nobackup/hzhang1/cs510/run_template/GEBCO_510x6x510_ver06_dig.bin .
ln -sf /nobackup/hzhang1/forcing/jra25/jra25_dlw* .
ln -sf /nobackup/hzhang1/cs510/ICBC_2009_iter26/jra25_xx_* .
ln -sf /nobackup/hbrix/ICBC/pco2a_blended_* .
ln -sf pco2a_blended_2012 pco2a_blended_2013
ln -sf pco2a_blended_2012 pco2a_blended_2014
ln -sf /nobackupp8/dmenemen/CMS/lou_brix_ECCO2/darwin_ag4/pick* .
ln -sf /nobackup/dmenemen/forcing/runoff/runoff-360x180x12.bin .
ln -sf /nobackup/hzhang1/cs510/run_template/tile00* .

# modify run_darwin_450 as needed
qsub run_darwin_450

##############################################################
# From ftp://ecco.jpl.nasa.gov/ECCO2/Darwin/CarbFlux_v2
# Global Attributes:
# description             = 'ECCO2-Darwin CS510 run ag4 monthly average'
# MITgcmUV_version        = 'checkpoint62w'
# IntegrationPeriod       = '2009-01-01 -- 2010-12-31'
# InitialConditionsPhysics= 'ECCO2-CS510 ver. 2 (adjoint), Jan09-Apr10,
#                            iter. 22.' 
# InitialConditionsBGC    = 'Linear combination of the initial conditions of
#                            five sensitivity experiments for all
#                            biogeochemical tracers (Greens functions run ag4)'
# SurfaceForcingPhysics   = 'ECCO2-CS510 ver. 2 (adjoint), Jan09-Apr10,
#                            iter. 22. JRA-25 after Apr10.'
# SurfaceForcingBGC       = 'Time and space variable atmospheric pCO2 blended
#                            from daily 2009 GEOS-Chem results provided by K.
#                            Bowman (JPL) on 3/3/2011 and global atmospheric
#                            monthly means from NOAA ESRL. Iron/dust
#                            climatological values from Mahowald (2006).'
# ModelParametersPhysics  = 'ECCO2-CS510 ver. 2 (adjoint), Jan09-Apr10,
#                            iter. 22.'
# ModelParametersBGC      = '5 phytoplankton species, set-up from Stephanie
#                            Dutkiewicz (MIT);  PIC phytoplankton production
#                            relative to organic carbon (R_PICPOC) = 0.133.'
# author                  = 'Holger Brix -- UCLA/JPL'
# date                    = '10-Feb-2012'
# contact                 = 'hbrix@ucla.edu'

# From pfe:/home5/hbrix/README_runs
# darwin_ag4:
# Starting beginning of 2009 from a weighted comination of the ptracers
# and dic pickup files from darwin runs 'a4v' 'a5v' 'a6v' 'a7v' 'a9v'
# with fact1 = 0.16033; fact2 = 0.56404; fact3 = -0.077899; fact4 = 0.10418;
# fact5 = 0.24935; Base run is a5v. Factors are from Green's
# function analysis of these four runs for 2009/10 using LDEO surface
# pCO2, Takahashi pCO2 climatology (variability, not absolute values and
# the global mean CO2 flux for 2010. pickup and pickup_seaice file are
# from 2009 adjoint integration.

# From pfe:/home5/hbrix/README_gf
# runs={'a4v' 'a5v' 'a6v' 'a7v' 'a9v'}; gfrun='ag4';
# fact1 = 0.16033; 
# fact2 = 0.56404; 
# fact3 = -0.077899; 
# fact4 = 0.10418;
# fact5 = 0.24935;
# piston velocity factor in dic_surfforcing.F from 0.337 to 0.34822
# R_PICPOC in darwin_generate_phyto.F from 0.04 to 0.133.
# optimization based on:
#    'a4v'    'a6v'    'a7v'    'a9v'    'a10v'    'a13v'
# Globflux_err = 0.01
# a5v -0.57061
# a4v 0.16033 +/- 0.010098
# a6v -0.077899 +/- 0.035022
# a7v 0.10418 +/- 0.0089553
# a9v 0.24935 +/- 0.0076539
# a10v -0.41537 +/- 0.047835
# -> piston factor = 0.34822
# a13v 1.55 +/- 0.064771
# -> R_PICPOC = 0.133
# global-mean flux: -2.5471
