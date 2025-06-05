# Experimental set-up for mangrove carbon export based on LLC270 jra55-do physical solution.
# Exports of C are estimated according to method in Savelli et al (in prep)
# Exports of C are an hourly climatology based on year 2000.
# EXF package has been modified to account for useYearlyFields and RepCycle at the same time 
# for repeating each year hourly C runoff over 1992-2020 period.
# To do that, data.cal has been change to 'NoLeapYear'.
#
# iter42/input is available at https://data.nas.nasa.gov/ecco/data.php?dir=/eccodata/llc_270/iter42/input
# ecco_darwin_v5/input is available at https://data.nas.nasa.gov/ecco/data.php?dir=/eccodata/llc_270/ecco_darwin_v5/input
# forcing/era_xx is available at https://ecco.jpl.nasa.gov/drive/files/Version5/Alpha/era_xx
# LOACv1.4.0_HJ available at https://ecco.jpl.nasa.gov/drive/files/ECCO2/LLC270/LOAC/LOACv1.4.0_HJ
# LOACriver_temp available at https://ecco.jpl.nasa.gov/drive/files/ECCO2/LLC270/LOAC/LOACriver_temp

#For details on nutrient loads, please refer to ecco_darwin/code_util/LOAC/GlobalNews/readme.txt

# Instructions for building and running ECCO-Darwin v05 with Darwin 3, JRA55-do, and mangrove carbon runoff

==============
# 1. Get code

git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
git clone https://github.com/darwinproject/darwin3
cd darwin3
git checkout 24885b71
mkdir build run
cd build

==============
# 2. Build executable

module purge
module load comp-intel mpi-hpe hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
../tools/genmake2 -of ../../ecco_darwin/v05/llc270/code/linux_amd64_ifort+mpi_ice_nas \
  -mo '../../ecco_darwin/v05/llc270_jra55do_mangroves/code_darwin ../../ecco_darwin/v05/llc270/code_darwin ../../ecco_darwin/v05/llc270_jra55do_mangroves/code ../../ecco_darwin/v05/llc270_jra55do/code ../../ecco_darwin/v05/llc270/code' -mpi
make depend
make -j 16

==============
# 3. Instructions for running ECCO-Darwin v05 with Darwin 3, JRA55-do, and globalNEWS nutrient runoff for 1992-2018 period

cd ../run
ln -sf ../build/mitgcmuv .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/nbp19_dmenemen_public_llc270/* .
ln -sf /nobackup/hzhang1/forcing/era_xx .
ln -sf /nobackup/hzhang1/forcing/jra55_do/LOACv1.4.0_HJ .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/input/19920101/to2021/xx*42.data .
ln -sf /nobackup/rsavelli/Model_mangrove_export/outputs/LLC_270_0MSL/bgc_runoff/* .
cp ../../ecco_darwin/v05/llc270/input/* .
cp ../../ecco_darwin/v05/llc270_jra55do/input/* .
cp ../../ecco_darwin/v05/llc270_jra55do_nutrients/input/* .
mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
# modify job_llc270_fdH as needed
qsub job_llc270_fdH

==============
# Additional important notes

---------------
--data.seaice--
---------------

Lines in data.seaice:
 # revert to old defaults before 2018-08-25
  SEAICEscaleSurfStress = .FALSE.,
  SEAICEaddSnowMass     = .FALSE.,
  SEAICE_useMultDimSnow = .FALSE.,
  SEAICEetaZmethod = 0,
  SEAICE_Olx       = 0,
  SEAICE_Oly       = 0,
  SEAICE_drag      = .002,
  SEAICE_waterDrag = .00534499514091350826,

are needed when switching from:
        cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co -D "11/28/17" MITgcm_code
to:
        git clone git://gud.mit.edu/darwin3-dev darwin3
        
because of following change to MITgcm:

checkpoint67d (2018/09/04)
o pkg/seaice: some new and hopefully more sensible defaults:
  - SEAICEscaleSurfStress = .TRUE.
  - SEAICEaddSnowMass = .TRUE.
  - SEAICEadvScheme = 77
  - SEAICE_useMultDimSnow = .TRUE.
  - SEAICE_Olx/y = Olx/y - 2
  - SEAICEetaZmethod = 3
  - SEAICE_drag = 0.001
  - CPP-flag SEAICE_ZETA_SMOOTHREG defined by default
  - change SEAICE_waterDrag default to 5.5e-3 and add reference density
    rhoConst to scale it. Add a warning if SEAICE_waterDrag > 1.
    This change also involves a slight re-organization (= simplification)
    in seaice_oceandrag_coeffs.F
    
----------------
