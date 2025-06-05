# iter42/input is available at https://data.nas.nasa.gov/ecco/data.php?dir=/eccodata/llc_270/iter42/input
# ecco_darwin_v5/input is available at https://data.nas.nasa.gov/ecco/data.php?dir=/eccodata/llc_270/ecco_darwin_v5/input
# forcing/era_xx is available at https://ecco.jpl.nasa.gov/drive/files/Version5/Alpha/era_xx

# Instructions for building and running ECCO-Darwin v05 with Darwin 3 and the RADIv1 metamodel.
The RADIv1 metamodel was later added to the official Darwin3 pkg.

#This solution is documented in:
#Carroll, D., Menemenlis, D., Dutkiewicz, S., Lauderdale, J. M., Adkins, J. F., Bowman, K. W., et al. (2022). 
#Attribution of space-time variability in global-ocean dissolved inorganic carbon. Global Biogeochemical Cycles, 
#36, e2021GB007162. https://doi.org/10.1029/2021GB007162

#RADIv1 is documented in:
#Sulpis, O., Humphreys, M. P., Wilhelmus, M. M., Carroll, D., Berelson, W. M., Menemenlis, D., 
#Middelburg, J. J., and Adkins, J. F. (2022). RADIv1: a non-steady-state early diagenetic model for ocean sediments in Julia and MATLAB/GNU Octave, 
#Geoscientific Model Development, 15, 2105–2131. https://doi.org/10.5194/gmd-15-2105-2022, 2022.

==============
# 1. Get code

git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
git clone https://github.com/darwinproject/darwin3
cd darwin3
git checkout 24885b71
mkdir build run run_clim99
cd build

==============
# 2. Build executable

module purge
module load comp-intel mpi-hpe hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12
../tools/genmake2 -of ../../ecco_darwin/v05/llc270/code/linux_amd64_ifort+mpi_ice_nas \
  -mo '../../ecco_darwin/v05/llc270_RADIv1/code_darwin ../../ecco_darwin/v05/llc270/code_darwin ../../ecco_darwin/v05/llc270/code' -mpi
make depend
make -j 16

==============
# 3. Instructions for running simulation (1992-2020 period)

cd ../run
ln -sf ../build/mitgcmuv .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/nbp19_dmenemen_public_llc270/* .
ln -sf /nobackup/hzhang1/forcing/era_xx .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/input/19920101/to2021/xx*42.data .
cp ../../ecco_darwin/v05/llc270/input/* .
cp ../../ecco_darwin/v05/llc270_RADIv1/input/* .
mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
# modify job_ECCO_darwin as needed
qsub job_ECCO_darwin

==============
# 8. Instructions for running simulation (1985-2022 period w/ repetitive year 1999)

cd ../run_clim99
ln -sf ../build/mitgcmuv .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/nbp19_dmenemen_public_llc270/* .
ln -sf /nobackup/dcarrol2/forcing/apCO2/1984_NOAA_MBL/apCO2* .
ln -sf /nobackup/dcarrol2/pub/LLC_270/v05_1985_on/pickup_ptracers_1985_on.data \
        pickup_ptracers.0000000001.data
ln -sf /nobackup/hzhang1/forcing/exf_1999 .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/input/19850101/pickup* .
cp ../../ecco_darwin/v05/llc270_RADIv1/input/* .
cp ../../ecco_darwin/v05/llc270_RADIv1/input_1999/* .
mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
# modify job_ECCO_darwin as needed
qsub job_ECCO_darwin


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
        
        
