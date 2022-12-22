# iter42/input is available at https://data.nas.nasa.gov/ecco/data.php?dir=/eccodata/llc_270/iter42/input
# ecco_darwin_v4/input is available at https://data.nas.nasa.gov/ecco/data.php?dir=/eccodata/llc_270/ecco_darwin_v4/input
# forcing/era_xx is available at https://ecco.jpl.nasa.gov/drive/files/Version5/Alpha/era_xx

# Instructions for building and running ECCO-Darwin v4 with Darwin 3

==============
# 1. Get code

git clone --depth 1 https://github.com/darwinproject/darwin3
git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
cd darwin3
mkdir build run
cd build

==============
# 2. Build executable for ECCO-Darwin v4 with Darwin 3

module purge
module load comp-intel/2020.4.304 mpi-hpe/mpt.2.25 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
../tools/genmake2 -of ../../ecco_darwin/v05/llc270_ship_track/code/linux_amd64_ifort+mpi_ice_nas \
  -mo '../../ecco_darwin/v05/llc270_ship_track/code_darwin ../../ecco_darwin/v05/llc270_ship_track/code'
make depend
make -j 16

==============
# 2. Instructions for running ECCO-Darwin v4 with Darwin 3 for 1992-2018 period

cd ../run
ln -sf ../build/mitgcmuv .
ln -sf /nobackupp19/dmenemen/public/llc_270/iter42/input/* .
ln -sf /nobackupp19/dmenemen/public/llc_270/ecco_darwin_v5/input/darwin_initial_conditions/* .
ln -sf /nobackupp19/dmenemen/public/llc_270/ecco_darwin_v5/input/darwin_forcing/* .
ln -sf /nobackup/hzhang1/forcing/era_xx .
cp ../../ecco_darwin/v05/llc270/input/* .
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