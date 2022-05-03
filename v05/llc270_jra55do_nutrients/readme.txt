# iter42/input is available at https://data.nas.nasa.gov/ecco/data.php?dir=/eccodata/llc_270/iter42/input
# ecco_darwin_v5/input is available at https://data.nas.nasa.gov/ecco/data.php?dir=/eccodata/llc_270/ecco_darwin_v5/input
# forcing/era_xx is available at https://ecco.jpl.nasa.gov/drive/files/Version5/Alpha/era_xx
# LOACv1.4.0_HJ available at https://ecco.jpl.nasa.gov/drive/files/ECCO2/LLC270/LOAC/LOACv1.4.0_HJ
# LOACriver_temp available at https://ecco.jpl.nasa.gov/drive/files/ECCO2/LLC270/LOAC/LOACriver_temp

# Instructions for building and running ECCO-Darwin v05 with Darwin 3, JRA55-do, and globalNEWS nutrient runoff

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
../tools/genmake2 -of ../../ecco_darwin/v05/llc270/code/linux_amd64_ifort+mpi_ice_nas \
  -mo '../../ecco_darwin/v05/llc270_jra55do_nutrients/code_darwin ../../ecco_darwin/v05/llc270_jra55do/code'
make depend
make -j 16

==============
# 3. Instructions for running ECCO-Darwin v05 with Darwin 3, JRA55-do, and globalNEWS nutrient runoff for 1992-2018 period

cd ../run
ln -sf ../build/mitgcmuv .
ln -sf /nobackupp2/dmenemen/public/llc_270/iter42/input/* .
ln -sf /nobackupp2/dmenemen/public/llc_270/ecco_darwin_v5/input/darwin_initial_conditions/* .
ln -sf /nobackupp2/dmenemen/public/llc_270/ecco_darwin_v5/input/darwin_forcing/* .
ln -sf /nobackup/hzhang1/forcing/era_xx .
ln -sf /nobackup/hzhang1/forcing/jra55_do/LOACv1.4.0_HJ .
ln -sf /nobackup/hzhang1/forcing/jra55_do/LOACriver_temp .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/input/19920101/to2021/xx*42.data .
ln -sf /nobackup/dcarrol2/LOAC/write_bin/globalNEWS/* .
cp ../../ecco_darwin/v05/llc270_jra55do/input/* .
cp ../../ecco_darwin/v05/llc270_jra55do_nutrients/input/data.darwin .
mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
# modify job_llc270_fdH as needed
qsub job_llc270_fdH
