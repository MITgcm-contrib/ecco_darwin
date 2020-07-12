# iter42/input is available at https://data.nas.nasa.gov/ecco/data.php?dir=/eccodata/llc_270/iter42/input
# forcing/era_xx is available at https://ecco.jpl.nasa.gov/drive/files/Version5/Alpha/era_xx
# LOACv1.4.0_HJ available at https://ecco.jpl.nasa.gov/drive/files/ECCO2/LLC270/LOAC/LOACv1.4.0_HJ
# LOACriver_temp available at https://ecco.jpl.nasa.gov/drive/files/ECCO2/LLC270/LOAC/LOACriver_temp


==============
# Build executable for forward-only llc270 iteration 42 optimized solution
git clone https://github.com/MITgcm-contrib/ecco_darwin
git clone https://github.com/MITgcm/MITgcm.git
cd MITgcm
git checkout `git rev-list -n 1 --first-parent --before="2017-11-28 00:00" master`

mkdir build
cd build
module purge
module load comp-intel/2016.2.181 mpi-sgi/mpt.2.14r19 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
MOD="../../ecco_darwin/v05/llc270_jra55do"
../tools/genmake2 -of ${MOD}/code/linux_amd64_ifort+mpi_ice_nas \
	          -mo ${MOD}/code
make depend
make -j 16


==============
# Instructions for running forward-only llc270 optimized solution (1992-2015)
cd ..
mkdir run
cd run
mkdir diags
ln -sf ../build/mitgcmuv .
ln -sf /nobackupp2/dmenemen/public/llc_270/iter42/input/* .
ln -sf /nobackup/hzhang1/forcing/era_xx .
ln -sf /nobackup/hzhang1/forcing/jra55_do/LOACv1.4.0_HJ .
ln -sf /nobackup/hzhang1/forcing/jra55_do/LOACriver_temp .
cp ${MOD}/input/* .
# modify job_llc270_fdH as needed
qsub job_llc270_fdH
