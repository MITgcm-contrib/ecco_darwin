# ECCO-Darwin: Mackenzie Delta configuration

You will find below the instruction to run ECCO-Darwin Mackenzie Delta simulations

**WARNING: These instructions works only for Pleïades users**

## 1. Get the code
```
git clone https://github.com/darwinproject/darwin3
git clone https://github.com/MITgcm-contrib/ecco_darwin
```

## 2. Build executable

```
cd darwin3
mkdir build run
cd build
module purge
module load comp-intel mpi-sgi hdf4 hdf5 netcdf
mv  ../../ecco_darwin/PATH/TO/SETUP/FILES/code/packages.conf ../../ecco_darwin/PATH/TO/SETUP/FILES/code/packages.conf_org
ln -s ../../ecco_darwin/PATH/TO/SETUP/FILES/code_darwin/packages.conf ../../ecco_darwin/PATH/TO/SETUP/FILES/code
../tools/genmake2 -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas -mpi -mo '../../ecco_darwin/PATH/TO/SETUP/FILES/code_darwin ../../ecco_darwin/PATH/TO/SETUP/FILES/code'
make depend
make -j 16
```
## 3. Run the setup
```
cd ../run
mkdir diags diags/daily diags/budget
ln -sf ../build/mitgcmuv .
```
### Make links to forcing files
```
ln -sf /nobackup/hzhang1/forcing/era_xx .
-- Freswater runnoff --
ln -sf /nobackup/cbertin/Forcing/river_runoff/Freswater/AGRO_interan .
ln -sf /nobackup/cbertin/Forcing/river_runoff/Temperature/Tokuda_Mac270modif .
-- biogeochemical runoff --
ln -sf /nobackup/cbertin/Forcing/river_runoff/Nutrients/Bertin_etal_21/Interannual/L20_R80/tDOCl .
ln -sf /nobackup/cbertin/Forcing/river_runoff/Nutrients/Bertin_etal_21/Interannual/L20_R80/tDOCr .
ln -sf /nobackup/cbertin/Forcing/river_runoff/Nutrients/Tank_etal_12/Interannual/tAlk .
ln -sf /nobackup/cbertin/Forcing/river_runoff/Nutrients/Tank_etal_12/Interannual/tDIC .
ln -sf /nobackup/cbertin/Forcing/river_runoff/Nutrients/GNW2_NutCim/tDON .
ln -sf /nobackup/cbertin/Forcing/river_runoff/Nutrients/GNW2_NutCim/tDOP .
ln -sf /nobackup/cbertin/Forcing/river_runoff/Nutrients/GNW2_NutCim/tDSi .
-- run files -- 
ln -sf /nobackup/cbertin/Forcing/run_template/PATH/* .
```
### Copy setup data files
```
cp ../../ecco_darwin/PATH/TO/SETUP/FILES/input/* .
cp ../../ecco_darwin/Pregions/mac_delta/llc270/job_mac270_bro .
```
### Run
```
qsub job_mac270_bro
```
