# Point Dume downscaled regional set-up

---
## Preliminary information
This regional set-up is extracted from L1 East Pacific set up (3-4 km) extracted from global LLC270 ECCO-Darwin solution following the method detailed in https://github.com/MITgcm-contrib/ecco_darwin/tree/master/regions/downscaling. The global solution is documented in:
Carroll, D., Menemenlis, D., Dutkiewicz, S., Lauderdale, J. M., Adkins, J. F., Bowman, K. W., et al. (2022). Attribution of space-time variability in global-ocean dissolved inorganic carbon. Global Biogeochemical Cycles, 36, e2021GB007162. https://doi.org/10.1029/2021GB007162

---
## 1. Get code

```
git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
git clone https://github.com/darwinproject/darwin3
cd darwin3
git checkout 24885b71
mkdir build run
cd build
```

---
## 2. Build executable

```
module purge
module load comp-intel mpi-hpe/mpt hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12
../tools/genmake2 -of ../../ecco_darwin/v05/llc270/code/linux_amd64_ifort+mpi_ice_nas \
  -mo '../../ecco_darwin/regions/PointDume/code_darwin ../../ecco_darwin/regions/PointDume/code' -mpi
make depend
make -j 16
```

---
## 3. Instructions for running simulation

Copy forcings, initial and boundary conditions and namelists:
```
cd ../run
ln -sf ../build/mitgcmuv .

ln -sf /nobackup/rsavelli/LOAC/PointDume/forcings/exf_L2 .
ln -sf /nobackup/rsavelli/LOAC/PointDume/forcings/GloFas* .
ln -sf /nobackup/rsavelli/PointDume_highres/grid/forcings/pickups/*0002419488* .
ln -sf /nobackup/rsavelli/PointDume_highres/grid/forcings/OBCS/* .
ln -sf /nobackup/rsavelli/PointDume_highres/grid/PointDume_20m_bathymetry_r8.bin .
ln -sf /nobackup/rsavelli/PointDume_highres/grid/delYFile .
ln -sf /nobackup/rsavelli/LOAC/PointDume/forcings/L2_wave .
cp ../../ecco_darwin/regions/PointDume/input/* .
cp ../../ecco_darwin/regions/PointDume/input_darwin/* .
```

Create diagnostic directory:
```
mkdir diags
mkdir diags/EtaN_day_snap
mkdir diags/EtaN_mon_mean
mkdir diags/SI_daily_snap
mkdir diags/TS_surf_daily_snap
mkdir diags/vel_surf_daily_snap
mkdir diags/TS_AW_daily_snap
mkdir diags/state_3D_mon_snap
mkdir diags/state_3D_mon_mean
mkdir diags/vel_3D_mon_snap
mkdir diags/vel_3D_mon_mean
mkdir diags/BGC_daily_consts
mkdir diags/BGC_daily_DO
mkdir diags/BGC_daily_PO
mkdir diags/BGC_daily_misc
mkdir diags/BGC_daily_cx
mkdir diags/BGC_daily_Chl
```

Submit batch file for running the model:
```
qsub job_pointdume
```
