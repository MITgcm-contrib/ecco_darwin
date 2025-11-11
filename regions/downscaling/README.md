# How to generate an ECCO-Darwin regional cut-out

Authors: Clément Bertin, Michael Wood, Dustin Carroll

## General information
This repository has been created to guide ECCO-Darwin users in generating their own regional configuration of the global-ocean ECCO-Darwin simulation. 
By following the instructions, you will be able to extract the information needed at the boundaries of your regional cut-out from ECCO-Darwin and generate the initial/boundary conditions with lon-lat coordinates.

## Main steps
The instruction files are organized as follows:
1. **README.md**: Requirements for getting started.
2. **STEP1.md**: Instructions to generate regional model grid files.
3. **STEP2.md**: Instructions to extract vectors from any LLC global configuration along the boundaries of the regional model using the ``diagnostic_vec`` package (credit: Michael Wood).
4. **STEP3.md**: Instructions to generate regional model set-up.

## Getting Started
To generate the regional configuration you will need:

1. Supercomputing capabilities to run the global-ocean v05 ECCO-Darwin model (*Along this guide we give an example by running v05 ECCO-Darwin (LLC 270) model on Pleidaes NASA supercomputer*).
2. Python3 with a dedicated anaconda environment (*This environment can be set up on the supercomputer or a local machine*). 

### 1. Get v05 ECCO-Darwin set-up & merge ``diagnostic_vec`` with MITgcm

Below, we detail the instructions to run ECCO-Darwin on the NASA Pleiades supercomputer:

> - clone darwin checkpointv67x github
```
mkdir downscalling 
cd downscalling
git clone https://github.com/darwinproject/darwin3
cd darwin3
git checkout 24885b71
```
> - clone diagnostic_vec github

**Note**: This package is not included in the official MITgcm realease but can be easly added to it. More information on ``diagnostic_vec`` package at https://github.com/mhwood/diagnostics_vec (*Credit*: Mike Wood). 
```
cd ..
git clone https://github.com/mhwood/diagnostics_vec.git
```
> - Merge diagnostic_vec package to darwin3
```
cd diagnostics_vec/utils/
python3 copy_doc_files_to_MITgcm.py -m ../../darwin3/
python3 copy_pkg_files_to_MITgcm.py -m ../../darwin3/
python3 copy_verification_files_to_MITgcm.py -m ../../darwin3/
```

> - Get v05 ECCO-Darwin set-up
```
cd ../.. (back to downscaling folder)
git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
mkdir config
cp -r ecco_darwin/v05/llc270/code config/.
cp -r ecco_darwin/v05/llc270/code_darwin config/.
cp -r ecco_darwin/v05/llc270/input config/.
```

### 2. Create the python3 anaconda environment

You can either create this environment on your supercomputer or a local machine.

> - Install anaconda following the instructions on this page: https://conda.io/projects/conda/en/latest/user-guide/install/index.html
> - Open a terminal and Create the environment:
```
conda create --name downscaling
conda activate downscaling
```
> - Install the following packages
```
conda install numpy matplotlib scipy pyproj netcdf4 xarray
conda install -c conda-forge xesmf 
```
> - Install [simplegrid](https://github.com/nasa/simplegrid) package, which is not available by `pip` or `conda install`. Instead, it must be cloned and then installed locally: 
```
git clone https://github.com/nasa/simplegrid.git
cd simplegrid
pip install .
```
> - install the [ecco_v4_py](https://github.com/ECCO-GROUP/ECCOv4-py). While a `pip` and `conda`  install of `ecco-v4-py` is available (`pip install ecco-v4-py`), it is recommended to manually generate a package in the conda environment site-packages to avoid some issues:
```
git clone https://github.com/ECCO-GROUP/ECCOv4-py.git
mkdir [conda dir]/envs/downscaling/lib/python3.12/site-packages/ecco_v4_py
cp ECCOv4-py/ecco_v4_py/* [conda dir]/envs/downscaling/lib/python3.12/site-packages/ecco_v4_py/.
```
> - Install the following additional packages
```
conda install -c conda-forge xgcm xmitgcm pyresample cartopy
```
- Install [MITgcm utils](https://github.com/MITgcm/MITgcm):
```
git clone https://github.com/MITgcm/MITgcm.git
cd [MITgcm dir]/utils/python/MITgcmutils
python setup.py install
```
