# How to generate an ECCO regional cut out

*Authors:* **Clément Bertin[^1], Michael Wood[^2], Dustin Carroll[^1][^2]**

---
## General information
This repository has been created to guide ECCO/ECCO-Darwin users in generating their own regional configuration of the global ECCO state estimate.\
<u>Notes:</u> 
> These instructions require a good level of understanding on how ECCO/MITgcm model works and run.\
> The following instructions show how to extract a regional cut out from the ECCO-Darwin v5 estimate, but this can be reproduce on any ECCO product.

## Main steps
The instruction files are organized as follows:
1. **README.md**: Requirements for getting started.
2. **STEP1.md**: Instructions to generate the input files STEP2.
3. **STEP2.md**: Instructions to extract vectors from any llc global configuration along the boundaries of the required regional model using ``diagnostic_vec`` package (credit: Michael Wood)
4. **STEP3.md**: Instructions to setup the regional model

---
## Getting Started
To generate the regional configuration you will need:

1. Supercomputing capabilities to run ECCO state estimate (*Along this guide we give an example by running ECCO-Darwin v5 [llc270] model on Pleidaes NASA supercomputer*)
2. python3 with a dedicated anaconda environment (*This environment can be set up on the supercomputer or a local machine*). 

### 1. Get ECCO-Darwin v5 setup & merge ``diagnostic_vec`` with MITgcm

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

<u>Note:</u>: This package is not included in the official MITgcm realease but can be easly merged to it. More information on ``diagnostic_vec`` package at https://github.com/mhwood/diagnostics_vec (*Credit*: Mike Wood). 
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

> - Get ECCO-Darwin v5 setup
```
cd ../.. (back to downscaling folder)
git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
mkdir regions/configs/parent_run/
cd regions/configs/parent_run/
cp -r ecco_darwin/v05/llc270/code .
cp -r ecco_darwin/v05/llc270/code_darwin .
cp -r ecco_darwin/v05/llc270/input .
```

### 2. Create the python3 anaconda environment

You can either create this environment on your supercomputer or a local machine.

> - Install anaconda following the instructions on this page: https://conda.io/projects/conda/en/latest/user-guide/install/index.html
> - Open a terminal and Create the environment:
```
conda config --add channels defaults
conda config --add channels conda-forge
conda create --name downscaling python=3.12
conda activate downscaling
```
> - Install the following packages
```
conda install numpy matplotlib scipy pyproj netcdf4 xarray xesmf xgcm xmitgcm pyresample cartopy
```
> - Install [simplegrid](https://github.com/nasa/simplegrid) package, which is not available by `pip` or `conda install`. Instead, it must be cloned and then installed locally: 
```
git clone https://github.com/nasa/simplegrid.git
cd simplegrid
pip install .
```
You will need to update computegrid.py routine of simplegrid here:
```
cd [conda dir]/envs/downscaling/lib/python3.12/site-packages/simplegrid/
vim computegrid.py
:%s/np.PZERO/0.0/g
:wq
```
> - install the [ecco_v4_py](https://github.com/ECCO-GROUP/ECCOv4-py). While a `pip` and `conda`  install of `ecco-v4-py` is available (`pip install ecco-v4-py`), it is recommended to manually generate a package in the conda environment site-packages to avoid some issues:
```
git clone https://github.com/ECCO-GROUP/ECCOv4-py.git
mkdir [conda dir]/envs/downscaling/lib/python3.12/site-packages/ecco_v4_py
cp ECCOv4-py/ecco_v4_py/* [conda dir]/envs/downscaling/lib/python3.12/site-packages/ecco_v4_py/.
```
- Install [MITgcm utils](https://github.com/MITgcm/MITgcm):
```
git clone https://github.com/MITgcm/MITgcm.git
cd [MITgcm dir]/utils/python/MITgcmutils
python setup.py install
```
---
[^1]: Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA, USA

[^2]: Moss Landing Marine Laboratories, San José State University, Moss Landing, CA, USA
