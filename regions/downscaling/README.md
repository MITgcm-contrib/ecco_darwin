# Generating ECCO-Darwin regional cut out configuration

Authors: Clément Bertin, Michael Wood, Dustin Carroll

## General information
This repository has been created to guide the ECCO-Darwin users in generating their own regional configuration of the global ECCO-Darwin simulation. By following the instructions you will be able to extract the information needed at the boundaries of your regional cut out from any llc global simulation and generate the boundary conditions with latitude/longitude coordinates.

## Main steps
The instruction files are organised as follows:
1. **README.md**: You will find here the basic requirements to getting started.
2. **STEP1.md**: you will find here the instructions to exctract vectors from any llc global configuration along the boundaries of the required regional model using ``diagnostic_vec`` package (credit Michael Wood: https://github.com/mhwood/diagnostics_vec) 
3. **STEP2.md**: you will find here the instructions to convert the exctracted vector into the boundary/intial conditions of your regional setup using the python3 codes provided.

## Getting Started
To generate the regional configuration you will need:
1. Capabilities to run any llc global model (Here we give an example by running ECCO-Darwin v5 [llc270] model on Pleidaes NASA superconputer)
2. python3 with a dedicated anaconda environment. Instructions to set the environment follow.

> - Install anaconda following the instructions on this page: https://conda.io/projects/conda/en/latest/user-guide/install/index.html
> - Open a terminal and Create the environment:
```
conda create --name downscaling
conda activate downscaling
```
> - Install the following packages
```
conda install numpy matplotlib scipy pyproj netcdf4
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
