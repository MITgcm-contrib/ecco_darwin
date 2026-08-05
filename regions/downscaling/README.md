# How to generate an ECCO regional cut out

*Authors:* **Clément Bertin[^1], Michael Wood[^2], Dustin Carroll[^1][^2]**

---
## General information
This repository has been created to guide ECCO/ECCO-Darwin users in generating their own regional configuration of the global ECCO state estimate.\
<u>Notes:</u> 
> These instructions require a good level of understanding on how ECCO/MITgcm model works and run.\
> The following instructions show how to extract a regional cutout from the ECCO-Darwin v5 estimate, but this can be reproduced on any ECCO product.

## Main steps
The instruction files are organized as follows:
1. **README.md**: Requirements for getting started — code, python environment, and the directory layout used throughout.
2. **STEP1.md**: Instructions to generate the regional grid, bathymetry, tiles, netcdf grid file, and the ``diagnostics_vec`` mask files needed by STEP2.
3. **STEP2.md**: Instructions to extract vectors from any llc global configuration along the boundaries of the required regional model using the ``diagnostics_vec`` package (credit: Michael Wood).
4. **STEP3.md**: Instructions to generate the initial and boundary conditions of the regional model, and to verify them.
5. **STEP4.md**: Instructions to build, configure, and run the downscaled regional set-up.
6. **STEP5.md**: Instructions to verify the downscaled run with a short experiment — that it integrated stably, and that it actually used the initial and boundary conditions generated in STEP3.

---
## Directory layout

Two separate directory trees are used throughout these instructions. Nearly every command below
is written relative to one of them, so it is worth setting them up as shown.

### 1. The working tree (on the supercomputer)

This is created in "Getting Started" below. Everything the model compiles and runs from lives here:

```
downscaling/                          <-- the top-level folder you create
├── darwin3/                          <-- MITgcm/darwin3 source (+ merged diagnostics_vec)
│   ├── pkg/diagnostics_vec/
│   ├── tools/
│   └── regions/                      <-- build/run directories are created here
│       ├── gen_ncgrid/{build,run}/       (STEP1 IV)
│       ├── parent_run/{build,run}/       (STEP1 V, STEP2 II)
│       └── downscaled/{build,run}/       (STEP4 IV)
├── diagnostics_vec/                  <-- Mike Wood's package, merged into darwin3
├── ecco_darwin/                      <-- this repository
└── regions/
    ├── <reg_nm>.mitgrid              <-- copied up from your grid directory
    ├── <reg_nm>_bathymetry.bin
    ├── tiles/
    └── configs/
        ├── parent_run/{code,code_darwin,input}/
        └── gen_ncgrid/{code,input,reg_example_files}/ + job_downsc_ivy
```

<u>Note:</u> the build and run directories under ``darwin3/regions/`` sit four levels below
``downscaling/``, which is why the ``-mo``/``-of`` paths in later steps are written as
``../../../../regions/...`` (four levels up, into ``downscaling/regions/``) while ``genmake2``
and ``tools/`` are reached with ``../../../`` (three levels up, into ``darwin3/``).

### 2. The configuration directory (`config_dir`, where the python utils work)

Every python script in ``utils/`` takes a ``-d``/``--config_dir`` argument. They all assume the
same layout underneath it, and several create their own subfolders. Point ``-d`` at the **same
directory every time**:

```
<config_dir>/                         <-- the -d argument of every utils script
├── <reg_nm>.mitgrid                  <-- created by gen_mitgrid.py       (STEP1 I)
├── <reg_nm>_bathymetry.bin           <-- created by gen_bathy.py         (STEP1 II)
├── <reg_nm>_ncgrid.nc                <-- created by stitch_ncgrid.py     (STEP1 IV.c)
├── delYFile                          <-- created by gen_delYFile.py      (STEP4 II.a)
├── tiles/                            <-- created by mitgrid2tiles.py     (STEP1 III)
├── mncs/                             <-- you copy the mnc_* folders here (STEP1 IV.c)
├── parent/
│   ├── inputs/                       <-- created by gen_dvmasks.py: <bnd>_BC_mask.bin
│   └── outputs/
│       ├── grid/                     <-- you copy the parent grid files here (STEP1 V)
│       ├── pickups/                  <-- you copy the parent pickup files here (STEP3 I)
│       └── OBCS/                     <-- you copy the diagnostics_vec output here (STEP3 I)
├── forcings/
│   ├── pickups/                      <-- created by gen_pickups.py       (STEP3 I)
│   ├── OBCS/                         <-- created by gen_obcs.py          (STEP3 II)
│   └── OBCS_corrected/               <-- created by check_obcs_transport.py -correct (STEP3 III)
└── diagnostics/                      <-- plots from the STEP3 III and STEP5 checks
```

<u>Note:</u> the folder is ``parent``, **singular**. The scripts build these paths literally
(e.g. ``gen_obcs.py`` reads ``parent/outputs/grid/``), so a ``parents/`` or ``parent/output/``
directory will not be found.

<u>Note:</u> ``config_dir`` can be anywhere convenient on the machine running the anaconda
environment. In the examples throughout these instructions it is written as
``ecco_darwin/regions/downscaling/`` for brevity, but a dedicated directory outside the
repository is usually cleaner.

---
## Getting Started
To generate the regional configuration, you will need:

1. Supercomputing capabilities to run ECCO state estimate (*Along this guide we give an example by running ECCO-Darwin v5 [llc270] model on the Pleiades NASA supercomputer*)
2. python3 with a dedicated anaconda environment (*This environment can be set up on the supercomputer or a local machine*). 

### 1. Get ECCO-Darwin v5 set-up & merge ``diagnostics_vec`` with MITgcm

Below, we detail the instructions to run ECCO-Darwin on the NASA Pleiades supercomputer.\
This section builds the working tree described in "Directory layout" above.

> - clone darwin checkpointv67x github
```
mkdir downscaling
cd downscaling
git clone https://github.com/darwinproject/darwin3
cd darwin3
git checkout 24885b71
```
> - clone diagnostics_vec github

<u>Note:</u> This package is not included in the official MITgcm release but can be easily merged to it. More information on the ``diagnostics_vec`` package at https://github.com/mhwood/diagnostics_vec (*Credit*: Mike Wood). 
```
cd ..
git clone https://github.com/mhwood/diagnostics_vec.git
```
> - Merge diagnostics_vec package into darwin3
```
cd diagnostics_vec/utils/
python3 copy_doc_files_to_MITgcm.py -m ../../darwin3/
python3 copy_pkg_files_to_MITgcm.py -m ../../darwin3/
python3 copy_verification_files_to_MITgcm.py -m ../../darwin3/
```

> - Get ECCO-Darwin v5 set-up

<u>Note:</u> the ``cp`` commands below are run from ``downscaling/``, so the destination is given
as a path rather than ``.``. Running them from inside ``regions/configs/parent_run/`` will fail,
because ``ecco_darwin/`` is not visible from there.
```
cd ../..
git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
mkdir -p regions/configs/parent_run/
cp -r ecco_darwin/v05/llc270/code regions/configs/parent_run/
cp -r ecco_darwin/v05/llc270/code_darwin regions/configs/parent_run/
cp -r ecco_darwin/v05/llc270/input regions/configs/parent_run/
```

### 2. Create the python3 anaconda environment

You can either create this environment on your supercomputer or a local machine.

> - Install anaconda following the instructions on this page: https://conda.io/projects/conda/en/latest/user-guide/install/index.html
> - Open a terminal and create the environment:
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
> - Install [MITgcm utils](https://github.com/MITgcm/MITgcm) (``MITgcmutils``, used by most of the ``utils`` scripts to read model binaries):
```
git clone https://github.com/MITgcm/MITgcm.git
cd [MITgcm dir]/utils/python/MITgcmutils
pip install .
```

<u>Note:</u> ``python setup.py install`` also works on older setuptools, but is deprecated; ``pip install .`` is the supported equivalent.

---
[^1]: Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA, USA

[^2]: Moss Landing Marine Laboratories, San José State University, Moss Landing, CA, USA
