#!/bin/bash
# Build serial 1-D darwin_ckpt68g executable on a compute node and set up run dirs.
# Called from PBS jobs; ROOT=/nobackup/dcarrol2/1-D
set -e
ROOT=/nobackup/dcarrol2/1-D
source /usr/share/modules/init/bash
module purge
module load comp-intel/2020.4.304 hdf4/4.2.12 hdf5/1.8.18_serial netcdf/4.4.1.1_serial python3/3.9.12
cd $ROOT/darwin3
# NAS FIPS OpenSSL blocks hashlib.md5() in darwin3's cog (see mitgcm-ecco build notes)
grep -q usedforsecurity tools/darwin/cogapp/cogapp.py || \
  sed -i 's/hashlib.md5()/hashlib.md5(usedforsecurity=False)/' tools/darwin/cogapp/cogapp.py
mkdir -p $ROOT/build && cd $ROOT/build
if [ ! -x mitgcmuv ]; then
  cp $ROOT/darwin3/tools/build_options/linux_amd64_ifort11 optfile
  NCI=$(nc-config --includedir); NFL=$(nf-config --flibs)
  printf '\nINCLUDES="$INCLUDES -I%s"\nLIBS="$LIBS %s"\n' "$NCI" "$NFL" >> optfile
  $ROOT/darwin3/tools/genmake2 -rootdir=$ROOT/darwin3 -mods=$ROOT/code -of=optfile > genmake.log 2>&1
  make depend > depend.log 2>&1
  make -j 8 > make.log 2>&1
fi
ls -la mitgcmuv
cd $ROOT
python3 scripts/make_relax.py sites/* > make_relax.log
