#!/bin/sh

# This script installs MLIP-3

# NOTE for macOS users:
# 1. macOS's `make` is GNU make 3.x, which is not supported by MLIP-3.
#    You need to install GNU make 4.x by `brew install make` and use `gmake` instead of `make` as
#    $ MAKE=gmake sh ./install_mlip3.sh
# 2. `gcc` (/usr/bin/gcc) does not mean GNU CC and hence the compilation will fail.
# To use this script, make sure that `gcc` command invokes GNU CC by, for example,
#   $ ln -s `which gcc-11` ./gcc
#   $ ln -s `which g++-11` ./g++
#   $ PATH=`pwd`:$PATH sh ./install_mlip3.sh

if [ -z ${MAKE} ]; then
  MAKE=make
fi

set -ue

URL=https://gitlab.com/ashapeev/mlip-3/-/archive/main/mlip-3-main.tar.gz
if [ ! -e mlip-3.tar.gz ]; then
  wget $URL -O mlip-3.tar.gz
fi
rm -rf mlip-3
mkdir mlip-3
tar zxf mlip-3.tar.gz -C mlip-3 --strip-components=1
cd mlip-3

./configure --no-mpi --compiler=gnu

cd make

GFORTRAN_VERSION=$(gfortran -dumpversion | cut -d. -f1)
if [ $GFORTRAN_VERSION -ge 10 ]; then
  echo "FFLAGS += -fallow-argument-mismatch" >> config.mk
fi

cd ..

$MAKE mlp