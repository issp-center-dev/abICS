#!/bin/sh

# This script installs aenet by using gfortran
# generate.x_serial, train.x_mpi, and predict.x_serial will be copied
# into ~/opt/aenet/bin
#
# NOTE:
# For macOS, `gcc` (/usr/bin/gcc) does not mean GNU CC and hence the compilation will fail.
# To use this script, make sure that `gcc` command invokes GNU CC by, for example,
#   ln -s `which gcc-11` ./gcc
#   PATH=`pwd`:$PATH sh ./install_aenet.sh

set -ue

VERSION=2.0.4

URL=https://github.com/atomisticnet/aenet/archive/refs/tags/v${VERSION}.tar.gz
wget $URL -O aenet.tar.gz
rm -rf aenet
mkdir aenet
tar zxf aenet.tar.gz -C aenet --strip-components=1

cd aenet/lib
make
cd ../src

GFORTRAN_VERSION=$(gfortran -dumpversion | cut -d. -f1)
if [ $GFORTRAN_VERSION -ge 10 ]; then
  LOCAL_FCFLAGS="-fallow-argument-mismatch -O2 -fexternal-blas \$(DEBUG)"
else
  LOCAL_FCFLAGS="-O2 -fexternal-blas \$(DEBUG)"
fi

cd makefiles
sed -i.orig "s/FCFLAGS *=.*/FCFLAGS = ${LOCAL_FCFLAGS}/" Makefile.gfortran_serial
sed -i.orig "s/FCFLAGS *=.*/FCFLAGS = ${LOCAL_FCFLAGS} -DPARALLEL/" Makefile.gfortran_mpi
cd ..

make -f makefiles/Makefile.gfortran_serial

cd ../bin
mkdir -p ~/opt/aenet/bin
cp generate.x-${VERSION}-gfortran_serial ~/opt/aenet/bin/generate.x_serial
cp predict.x-${VERSION}-gfortran_serial ~/opt/aenet/bin/predict.x_serial

cd ../src
make clean
make -f makefiles/Makefile.gfortran_mpi
cd ../bin
cp train.x-${VERSION}-gfortran_mpi ~/opt/aenet/bin/train.x_mpi
