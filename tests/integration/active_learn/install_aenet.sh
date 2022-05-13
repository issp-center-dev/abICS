#!/bin/sh

# This script installs aenet by using gfortran
# generate.x_serial, train.x_mpi, and predict.x_serial will be copied
# into ~/opt/aenet/bin

VERSION=2.0.4

URL=https://github.com/atomisticnet/aenet/archive/refs/tags/v${VERSION}.tar.gz
wget $URL -O aenet.tar.gz
mkdir aenet
tar zxf aenet.tar.gz -C aenet --strip-components=1

cd aenet/lib
make
cd ../src

GFORTRAN_VERSION=$(gfortran -dump-version | cut -d. -f1)
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
make -f makefiles/Makefile.gfortran_mpi

cd ../bin
mkdir -p ~/opt/aenet/bin
cp generate.x-${VERSION}-gfortran_serial ~/opt/aenet/bin/generate.x_serial
cp train.x-${VERSION}-gfortran_mpi ~/opt/aenet/bin/train.x_mpi
cp predict.x-${VERSION}-gfortran_serial ~/opt/aenet/bin/predict.x_serial
