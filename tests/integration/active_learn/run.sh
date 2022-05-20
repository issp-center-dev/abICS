#!/bin/sh

export OMP_NUM_THREADS=1

set -e

rm -f ALloop.progress
rm -rf AL0
rm -rf AL1
rm -rf MC0
rm -rf MC1

sh ./AL.sh
sh ./MC.sh
sh ./AL.sh
sh ./MC.sh

if [ -e MC1/kTs.npy ] ; then
  echo OK
  exit 0
else
  echo FAILED
  exit 1
fi
