#!/usr/bin/env bash

set -e # exit on first error
set -u # treat unset variables as error
set -o pipefail # prevents errors in a pipeline from being masked

# Setup bash module environment
set +x # echo commands off for module
. /usr/local/pkgs/modules/init/bash
module purge
source acceptance/boot-modules.sh
module load cmake_3.15.5
module list
set -x # echo commands

log=`pwd`/../log-build.txt
trap "cat $log" EXIT
mkdir -p build
export CMAKE_PREFIX_PATH="${mpi_path}:${parmetis_path}:${egads_path}:${opencascade_path}"
( cd build && \
      cmake \
	  -DCMAKE_INSTALL_PREFIX=`pwd` \
      -DCMAKE_C_FLAGS="-g -O2" \
      -DCMAKE_C_COMPILER=gcc \
      .. > $log 2>&1\
      && make -j 8 install >> $log 2>&1\
      && ctest --output-on-failure >> $log 2>&1) || exit 1
trap - EXIT

log=`pwd`/../log-bootstrap.txt
trap "cat $log" EXIT
export PATH=${PATH}:`pwd`/build/bin
( cd acceptance/hemisphere/uniform && \
      ./generate.sh >> $log 2>&1 ) || exit 1
trap - EXIT

cat `pwd`/../log-build.txt
cat `pwd`/../log-bootstrap.txt
