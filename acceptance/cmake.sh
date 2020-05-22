#!/usr/bin/env bash

set -e # exit on first error
set -u # Treat unset variables as error

set +x # echo commands off for module

# Setup bash module environment
. /usr/local/pkgs/modules/init/bash
module purge
source acceptance/boot-modules.sh
module load cmake_3.15.5
module list

set -x # echo commands

log=`pwd`/../log-build.txt

export CMAKE_PREFIX_PATH=${mpi_path}:${egads_path}:${opencascade_path}
trap "cat $log" EXIT
mkdir -p build
( cd build && \
      cmake \
	  -DCMAKE_INSTALL_PREFIX=`pwd` \
      -DCMAKE_C_FLAGS="-g -O2" \
      -DCMAKE_C_COMPILER=gcc \
      .. > $log 2>&1\
      && make -j 8 install >> $log 2>&1) || exit 1
trap - EXIT
