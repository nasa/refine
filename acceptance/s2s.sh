#!/usr/bin/env bash

set -e # exit on first error
set -u # Treat unset variables as error

# Setup bash module environment
set +x # echo commands off for module
source /usr/local/pkgs/modules/init/bash
module purge
source acceptance/s2s-modules.sh
module list
set -x # echo commands

log=`pwd`/../log-build.txt

./bootstrap > $log 2>&1
mkdir -p egads
( cd egads && \
      ../configure \
	  --prefix=`pwd` \
          --with-mpi=${mpi_path} \
          --with-parmetis=${parmetis_path} \
	  --with-EGADS=${egads_path} \
          --with-OpenCASCADE=${opencascade_path} \
	  CFLAGS="-g -O2" \
	  CC=icc >> $log 2>&1 \
      && make -j 8 >> $log 2>&1 \
      && make install >> $log 2>&1 \
    ) \
    || exit 1

export PATH=`pwd`/egads/src:${PATH}

cd ../tutorials

./continuous-integration/test-refine-cases-on-k.sh

