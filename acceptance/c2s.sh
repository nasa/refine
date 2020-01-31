#!/usr/bin/env bash

set -x # echo commands

set -e # exit on first error
set -u # Treat unset variables as error

# Setup bash module environment
. /usr/local/pkgs/modules/init/bash

module purge
source acceptance/c2s-modules.sh

log=`pwd`/../log-build.txt

./bootstrap > $log 2>&1
mkdir -p egads
( cd egads && \
      ../configure \
	  --prefix=`pwd` \
	  --with-EGADS=${egads_path} \
	  CFLAGS="-g -O2" \
	  CC=gcc >> $log 2>&1 \
      && make >> $log 2>&1 \
      && make install >> $log 2>&1 \
    ) \
    || exit 1

mkdir -p parmetis
( cd parmetis && \
    ../configure \
	--prefix=`pwd` \
	--with-metis=${metis_path} \
	--with-parmetis=${parmetis_path} \
	--with-EGADS=${egads_path} \
	--enable-lite \
	CFLAGS="-DHAVE_MPI -g -O2" \
	CC=mpicc >> $log 2>&1 \
      && make >> $log 2>&1 \
      && make install >> $log 2>&1 \
    ) \
    || exit 1

export PATH=`pwd`/egads/src:${PATH}

cd ../acceptance

dir=`pwd`/om6
log=`pwd`/../log-om6-run.txt
(cd $dir && ./run.sh  > $log 2>&1 || touch $dir/FAILED ) &

dir=`pwd`/jsm-nac/geom
log=`pwd`/../log-jsm-nac-geom-init-grid.txt
(cd $dir && ./nac-box.sh  > $log 2>&1 || touch $dir/FAILED ) &

wait

find . -name FAILED

echo -e \\n\
# Build has failed if any failed cases have been reported
exit `find . -name FAILED | wc -l`

