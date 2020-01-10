#!/usr/bin/env bash

set -x # echo commands

set -e # exit on first error
set -u # Treat unset variables as error

# Setup bash module environment
. /usr/local/pkgs/modules/init/bash

module purge
source acceptance/boot-modules.sh

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

cd ../C2S

cases=`find . -name accept.sh`
puts ${cases}
i=0
for case in ${cases}
    do
	dir=`dirname -- ${case}`
	(cd ${dir}  && ./accept.sh > accept-out.txt 2>&1 || touch $dir/FAILED ) &
	((i=i+1))
	if [ $((i%8)) -eq 0 ];
	then
	    sleep 30
	fi
    done
    
wait

find . -name FAILED

echo -e \\n\
# Build has failed if any failed cases have been reported
exit `find . -name FAILED | wc -l`

