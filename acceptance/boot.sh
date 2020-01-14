#!/usr/bin/env bash

set -e # exit on first error
set -u # Treat unset variables as error

set +x # echo commands off for module

# Setup bash module environment
. /usr/local/pkgs/modules/init/bash
module purge
source acceptance/boot-modules.sh
module list

set -x # echo commands

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

export PATH=`pwd`/egads/src:${PATH}

output=`pwd`/../log-bootstrap-status.txt

cd ../C2S

sketches=`find . -name accept.sh`
echo ${sketches}
i=0
for sketch in ${sketches}
do
    dir=`dirname -- ${sketch}`
    (cd ${dir}  && ./accept.sh > accept-out.txt 2>&1 || touch ./FAILED ) &
    ((i=i+1))
    if [ $((i%12)) -eq 0 ];
    then
	sleep 30
    fi
done
    
wait

statuses=`find . -name accept-bootstrap-status.txt | sort`
for status in ${statuses}
do
    dir=`dirname -- ${status}`
    echo "`cat ${status}` ${dir}" >> ${output}
done

cat ${output}
find . -name FAILED

echo -e \\n\
# Build has failed if any failed sketches have been reported
exit `find . -name FAILED | wc -l`

