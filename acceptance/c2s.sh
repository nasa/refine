#!/usr/bin/env bash

set -e # exit on first error
set -u # Treat unset variables as error

# Setup bash module environment
. /usr/local/pkgs/modules/init/bash

source acceptance/c2s-modules.sh

set -x # echo

./bootstrap
mkdir -p egads
( cd egads && \
    ../configure \
    --prefix=`pwd` \
    --with-EGADS="/u/shared/fun3d/fun3d_users/modules/ESP/114/EngSketchPad" \
    CFLAGS="-g -O2" \
    && make -j \
    && make install \
    ) \
    || exit 1

export PATH=${PATH}:`pwd`/egads/src

cd ../acceptance

dir=`pwd`/om6/geom
log=`pwd`/../log-om6-geom-init-grid.txt
(cd $dir && ./init-grid.sh  > $log 2>&1 || touch $dir/FAILED ) &

wait

find \$c2s -name FAILED

echo -e \\n\
# Build has failed if any failed cases have been reported
exit `find \$c2s -name FAILED | wc -l`

