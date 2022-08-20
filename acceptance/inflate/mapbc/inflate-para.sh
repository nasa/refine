#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

mpiexec -np 8 ${src}/refmpi collar normal \
    ../box.meshb \
    10 \
    0.01 \
    1.0 \
    1.68 \
    --usm3d-mapbc box.mapbc cigar 3

exit
