#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

mpiexec -np 8 ${src}/refmpi collar radial  \
    ../box.meshb \
    10 \
    0.005 \
    0.10 \
    1.68 \
    --fun3d-mapbc inflated.mapbc

exit
