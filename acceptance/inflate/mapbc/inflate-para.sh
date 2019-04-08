#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    one=$1/one
    two=$1/src
else
    one=${HOME}/refine/parmetis/one
    two=${HOME}/refine/parmetis/src
fi

mpiexec -np 8 ${two}/ref_inflatable  \
    ../box.meshb \
    10 \
    0.01 \
    1.0 \
    1.68 \
    --mapbc box.mapbc cigar 3

exit
