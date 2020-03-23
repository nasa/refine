#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    one=$1/one
    two=$1/src
else
    one=${HOME}/refine/egads/one
    two=${HOME}/refine/egads/src
fi

${two}/ref_inflatable  \
    ../box.meshb \
    -10 \
    0.005 \
    0.10 \
    1.68 \
    4 5 6 \
    --shift 2 0 5

exit
