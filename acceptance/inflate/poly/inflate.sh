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

serveCSM -batch poly.csm 

${two}/ref bootstrap poly.egads

${two}/ref_inflatable  \
    poly-vol.meshb \
    10 \
    0.1 \
    2.0 \
    1.68 \
    2 3 4 5 6 7 8 9

exit
