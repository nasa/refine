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

geomfile=cylinder.egads

${two}/ref_geom_test --tetgen ${geomfile} cylinder-vol.meshb 1 1 90 \
    | tee accept-revolve-cylinder-gen-vol.status

