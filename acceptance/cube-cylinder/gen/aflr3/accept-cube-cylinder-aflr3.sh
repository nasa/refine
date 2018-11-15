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

geomfile=../csm/cube-cylinder.egads

ref_geom_test --tess ${geomfile} cube-cylinder.meshb 0.2 0.01 15.0


