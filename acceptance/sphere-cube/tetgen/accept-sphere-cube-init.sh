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

${two}/ref_geom_test --tetgen ../geom/sphere-cube.egads sphere-cube-init.meshb \
      0.1 0.01 15 \
      --surf | tee sphere-cube-init.out

${two}/ref_driver -i sphere-cube-init.meshb -g ../geom/sphere-cube.egads \
      -x sphere-cube-crv.meshb -r 2 -t | tee sphere-cube-crv.out

