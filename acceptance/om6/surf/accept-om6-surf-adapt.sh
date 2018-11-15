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

${two}/ref_geom_test --tetgen ../recon/onera-m6-sharp-te.egads om6-surf.meshb \
      80 0.3 15 \
      --surf

