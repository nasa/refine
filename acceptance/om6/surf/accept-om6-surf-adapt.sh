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

${two}/ref_geom_test --tess ../recon/onera-m6-sharp-te.egads om6-init.meshb \
      80 0.3 15 \
      --surf

${two}/ref_driver -i om6-init.meshb -g ../recon/onera-m6-sharp-te.egads \
      -x om6-crv.meshb -r 2 -t

