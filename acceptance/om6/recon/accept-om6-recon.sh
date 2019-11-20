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

geomfile=onera-m6-sharp-te.egads

${two}/ref_egads_test --recon om6-curv-no-geom.meshb ${geomfile} > accept-om6-recon.status

../../recon.rb accept-om6-recon.status 1.0e-8
