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

geomfile=cone-cone.egads

${two}/ref_egads_test --recon cone-cone-no-geom.meshb ${geomfile} > accept-cone-cone-recon.status

../../recon.rb accept-cone-cone-recon.status
