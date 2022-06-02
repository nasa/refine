#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

geomfile=hemisphere.egads

if [ ! -f hemisphere-vol.meshb ]
then
    serveCSM -batch -skipTess ../hemisphere.csm
    ${src}/ref bootstrap ${geomfile}
fi

${src}/ref_acceptance hemisphere-vol.meshb hemisphere-vol-metric.solb 0.1
${src}/ref adapt hemisphere-vol.meshb \
      --egads ${geomfile} \
      --metric hemisphere-vol-metric.solb \
      --fun3d-mapbc hemisphere-vol.mapbc \
      --prism \
      -x hemihair1.meshb \
      -x hemihair1.plt \
      -x hemihair1.lb8.ugrid \
      -s 5

