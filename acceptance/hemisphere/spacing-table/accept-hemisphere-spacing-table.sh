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

serveCSM -batch -skipTess ../hemisphere.csm
${src}/ref bootstrap ${geomfile}

${src}/ref_acceptance hemisphere-vol.meshb hemisphere-vol-metric.solb 0.1
${src}/ref adapt hemisphere-vol.meshb \
      --egads ${geomfile} \
      --fun3d-mapbc hemisphere-vol.mapbc \
      --spacing-table spacing-table.txt \
      -x hemitab1.meshb \
      -x hemitab1.plt \
      -s 10

