#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

geomfile=cigar.egads

serveCSM -batch -skipTess cigar.csm
${src}/ref bootstrap ${geomfile} -s 5

${src}/ref adapt cigar-vol.meshb \
      -s 8 \
      --egads ${geomfile} \
      --uniform cyl ceil 0.05 1  0 0 0  1 0 0  1 1 \
      -x cigar.meshb \
      -x cigar.plt

${src}/ref collar interpolated  \
    cigar.meshb \
    10 \
    0.1 \
    5.0 \
    1.4 \
    --usm3d-mapbc cigar-usm3d.mapbc cigar 3 \
    --origin 0 0 0.5 \
    -x cigar-inflated.meshb \
    -x cigar-inflated.plt

