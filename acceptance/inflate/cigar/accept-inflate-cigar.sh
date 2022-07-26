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
${src}/ref bootstrap ${geomfile}

${src}/ref adapt cigar-vol.meshb \
      --egads ${geomfile} \
      -x cigar.meshb \
      -x cigar.plt

