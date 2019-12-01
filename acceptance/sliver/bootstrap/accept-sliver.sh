#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

serveCSM -batch sliver.csm
${src}/ref boostrap sliver.egads

geomfile=sliver.egads

${src}/ref_driver -i sliver-vol.meshb -g ${geomfile} -o sliver01 -t \
      -f sliver01-final.tec 
mv ref_gather_movie.tec sliver01-movie.tec

${src}/ref_histogram_test sliver01.meshb sliver01-final-metric.solb \
 > accept-sliver-curvature.status

cat accept-sliver-curvature.status
../../check.rb accept-sliver-curvature.status 0.01 2.0



