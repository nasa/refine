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

geomfile=../geom/cube-wire.egads
cp ../bootstrap/cube-wire-vol.meshb cube-wire01.meshb

${two}/ref_driver -i cube-wire01.meshb -g ${geomfile} -o cube-wire02 -r 1 -t \
      -f cube-wire02-final.tec 
mv ref_gather_movie.tec cube-wire01-movie.tec

${two}/ref_driver -i cube-wire02.meshb -g ${geomfile} -o cube-wire03 -r 10 -t \
      -f cube-wire03-final.tec 
mv ref_gather_movie.tec cube-wire02-movie.tec

${two}/ref_histogram_test cube-wire03.meshb cube-wire03-final-metric.solb \
 > accept-cube-wire-curvature.status

cat accept-cube-wire-curvature.status
../../check.rb accept-cube-wire-curvature.status 0.2 2.2



