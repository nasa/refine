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

# ${two}/ref_geom_test ${geomfile} cone-cone.meshb

${two}/ref_acceptance cone-cone.meshb cone-cone.metric 0.1
${two}/ref_driver -i cone-cone.meshb -g ${geomfile} -m cone-cone.metric -o ref_driver1 -t -r 1
mv ref_gather_movie.tec ref_driver1_movie.tec
${two}/ref_acceptance ref_driver1.meshb ref_driver1.metric 0.1
${two}/ref_metric_test ref_driver1.meshb ref_driver1-final-metric.solb > accept-cone-cone-uniform-two-01.status

${two}/ref_driver -i ref_driver1.meshb -g ${geomfile} -m ref_driver1.metric -o ref_driver2 -t -r 1
mv ref_gather_movie.tec ref_driver2_movie.tec
${two}/ref_acceptance ref_driver2.meshb ref_driver2.metric 0.1
${two}/ref_metric_test ref_driver2.meshb ref_driver2-final-metric.solb > accept-cone-cone-uniform-two-02.status

cat accept-cone-cone-uniform-two-02.status
../../../check.rb accept-cone-cone-uniform-two-02.status 0.20 6.0

