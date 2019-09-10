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

geomfile=annulus.egads

# ${two}/ref_geom_test ${geomfile} annulus.meshb

${two}/ref_acceptance -ugawg polar-2 annulus.meshb annulus-metric.solb
${two}/ref_driver -i annulus.meshb -g ${geomfile} -m annulus-metric.solb -o ref_driver1 -t
mv ref_gather_movie.tec ref_driver1_movie.tec
${two}/ref_acceptance -ugawg polar-2 ref_driver1.meshb ref_driver1-metric.solb
${two}/ref_metric_test ref_driver1.meshb ref_driver1-metric.solb > accept-annulus-uniform-01.status

${two}/ref_driver -i ref_driver1.meshb -g ${geomfile} -m ref_driver1-metric.solb -o ref_driver2 -t
mv ref_gather_movie.tec ref_driver2_movie.tec
${two}/ref_acceptance -ugawg polar-2 ref_driver2.meshb ref_driver2-metric.solb
${two}/ref_metric_test ref_driver2.meshb ref_driver-metric.solb > accept-annulus-uniform-02.status

cat accept-annulus-uniform-02.status
../../check.rb accept-annulus-uniform-02.status 0.3 3.0



