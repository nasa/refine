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
${two}/ref adapt annulus.meshb -g ${geomfile} -m annulus-metric.solb -x annulus1.meshb -t -f annulus1-final.tec
mv ref_gather_movie.tec annulus1_movie.tec
${two}/ref_acceptance -ugawg polar-2 annulus1.meshb annulus1-metric.solb
${two}/ref_metric_test annulus1.meshb annulus1-metric.solb > accept-annulus-uniform-01.status

${two}/ref adapt annulus1.meshb -g ${geomfile} -m annulus1-metric.solb -x annulus2.meshb -t -f annulus2-final.tec
mv ref_gather_movie.tec annulus2_movie.tec
${two}/ref_acceptance -ugawg polar-2 annulus2.meshb annulus2-metric.solb
${two}/ref_metric_test annulus2.meshb annulus2-metric.solb > accept-annulus-uniform-02.status

cat accept-annulus-uniform-02.status
../../check.rb accept-annulus-uniform-02.status 0.3 3.0



