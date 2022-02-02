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

geomfile=ega.egads

# ${two}/ref_geom_test ${geomfile} ega.meshb

${two}/ref_acceptance ega.meshb ega.metric 0.1
${two}/ref adapt ega.meshb -g ${geomfile} -m ega.metric -o ega1.meshb -t
mv ref_gather_movie.tec ega1_movie.tec
${two}/ref_acceptance ega1.meshb ega1.metric 0.1
${two}/ref_metric_test ega1.meshb ega1.metric > accept-cube-cylinder-uniform-01.status

${two}/ref adapt ega1.meshb -g ${geomfile} -m ega1.metric -o ega2.meshb -t
mv ref_gather_movie.tec ega2_movie.tec
${two}/ref_acceptance ega2.meshb ega2.metric 0.1
${two}/ref_metric_test ega2.meshb ega2.metric > accept-cube-cylinder-uniform-02.status

cat accept-cube-cylinder-uniform-02.status
../../check.rb accept-cube-cylinder-uniform-02.status 0.3 3.0



