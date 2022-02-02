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

geomfile=cube-cylinder.egads

# ${two}/ref_geom_test ${geomfile} cube-cylinder.meshb.meshb

${two}/ref_acceptance cube-cylinder.meshb cube-cylinder.metric 0.1
${two}/ref adapt cube-cylinder.meshb -g ${geomfile} -m cube-cylinder.metric -x cube-cylinder1.meshb
${two}/ref_acceptance cube-cylinder1.meshb cube-cylinder1.metric 0.1
${two}/ref_metric_test cube-cylinder1.meshb cube-cylinder1.metric > accept-cube-cylinder-uniform-01.status

${two}/ref adapt cube-cylinder1.meshb -g ${geomfile} -m cube-cylinder1.metric -x cube-cylinder2.meshb
${two}/ref_acceptance cube-cylinder2.meshb cube-cylinder2.metric 0.1
${two}/ref_metric_test cube-cylinder2.meshb cube-cylinder2.metric > accept-cube-cylinder-uniform-02.status

cat accept-cube-cylinder-uniform-02.status
../../check.rb accept-cube-cylinder-uniform-02.status 0.16 1.6



